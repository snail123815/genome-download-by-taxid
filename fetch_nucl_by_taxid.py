from email.errors import InvalidMultipartContentTransferEncodingDefect
from Bio import Entrez
import os
import pickle
import lzma
from datetime import date
from time import strftime, sleep
from urllib.error import HTTPError
from http.client import IncompleteRead
import sys
from Bio import SeqIO
from Bio import Seq
import subprocess
import logging


from sklearn.utils import all_estimators
from calHash_on_args import calHash
import pickle


class RecordIO:
    def __init__(self, file, prefix='IDs from source id groups '):
        self.file = file
        self.pf = prefix.strip() + ' '
        self._completeMark = "ALL"
        self.logger = logging.getLogger()
        try:
            with open(self.file, 'r') as fh:
                headerLine = fh.readline()
        except FileNotFoundError:
            headerLine = ""
        self.gn, self.nGroups, self.completed = self._readHeader(headerLine)
        if self.nGroups == -1 and not self.completed:
            with open(self.file, 'w') as fh:
                fh.write(self.pf + f'{self.gn}/{self.nGroups}' + '\n')
                # done/total
    
    def _readHeader(self, headerLine):
        pf_from_file = ' '.join(headerLine.split(' ')[:-1])
        gn = 0
        nGroups = -1
        completed = False
        if pf_from_file == self.pf[:-1]:
            mark = headerLine.split(' ')[-1].strip()
            if mark == self._completeMark:
                gn = -1
                nGroups = -1
                completed = True
            else:
                gn, nGroups = mark.split('/')
                gn = int(gn)
                nGroups = int(nGroups)
                assert gn < nGroups or nGroups == -1, f'File content error, {gn} should < {nGroups}'
                completed = False
        elif pf_from_file == '':
            pass
        else:
            self.logger.warning(f'File prefix changed, file will be overwritten.')
        return gn, nGroups, completed
    
    def writeList(self, l):
        self.logger.info(f'   Write whole list to file {self.file}.')
        with open(self.file, 'w') as fh:
            fh.write(self.pf + self._completeMark + "\n")
            fh.writelines(f'{i}\n' for i in l)
        self.completed = True

    def readList(self):
        l = []
        with open(self.file, 'r') as fh:
            headerLine = fh.readline()
            self.gn, nGroups, self.completed = self._readHeader(headerLine)
            if self.nGroups == -1:
                self.nGroups = nGroups
            else:
                assert self.nGroups == nGroups
            self.logger.info(f'   Reading from previous file {headerLine.strip()}.')
            [l.append(i.strip()) for i in fh]
        return l

    def appendList(self, l, groupIdx=None):
        groupIdx = groupIdx or self.gn
        assert groupIdx == self.gn, f'Index not correct: should be {self.gn}, but {groupIdx}'
        with open(self.file, 'r') as fh:
            headerLine = fh.readline()
        mark = headerLine.split(' ')[-1].strip()
        if mark == self._completeMark:
            self.logger.error(f'   The file is already full! Want to start again?')
            sys.exit()
        else:
            gn, nGroups = mark.split('/')
            gn = int(gn)
            nGroups = int(nGroups)
            if nGroups == -1 and self.nGroups == -1 and gn == 0: # only warn once
                self.logger.warning('   Total group number not known. Please do RecordIO.nGroups = n before appending list.')
            elif self.nGroups == -1:
                self.nGroups = nGroups
            elif nGroups == -1:
                nGroups = self.nGroups
            assert nGroups == self.nGroups, 'Target file and target list are not from the same source.'
            assert gn == self.gn, 'Missing group or previous "sed" not successful'
            self.gn += 1
            if self.gn == self.nGroups:
                self.completed = True
            newHeader = self.pf + (f'{gn+1}/{nGroups}' if not self.completed else self._completeMark)
            # change the first line of target file
            cmd = ['sed', '-i', '', '-e', '1,1s/.*/' + newHeader.replace('/', '\/') + '/g', self.file]
            subprocess.call(cmd)
        with open(self.file, 'a') as fh:
            fh.writelines(f'{i}\n' for i in l)


def splitGroup(listObj, step):
    ll = [] # list of list with maximum {step} of items each
    if len(listObj) > step:
        nGroups = len(listObj)//step + (1 if len(listObj)%step>0 else 0)
        logger.info(f'   {nGroups} groups.')
        for i in range(nGroups):
            try:
                ll.append(listObj[i*step:i*step+step])
            except IndexError:
                ll.append(listObj[i*step:])
    else:
        ll.append(listObj)
    return ll


def fetch_nucl_by_taxID(targetTx, minLength, targetDir, api, email):
    Entrez.api_key = api
    Entrez.email = email
    txHash = calHash(targetTx)
    nuclIdHash = calHash(targetTx, minLength)

    txids_file = os.path.join(targetDir, f'txid_{txHash}.txt')
    txids_io = RecordIO(txids_file)
    asmids_file = os.path.join(targetDir, f'asmb_{txHash}.txt')
    asmids_io = RecordIO(asmids_file)
    nuclids_file = os.path.join(targetDir, f'nucl_{nuclIdHash}.txt')
    nuclids_io = RecordIO(nuclids_file)


    # TODO what if you want one single tax? 0 return for the first elink? or else?

    logger.info(f'{" Start Fetch Nucleotides by Taxonomy ID ":=^80}')
    # 1st step
    # Get all child taxonomy from database
    logger.info(f'1. Fetch childs of taxid {targetTx}.')
    allTax = txids_io.readList()
    if not txids_io.completed:
        handle = Entrez.elink(
            db='taxonomy',
            dbfrom='taxonomy',
            id=targetTx,
            linkname='taxonomy_taxonomy_child',
        )
        record_t = Entrez.read(handle)
        handle.close()

        allTax = [x['Id'] for x in record_t[0]['LinkSetDb'][0]['Link']]
        # Only valid if there is only one record, if expecting more records, use the method from 2nd step
        txids_io.writeList(allTax)
    logger.info(f'   Total number of child taxids is {len(allTax)}.')
    logger.info(f'   The first 5 are: {allTax[:5]}')

    # 2nd step
    # Get assemblies from taxid
    logger.info(f'2. Fetch linked assembly ids of fetched taxids')
    allTax = allTax[:1400]
    eachGroup_ass = 500 # submit this number of ids in each request. Prevent dead response.
    taxGroups = splitGroup(allTax, eachGroup_ass)
    allAssemblies = asmids_io.readList()
    if asmids_io.nGroups == -1:
        asmids_io.nGroups = len(taxGroups)
    i = asmids_io.gn
    while not asmids_io.completed:
        assemblies = []
        logger.info(f'   Getting assembly ids for group {i+1}/{asmids_io.nGroups}')
        handle = Entrez.elink(
            db='assembly',
            dbfrom='taxonomy',
            id=taxGroups[i],
            linkname='taxonomy_assembly',
        )
        record_a = Entrez.read(handle)
        handle.close()
        [assemblies.extend([l['Id'] for l in r['LinkSetDb'][0]['Link']]) for r in record_a if len(r['LinkSetDb']) > 0]
        asmids_io.appendList(assemblies, i)
        allAssemblies.extend(assemblies)
        i += 1
        # not all taxid have corresponding assembly, many have more than 1 assembly.
    logger.info(f'   Total number of assemblies is {len(allAssemblies)}.')
    logger.info(f'   The first 5 are: {allAssemblies[:min(len(allAssemblies), 5)]}')

    # 3nd step
    # Get nucleic id from refseq
    logger.info(f'3. Fetch nucl ids of fetched assemblies, only refseq.')
    eachGroup_nucl = 200 # submit this number of ids in each request. Prevent dead response.
    assGroups = splitGroup(allAssemblies, eachGroup_nucl)
    allNucl = nuclids_io.readList()
    if nuclids_io.nGroups == -1:
        nuclids_io.nGroups = len(assGroups)
    i = nuclids_io.gn
    while not nuclids_io.completed:
        nucls = []
        logger.info(f'   Getting nucleotide ids for group {i+1}/{nuclids_io.nGroups}')
        handle = Entrez.elink(
            db = 'nuccore',
            dbfrom = 'assembly',
            id=assGroups[i],
            linkname='assembly_nuccore_refseq', # will remove most irrelevent sequences
            term=f'{minLength}:99999999[SLEN]'
        )
        record_n = Entrez.read(handle)
        handle.close()
        [nucls.extend([l['Id'] for l in r['LinkSetDb'][0]['Link']]) for r in record_n if len(r['LinkSetDb']) > 0]
        nuclids_io.appendList(nucls)
        allNucl.extend(nucls)
        i += 1
    logger.info(f'   Total number of nucleotides is {len(allNucl)}.')
    logger.info(f'   The first 5 are: {allNucl[:min(len(allNucl), 5)]}')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('api', type=str, help='api key from NCBI')
    parser.add_argument('email', type=str, help='Email address for identify yourself')

    args = parser.parse_args()

    targetDir = './test_output'
    os.makedirs(targetDir, exist_ok=True)
    logging.basicConfig(filename=os.path.join(targetDir, 'Fetch_taxonomy_nucl.log'), level=logging.INFO)

    logger = logging.getLogger()

    fetch_nucl_by_taxID(targetTx='201174', minLength=10000, targetDir=targetDir, api=args.api, email=args.email)
    # Actinobacteria, some super tax

        ########################## do not work for elink
        # count = 0
        # start = 0
        # webenv = None
        # retmax = 2000
        # while start < count or webenv is None:
        #     handle = Entrez.elink(
        #         db='assembly',
        #         dbfrom='taxonomy',
        #         id=taxs,
        #         usehistory=True,
        #         webenv=webenv,
        #         linkname='taxonomy_assembly',
        #         retmax=retmax
        #     )
        #     record_a = Entrez.read(handle)
        #     handle.close()
        #     count = int(record_a['Count'])
        #     if count == '0':
        #         logger.warning(f'   Zero record in this group')
        #         break
        #     if 'WarningList' in record_a:
        #         logger.warning(record_a['WarningList']['OutputMessage'])
        #     if start == 0:
        #         webenv = record_a['WebEnv']
        #         logger.info(f'   Web accession {webenv}')
        #         logger.info(f'   Total ids {count}.')
        #     start += retmax
        #     [allAssemblies.extend([l['Id'] for l in r['LinkSetDb'][0]['Link']]) for r in record_a if len(r['LinkSetDb']) > 0]
        ##########################