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
    def __init__(self, file):
        self.file = file
    
    def writeList(self, l):
        logger.info(f'   Write list to file {self.file}.')
        with open(self.file, 'w') as fh:
            fh.writelines(f'{i}\n' for i in l)
    
    def readList(self):
        l = []
        with open(self.file, 'r') as fh:
            logger.info('   Reading from previous.')
            [l.append(i.strip()) for i in fh]
        return l


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
    try:
        allTax = txids_io.readList()
    except FileNotFoundError:
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
    eachGroup_ass = 200 # submit this number of ids in each request. Prevent dead response.
    try:
        allAssemblies = asmids_io.readList()
    except FileNotFoundError:
        if len(allTax) > eachGroup_ass:
            taxGroups = []
            nGroups = len(allTax)//eachGroup_ass + (1 if len(allTax)%eachGroup_ass>0 else 0)
            logger.info(f'   {nGroups} groups.')
            for i in range(nGroups):
                try:
                    taxGroups.append(allTax[i*eachGroup_ass:i*eachGroup_ass+eachGroup_ass])
                except IndexError:
                    taxGroups.append(allTax[i*eachGroup_ass:])
        allAssemblies = []
        for i, taxs in enumerate(taxGroups):
            logger.info(f'   Getting assembly ids for group {i+1}/{nGroups}')
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
            handle = Entrez.elink(
                db='assembly',
                dbfrom='taxonomy',
                id=taxs,
                linkname='taxonomy_assembly',
            )
            record_a = Entrez.read(handle)
            handle.close()
            [allAssemblies.extend([l['Id'] for l in r['LinkSetDb'][0]['Link']]) for r in record_a if len(r['LinkSetDb']) > 0]
            # not all taxid have corresponding assembly, many have more than 1 assembly.
        asmids_io.writeList(allAssemblies)
    logger.info(f'   Total number of assemblies is {len(allAssemblies)}.')
    logger.info(f'   The first 5 are: {allAssemblies[:min(len(allAssemblies), 5)]}')

    # 3nd step
    # Get nucleic id from refseq
    logger.info(f'3. Fetch nucl ids of fetched assemblies, only refseq.')
    eachGroup_nucl = 200 # submit this number of ids in each request. Prevent dead response.
    try:
        allNucl = nuclids_io.readList()
    except FileNotFoundError:
        if len(allAssemblies) > eachGroup_nucl:
            assGroups = []
            nGroups = len(allAssemblies)//eachGroup_nucl + (1 if len(allAssemblies)%eachGroup_nucl>0 else 0)
            logger.info(f'   {nGroups} groups.')
            for i in range(nGroups):
                try:
                    assGroups.append(allAssemblies[i*eachGroup_nucl:i*eachGroup_nucl+eachGroup_nucl])
                except IndexError:
                    assGroups.append(allAssemblies[i*eachGroup_nucl:])
        allNucl = []
        for i, asses in enumerate(assGroups):
            logger.info(f'   Getting nucleotide ids for group {i+1}/{nGroups}')
            handle = Entrez.elink(
                db = 'nuccore',
                dbfrom = 'assembly',
                id=asses,
                linkname='assembly_nuccore_refseq', # will remove most irrelevent sequences
                term=f'{minLength}:99999999[SLEN]'
            )
            record_n = Entrez.read(handle)
            handle.close()
            [allNucl.extend([l['Id'] for l in r['LinkSetDb'][0]['Link']]) for r in record_n if len(r['LinkSetDb']) > 0]
        nuclids_io.writeList(allNucl)
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

    logger = logging.getLogger(__name__)

    fetch_nucl_by_taxID(targetTx='201174', minLength=10000, targetDir=targetDir, api=args.api, email=args.email)
    # Actinobacteria, some super tax
