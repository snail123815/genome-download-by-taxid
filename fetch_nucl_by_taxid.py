from Bio import Entrez
import os
import pickle
from time import strftime, sleep
import sys
import subprocess
import logging
import pickle
from modules.calHash_on_args import calHash
import threading
from concurrent.futures import ThreadPoolExecutor


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
            self.logger.info(f'   Reading from previous file, header: {headerLine.strip()}.')
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
            with open(self.file, 'r') as fh:
                lines = fh.readlines()
            lines[0] = newHeader + '\n'
            with open(self.file, 'w') as fh:
                fh.writelines(lines)
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


def elink_by_step(RecIO, allIds, db, dbfrom, linkname, term=None, step=500):
    groups = splitGroup(allIds, step)
    allTarIds = RecIO.readList()
    if RecIO.nGroups == -1:
        RecIO.nGroups = len(groups)
    i = RecIO.gn
    while not RecIO.completed:
        ids = []
        logger.info(f'   Getting nucleotide ids for group {i+1}/{RecIO.nGroups}')
        handle = Entrez.elink(
            db = db,
            dbfrom = dbfrom,
            id=groups[i],
            linkname=linkname, # will remove most irrelevent sequences
            term=term
        )
        record = Entrez.read(handle)
        handle.close()
        [ids.extend(
            [l['Id'] for l in r['LinkSetDb'][0]['Link']]
        ) for r in record if len(r['LinkSetDb']) > 0]
        RecIO.appendList(ids)
        allTarIds.extend(ids)
        i += 1
    return allTarIds


def fetch_nucl_by_taxID(targetTx, minLen, targetDir, api, email, isTest=False):
    Entrez.api_key = api
    Entrez.email = email
    txHash = calHash(targetTx, isTest)
    nuclIdHash = calHash(targetTx, minLen, isTest)

    txids_file = os.path.join(targetDir, f'txid_{txHash}.txt')
    txids_io = RecordIO(txids_file)
    asmids_file = os.path.join(targetDir, f'asmb_{txHash}.txt')
    asmids_io = RecordIO(asmids_file)
    nuclids_file = os.path.join(targetDir, f'nucl_{nuclIdHash}.txt')
    nuclids_io = RecordIO(nuclids_file)

    # TODO what if you want one single tax? 0 return for the first elink? or else?

    logger.info(f'{" Start Getting Nucleotide IDs by Taxonomy ID ":=^80}')
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

    if isTest:
        allTax = allTax[:min(len(allTax),1600)]

    # 2nd step
    # Get assemblies from taxid
    logger.info(f'2. Fetch linked assembly ids of fetched taxids')
    allAssemblies = elink_by_step(asmids_io, allTax,
                                  db='assembly',
                                  dbfrom='taxonomy',
                                  linkname='taxonomy_assembly')
    logger.info(f'   Total number of assemblies is {len(allAssemblies)}.')
    logger.info(f'   The first 5 are: {allAssemblies[:min(len(allAssemblies), 5)]}')

    if isTest:
        allAssemblies = allAssemblies[:min(len(allAssemblies), 1100)]

    # 3nd step
    # Get nucleic id from refseq
    logger.info(f'3. Fetch nucl ids of fetched assemblies, only refseq.')
    allNucl = elink_by_step(nuclids_io, allAssemblies,
                            db='nuccore',
                            dbfrom='assembly',
                            linkname='assembly_nuccore_refseq',
                            term=f'{minLen}:99999999[SLEN]')
    logger.info(f'   Total number of nucleotides is {len(allNucl)}.')
    logger.info(f'   The first 5 are: {allNucl[:min(len(allNucl), 5)]}')
    return allNucl, nuclids_io.file


def fetch_nucl_by_step(targetDir, idsFile, step=200, isTest=False, maxconnections=10):
    # NOTE maxconnections>10 do not work as expected. Because of the Semaphore settings.
    # TODO make maxconnections>10 work!
    logger=logging.getLogger()
    logger.info(f'{" Start fetching nucleotide data ":=^80}')
    fetchHash = calHash(idsFile, step, isTest)
    infoFile = os.path.join(targetDir, f'fetch_{fetchHash}.pkl')
    logger = logging.getLogger('Fetching_nucleotide')
    ids = RecordIO(idsFile).readList()
    logger.info(f'Total number of ids {len(ids)}, first 5: {ids[:min(len(ids),5)]}')
    if isTest:
        ids = ids[:min(len(ids), 500)]

    idGroups = splitGroup(ids, step)
    toFetchGroupIdxs = list(range(len(idGroups)))
    finishedGroupIdxs = []

    infoLock = threading.Lock()
    def updateInfo(n, finish=True):
        with infoLock:
            if n != -1 and finish:
                finishedGroupIdxs.append(n)
                finishedGroupIdxs.sort()
                toFetchGroupIdxs.remove(n)
            elif not finish:
                toFetchGroupIdxs.append(n)
                toFetchGroupIdxs.sort()
                finishedGroupIdxs.remove(n)
            with open(infoFile, 'wb') as fh:
                pickle.dump(finishedGroupIdxs, fh)

    if os.path.isfile(infoFile):
        with open(infoFile, 'rb') as fh:
            finishedGroupIdxs = pickle.load(fh)
        [toFetchGroupIdxs.remove(n) for n in finishedGroupIdxs]
    else:
        updateInfo(-1)

    maxDig = len(str(len(idGroups)+1))

    sema = threading.Semaphore(10)
    # API key allow only 10 connections per second, will sleep 1 sec within the function to comply.
    # This may be redundant as the Entrez module has _open.previous variable to track timming:
    # https://github.com/biopython/biopython/blob/e001f2eed6a9cc6260f60f87ed90fcf8ea2df4ca/Bio/Entrez/__init__.py#L561
    def fetchGbk(n):
        sema.acquire()
        file = os.path.join(targetDir, f'nucl_{str(n+1).zfill(maxDig)}.gbk')
        idFile = os.path.join(targetDir, f'nucl_{str(n+1).zfill(maxDig)}.txt')
        trials = 3
        success = False
        logger.info(f'Fetching group No. {n+1}, {file}')
        while not success and trials > 0:
            try:
                handle = Entrez.efetch(
                    db='nuccore',
                    id=idGroups[n],
                    rettype='gbwithparts',
                    retmode='text'
                )
                with open(file, 'w') as fh:
                    record = handle.read()
                    # Example error message if rates are exceeded:
                    # {"error":"API rate limit exceeded","count":"11"}
                    # However, I don't know what is the return from error message.
                    # Str or byte? Or maybe directly HTTPError?
                    if not isinstance(record, str):
                        logger.warning('It seems you have some API problem:')
                        rec = record.decode()
                        logger.warning(rec[:min(100, len(rec))])
                        raise # go to next trial
                    else:
                        fh.write(record)
                handle.close()
                success = True
            except Exception as e:
                logger.warning(f'Error fetching group {n+1}:')
                logger.warning(f'{e}')
                sleep(2)
                trials -= 1
        if success:
            with open(idFile, 'w') as fh:
                [fh.write(str(i) + '\n') for i in idGroups[n]]
            logger.info(f"Finished group {n+1}: {os.stat(file).st_size/1024/1024:.2f} MB {os.path.split(file)[-1]}")
            updateInfo(n) # Number finished, update after fetching
        else:
            logger.info(f"3 trials not succeed for group {n+1}")
        sema.release()
        sleep(1)
    
    with ThreadPoolExecutor(max_workers=maxconnections) as worker:
        for n in toFetchGroupIdxs:
            worker.submit(fetchGbk, n)

    return toFetchGroupIdxs




if __name__ == '__main__':
    import argparse

    scriptDir = os.path.split(os.path.realpath(__file__))[0]
    testOutputDir = os.path.join(scriptDir, 'test_output')

    parser = argparse.ArgumentParser(exit_on_error=False)
    parser.add_argument('api', type=str, nargs='?', default="", help='API key from NCBI')
    parser.add_argument('email', type=str, nargs='?', default="", help='Email address for identify yourself')
    parser.add_argument('-t', action='store_true', help='Set if you want to run a test. Only few records will be collected')
    parser.add_argument('--taxIds', type=str, nargs='+', default=['201174'], help='Target taxonomyIds')
    parser.add_argument('--minLen', type=int, default=10000, help='Discard nucleotides less than this length')
    parser.add_argument('--outputDir', type=str, default=testOutputDir, help='Output dir')
    
    args = parser.parse_args()
    os.makedirs(args.outputDir, exist_ok=True)
    logging.basicConfig(filename=os.path.join(args.outputDir, 'Fetch_taxonomy_nucl.log'), level=logging.INFO)
    logger = logging.getLogger()
    argsFile = os.path.join(args.outputDir, 'args.txt')

    if args.api == "" or args.email == "":
        assert args.api == "" and args.email == "", 'Leave both api and email empty if you want to restart a previsou run.'
        if not os.path.isfile(argsFile):
            raise FileNotFoundError(f'File {argsFile} not exist. Cannot continue.')
        with open(argsFile, 'r') as fh:
            logger.info(f'Getting arguments from previous run in {argsFile}')
            argLine = fh.readline()
            try:
                args = parser.parse_args(argLine.split(' ')[1:])
            except argparse.ArgumentError as e:
                print('catch')
                logger.info('Argument parsing error. Argument line in file:')
                logger.info(argLine)
                raise e
    else:
        with open(argsFile, 'w') as fh:
            fh.write(' '.join(sys.argv))

    # Check if test is correctly set or unset
    if args.t:
        logger.info(f'Doing test fetching {args.taxIds} with minimum length {args.minLen} to output dir {testOutputDir}.')
        logger.info('Only fetch first 500 ids')
        if args.taxIds != ['201174'] or args.minLen != 10000 or os.path.realpath(args.outputDir) != testOutputDir:
            parser.error('If test (-t), please leave "taxIds", "minLen", and "outputDir" default.')
    elif args.taxIds == ['201174'] and args.minLen == 10000 and os.path.realpath(args.outputDir) == testOutputDir:
        parser.error('Arguments --taxIds, --minLen, --outputDir are required!')
    elif os.path.realpath(args.outputDir) == testOutputDir:
        parser.error(f'Please do not set output dir as {testOutputDir}')

    os.makedirs(args.outputDir, exist_ok=True)
    allNucl, nuclIdsFile = fetch_nucl_by_taxID(targetTx=args.taxIds,
                                               minLen=args.minLen,
                                               targetDir=args.outputDir,
                                               api=args.api, email=args.email,
                                               isTest=args.t)
    missingGroups = fetch_nucl_by_step(args.outputDir, nuclIdsFile, isTest=args.t)
    if not args.t:
        trials = 3
        while len(missingGroups) > 0 and trials > 0:
            trials -= 1
            retryStr = f'Retrying {len(missingGroups)} missed groups. Trial No.{3-trials}'
            logger.info(f'{retryStr:=^80}')
            missingGroups = fetch_nucl_by_step(args.outputDir, nuclIdsFile)
        if len(missingGroups) > 0:
            logger.warning(f'Finish with missing groups: {missingGroups}')

    logger.info(f'{" DONE ":=^80}')
