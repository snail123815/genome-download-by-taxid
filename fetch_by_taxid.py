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


Entrez.api_key = input("Your api key for EntrezSearch:")
Entrez.email = input("Your email address for EntrezSearch:")

targetTax = '201174' # Actinobacteria, some super tax
# TODO what if you want one single tax? 0 return for the first elink? or else?

handle = Entrez.elink(
    db='taxonomy',
    dbfrom='taxonomy',
    id=targetTax,
    linkname='taxonomy_taxonomy_child',
)
record_t = Entrez.read(handle)
handle.close()

print(len(record_t))

print(record_t[0]['LinkSetDb'][0].keys())
print(record_t[0]['LinkSetDb'][0]['LinkName'])
x = record_t[0]['LinkSetDb'][0]['Link']
allTax = [i['Id'] for i in x]
print(len(allTax))

print(allTax[:5])


handle = Entrez.elink(
    db='assembly',
    dbfrom='taxonomy',
    id=allTax[:100],
    linkname='taxonomy_assembly',
)
record_a = Entrez.read(handle)
handle.close()

print(len(record_a))
assemblies = [r['LinkSetDb'][0]['Link'] if len(r['LinkSetDb']) > 0 else [] for r in record_a]
lenAssemblies = [len(x) for x in assemblies]
print(lenAssemblies)
allAssemblies = []
[allAssemblies.extend([l['Id'] for l in r['LinkSetDb'][0]['Link']]) for r in record_a if len(r['LinkSetDb']) > 0]
print(len(allAssemblies))
n = lenAssemblies.index(209)
print(record_a[n])


handle = Entrez.elink(
    db = 'nuccore',
    dbfrom = 'assembly',
    id=allAssemblies,
    linkname='assembly_nuccore_refseq', # will remove most irrelevent sequences
    # linkname='assembly_nuccore',
    term='10000:99999999[SLEN]'
)
record_n = Entrez.read(handle)
handle.close()

x = record_n[1]['LinkSetDb']
print(allAssemblies[1])
print(x)
r = [r['LinkSetDb'][0]['Link'] if len(r['LinkSetDb']) > 0 else [] for r in record_n]
rl = [len(a) for a in r]
allnucl = []
[allnucl.extend([l['Id'] for l in a]) for a in r if len(r)>0]

print(rl)
print(len(allnucl))
