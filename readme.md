# Fetch full genbank file by taxonomy ID

Given a taxonomy ID, the script can download all its children taxonomy genome files.

```mermaid
flowchart LR
  TaxID --> Child[Children TaxIDs] --> Assemblies --> NuclIDs[Nucleotide IDs] --> files[Genbank files]
```

[TOC]

## Usage

Use `-h` to get help:

Note `api` and `email` are essential. Make sure you have an Entrez API key before you start.

```console
» python3 fetch_nucl_by_taxid.py -h
usage: fetch_nucl_by_taxid.py [-h] [-t] [--taxIds TAXIDS [TAXIDS ...]] [--minLen MINLEN] [--outputDir OUTPUTDIR] [api] [email]

positional arguments:
  api                   API key from NCBI
  email                 Email address for identify yourself

optional arguments:
  -h, --help            show this help message and exit
  -t                    Set if you want to run a test. Only few records will be collected
  --taxIds TAXIDS [TAXIDS ...]
                        Target taxonomyIds
  --minLen MINLEN       Discard nucleotides less than this length
  --outputDir OUTPUTDIR
                        Output dir
```

You can monitor the status by:

```console
tail -F test_output/Fetch_taxonomy_nucl.log
```

## Test

```console
python3 fetch_nucl_by_taxid.py -t NCBI_ENTREZ_API_KEY e@mail.com
```

It will create `test_output` dir and start download few files. Check the result log file with:

```console
tail -F test_output/Fetch_taxonomy_nucl.log
```

Until it says `DONE` in the last line.
