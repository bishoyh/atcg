#SeaSpider is a GO-conscious, sequence annotation tool developed for MetaFishNet

# SeaSpider #
SeaSpider runs on Python2.5+/LINUX, requiring a local BLAST program and MySQL database supports. A trimmed version, seaspider-lite is also provided, which does not require MySQL databases.


### What it does ###

  * search homology in MetaFishNet sequences;
  * annotate new fish sequences with gene ontology (requiring additional database support);
  * search NCBI remotely for plain gene annotations.

### Requirements ###

Requires Python 2.5 and up, which is included in most newer
Linux distributions. SeaSpider was not tested on Windows or Mac,
though the porting may be trivial.

Aslo requires a local BLAST
program, and
three sequence databases (see below).

MySQL database supports for
Gene Ontology
MetaFishNet database
a small Ensembl port


### Install and configure BLAST ###

1) Download BLAST from
[ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/)
SeaSpider was tested on version 2.2.18 and 2.2.22.
The doc/index.html in downloaded package contains how to set
up your local BLAST.

2) Download MetaFishNet sequences from
http://code.google.com/p/atcg
Unpack the file to your blast/mydb/ folder.
Run "formatdb" to create database files, e.g.
$ ../bin/formatdb -i metafishnet\_seq.fasta -p F -o T

3) Edit the configuration section in config\_spider.py,
to make "blast\_path" and "metafishnet\_seqdb" point to
your BLAST paths.



# seaspider-lite #

This is a trimmed version for mapping new sequences to
MetaFishNet genes, returning the Ensembl gene IDs for hits.
seaspider-lite requires BLAST but not relational databases.

USE: python spider3.py inputfile outputdir

### Install and configure BLAST ###

1) Download BLAST from
[ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/)
SeaSpider was tested on version 2.2.18 and 2.2.22.
The doc/index.html in downloaded package contains how to set
up your local BLAST.

2) Download MetaFishNet sequences from
http://code.google.com/p/atcg
Unpack the file to your blast/mydb/ folder.
Run "formatdb" to create database files, e.g.
$ ../bin/formatdb -i metafishnet\_seq.fasta -p F -o T

3) Edit the configuration section in seaspider-lite.py,
to make "blast\_path" and "metafishnet\_seqdb" point to
your BLAST paths.