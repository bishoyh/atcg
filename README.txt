===== SeaSpider =====
This file is to be updated.


Check Google Code SVN for latest version,
http://code.google.com/p/atcg 

Features
    [0] search homology in MetaFishNet sequences only
    [1] annotate new fish sequences with gene ontology
        (requiring addtional database support);
    [2] search NCBI remotely for plain gene annotations.

This program is licensed under a permissive MIT license.


== REQUIREMENTS ==

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



== Install and configure BLAST ==

1) Download BLAST from
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/
SeaSpider was tested on version 2.2.18 and 2.2.22.
The doc/index.html in downloaded package contains how to set
up your local BLAST.

2) Download MetaFishNet sequences from
http://code.google.com/p/atcg
Unpack the file to your blast/mydb/ folder.
Run "formatdb" to create database files, e.g.
$ ../bin/formatdb -i metafishnet_seq.fasta -p F -o T

3) Edit the configuration section in seaspider-lite.py,
to make "blast_path" and "metafishnet_seqdb" point to 
your BLAST paths.


