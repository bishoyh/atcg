"""
config_spider.py
configurations for spider3 (SeaSpider)
Edit environmental parameters if you need to.

last update: 12/20/2009
"""


## ----------------------------------------------------------------------
## Local installation parameters of BLAST (blastall, blastcl3)
## consistent with .ncbirc file
##

blast_path = '/home/habpi/blast/bin/'
path_blastall = blast_path + 'blastall'
# blastcl3 is deprecated
path_blastcl3 = blast_path + 'blastcl3'

## ----------------------------------------------------------------------
## locally compiled sequence DB for BLAST use, done by formatdb

metafishnet_seqdb = '/home/habpi/blast/mydb/metafishnet_seq.fasta'
zebrafish_seqdb = '/home/habpi/blast/mydb/dre_seq_en56.fasta'
generic_seqdb = '/home/habpi/blast/mydb/goseq.fasta'


## ----------------------------------------------------------------------
## LOCAL_GO_DB is a local copy of Gene Ontology database
## e.g. MySQLdb.connect(user='guest', db='go1123')
##

LOCAL_GO_DB = 'go20091220'
LOCAL_ENSEMBL_DB = 'ensembl56'
LOCAL_MFN_DB = 'fmdb'

## ----------------------------------------------------------------------
## temporary files names, and output names, all under outputdir

TMPT_IN = 'blast_tmpt_input'
SPIDER_LOG = 'spider_log'
ANNOTATION_RESULT = 'seaspider_annotation_result.txt'




