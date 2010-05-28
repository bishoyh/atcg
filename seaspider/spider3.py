#!/usr/bin/env python
# Copyright (c) 2008-2010 Shuzhao Li.
#
# Licensed under the MIT License.
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.



"""
spider3.py
The annotator with MetaFishNet/SeaSpider, three options:
    [0] search homology in MetaFishNet sequences only
        (metabolic genes, from 5 fish genomces)
    [1] fish ontology,
        GO annotation by finding homologs in zebrafish then GO generic, locally
    [2] search NCBI remotely, no GO, descriptive only. 

Use:
python spider3.py inputfile outputdir [annotation_option]

Note:
SeaSpider needs support of sequence databases and a few MySQL databases 
to accomplish its tasks. The links to download databases are given at 
http://metafishnet.appspot.com .
A shelved database is used to store data while spider3 is working.

05/28/2010
Fully functional with BLAST ver 2.1.18. However, blastcl3 
may be deprecated in the latest BLAST.
The remote search option may need rewrite.
"""

import sys, os, shelve
from spider_functions import *

BLAST_BATCH_SIZE = 200

class spider3:
    """
    the annotation tool in SeaSpider
    """
    def dispatcher(self, args):
        """
        dispatcher subprogram according to args
        """
        workplace_db = shelve.open(os.path.join(args[1], 'workplacedb.dat'),
                                     writeback=True)
        workplace_db['records'] = []
        records = workplace_db['records']
        alist = read_fasta(args[0])
        for item in alist:
            entry = work_entry(input_id=item[0], sequence=item[1])
            records.append(entry)
        
        if len(args) == 2 or args[2] == '0':
            self.mfnmap(records, args[1])
        elif args[2] == '1':
            self.fish_ontology(records, args[1])
        elif args[2] == '2':
            self.remote_desc(records, args[1])
        else:
            print "annotation options should be one of [0, 1, 2]."

        workplace_db.close()

    def mfnmap(self, shelve_list, outputdir):
        """
        This option searches against the metabolic genes in MetaFishNet.
        Without GO analysis, mapping to Ensembl_id and descriptions.
        """
        tmpt_input = os.path.join(outputdir, TMPT_IN)
        pcmd = path_blastall + ' -p blastn -d ' + metafishnet_seqdb \
                + ' -i ' + tmpt_input + ' -e 1E-5 -m 7'
        
        worklist = [entry for entry in shelve_list if not entry.homolog]
        workdict = {}
        batchnum = int(len(worklist)/BLAST_BATCH_SIZE) + 1
        for ii in range(batchnum):
            ## process in batch of BLAST_BATCH_SIZE entries
            self.make_input(worklist[ii*BLAST_BATCH_SIZE:
                    ii*BLAST_BATCH_SIZE+BLAST_BATCH_SIZE], tmpt_input)
            print "processing metafishnet blast batch %d of %d" %(ii+1, batchnum)
            dii = self.blast(pcmd)
            workdict.update(dii)
        
        dbc = MySQLdb.connect(user='guest', db=LOCAL_MFN_DB)
        cursor = dbc.cursor()
        print "updating annotations from MetaFishNet database......."
        for entry in worklist:
            try:
                hit = workdict[entry.input_id][0]['Hit_accession']
                entry.homolog = hit
                entry.description = self.get_mfn_desc(cursor, hit)
            except (KeyError, IndexError):
                pass
        
        dbc.close()
        print "Finished metafishnet search."
        
        out = open(os.path.join(outputdir, ANNOTATION_RESULT), 'w')
        unfound = []
        for entry in shelve_list:
            if entry.homolog:
                out.write("\t".join( [entry.input_id, entry.homolog, entry.description,
                                        ] ) + "\n")
            else:
                unfound.append(entry.input_id)
        out.close()
        out = open(os.path.join(outputdir, 'unfound.txt'), 'w')
        out.write("\n".join(unfound) + "\n")
        out.close()
    
    def fish_ontology(self, shelve_list, outputdir, m=True):
        """
        hardwired process of doing fish gene ontology annotation
        follow zebrafish first, then generic GO
        note the seq_db and relational db are hard coded in this section!
        m for metabolic genelist output option.
        """

        #out_log = open(os.path.joint(outputdir, SPIDER_LOG), 'a')
        tmpt_input = os.path.join(outputdir, TMPT_IN)
        
        #
        # zebrafish
        #
        pcmd = path_blastall + ' -p blastn -d ' + zebrafish_seqdb \
                + ' -i ' + tmpt_input + ' -e 1E-5 -m 7'

        worklist = [entry for entry in shelve_list if not entry.raw_go_terms]
        workdict = {}
        batchnum = int(len(worklist)/BLAST_BATCH_SIZE) + 1
        for ii in range(batchnum):
            ## process in batch of BLAST_BATCH_SIZE entries
            self.make_input(worklist[ii*BLAST_BATCH_SIZE:
                    ii*BLAST_BATCH_SIZE+BLAST_BATCH_SIZE], tmpt_input)
            print "processing zebrafish blast batch %d of %d" %(ii+1, batchnum)
            dii = self.blast(pcmd)
            workdict.update(dii)
        
        dbc = MySQLdb.connect(user='guest', db=LOCAL_ENSEMBL_DB)
        cursor = dbc.cursor()
        print "updating annotations from database......."
        for entry in worklist:
            try:
                hit = workdict[entry.input_id][0]['Hit_accession']
                entry.homolog = hit
                returned_go = self.get_zebrafish_go(cursor, hit)
                entry.raw_go_terms = [x[0] for x in returned_go if x[0]]
                entry.description = self.get_zebrafish_desc(cursor, hit)
            except (KeyError, IndexError):
                pass
        
        dbc.close()
        print "Finished zebrafish search."
        
        #
        # generic
        #
        pcmd = path_blastall + ' -p blastx -d ' + generic_seqdb \
                + ' -i ' + tmpt_input + ' -e 1E-5 -m 7'

        worklist = [entry for entry in shelve_list if not entry.raw_go_terms]
        workdict = {}
        batchnum = int(len(worklist)/BLAST_BATCH_SIZE) + 1
        for ii in range(batchnum):
            ## process in batch of BLAST_BATCH_SIZE entries
            self.make_input(worklist[ii*BLAST_BATCH_SIZE: 
                    ii*BLAST_BATCH_SIZE+BLAST_BATCH_SIZE], tmpt_input)
            print "processing generic_go blast batch %d of %d" %(ii+1, batchnum)
            dii = self.blast(pcmd)
            workdict.update(dii)
        
        # local GO db connect for both raw_go and next step of collective results
        dbc = MySQLdb.connect(user='guest', db=LOCAL_GO_DB)
        cursor = dbc.cursor()
        print "updating annotations from GO database......."
        for entry in worklist:
            try:
                hit = workdict[entry.input_id][0]['Hit_accession']
                entry.homolog = hit
                returned_go = self.get_generic_go(cursor, hit)
                entry.raw_go_terms = [x[0] for x in returned_go if x[0]]
                entry.description = self.get_go_desc(cursor, hit)
            except (KeyError, IndexError):
                pass
        
        print "Finished generic GO seq db search."
        
        self.write_result(shelve_list, outputdir, cursor, m)
        dbc.close()


    def write_result(self, shelve_list, outputdir, cursor, m=True):
        """
        export collective results; if m:
        also write a separate file for metabolic genes
        """
        out = open(os.path.join(outputdir, ANNOTATION_RESULT), 'w')
        for entry in shelve_list:
            if entry.raw_go_terms:
                entry.go_analyze(cursor)#   cursor already openned for LOCAL_GO_DB
                out.write("\t".join( [entry.input_id, entry.homolog, entry.description,
                             ",".join(entry.go_sleek_terms),
                             "$".join(entry.go_sleek_desc)] ) + "\n")
        out.close()
        print "Wrote ", ANNOTATION_RESULT
        if m:
            out = open(os.path.join(outputdir, "metabolic_genes"), 'w')
            for entry in shelve_list:
                if entry.is_metabolic:
                    out.write("\t".join( [entry.input_id, entry.homolog, entry.description,
                                ",".join(entry.go_sleek_terms),
                                "$".join(entry.go_sleek_desc)] ) + "\n")
            out.close()
        print "Wrote ", "metabolic_genes"

    
    def remote_desc(self, shelve_list, outputdir):
        """
        annotation by querying NCBI remotely;
        relaxed to E-value -5, no check_homology.
        !!Since blastcl3 is deprecated now, this function needs update!!
        """
        tmpt_input = os.path.join(outputdir, TMPT_IN)
        pcmd = path_blastcl3 + ' -p blastn -d "nr" -i ' + tmpt_input \
                             + ' -e 1E-5 -v 5 -b 5 -m 7'
        
        workdict = {}
        batchnum = int(len(shelve_list)/BLAST_BATCH_SIZE) + 1
        for ii in range(batchnum):
            ## process in batch of BLAST_BATCH_SIZE entries
            self.make_input(shelve_list[ii*BLAST_BATCH_SIZE: 
                    ii*BLAST_BATCH_SIZE+BLAST_BATCH_SIZE], tmpt_input)
            print "blasting batch %d of %d" %(ii+1, batchnum)
            dii = self.blast(pcmd, False, True)
            workdict.update(dii)
            
        print "Finished queries. Processing..."
        for entry in shelve_list:
            try:
                hits = workdict[entry.input_id]
                for h in hits:
                    if check_meaningful(h['Hit_def']):
                        entry.homolog = h['Hit_accession']
                        entry.description = h['Hit_def']
                        entry.evalue = h['E_value']
                        break
            except KeyError, IndexError:
                pass

        out = open(os.path.join(outputdir, ANNOTATION_RESULT), 'w')
        for entry in shelve_list:
            if entry.description:
                out.write("\t".join( [entry.input_id, entry.homolog,
                                 entry.description, entry.evalue] ) + "\n")
        out.close()
        print "Done."


    def make_input(self, entry_list, tmpt_input):
        """
        write a temporary input file for blast
        """
        t = open(tmpt_input, 'w')
        for entry in entry_list:
            t.write('>' + entry.input_id + '\n' + entry.sequence + '\n')
        t.close()

    def blast(self, pcmd, check_homology=True, cl3=False):
        """
        this is a BLAST wrapper, run pcmd and return hits dict.
        cl3 uses old style xml, parsed differently
        """
        p=Popen(pcmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        # use Popen.communicate() to avoid OS pipe buffer problems 
        (blastxml, blasterr) = p.communicate()
        # if blastxml, check error in blast run
        if cl3:
            return mp.parse(blastxml, check_homology, True)
        else:
            return sp.parse(blastxml, check_homology, False)

    #
    # database functions for zebrafish and generic_go
    #

    def get_zebrafish_go(self, cursor, hit):
        """
        get GO terms for hits from local ensembl (simplified) db
        """
        query_templ = "select gene_ontology_id from dre_go where ensembl_id='%s'"
        cursor.execute(query_templ % hit)
        return cursor.fetchall()
    
    def get_zebrafish_desc(self, cursor, hit):
        """find the description from local ensembl (simplified) db"""
        query_templ = "select description from dre_desc where ensembl_id = '%s'"
        cursor.execute(query_templ % hit)
        return cursor.fetchall()[0][0] or ""

    def get_mfn_desc(self, cursor, hit):
        """find the description from MetaFishNet db"""
        query_templ = "select description from gene where ensembl_id = '%s'"
        cursor.execute(query_templ % hit)
        return cursor.fetchall()[0][0] or ""

    def get_generic_go(self, cursor, seq_id):
        """find the GO terms by seq_id in local GO db"""
        # will make a temporary table to speed this up
        #
        query_templ = "select term.acc from gene_product \
                inner join gene_product_seq on (gene_product.id=gene_product_seq.gene_product_id) \
                inner join seq on (seq.id=gene_product_seq.seq_id) \
                inner join association on (gene_product.id=association.gene_product_id) \
                inner join term on (association.term_id=term.id) \
                where seq.id = %s"
        cursor.execute(query_templ % seq_id)
        return cursor.fetchall()

    def get_go_desc(self, cursor, seq_id):
        """find the seq.description by seq_id in local GO db"""
        query_templ = "select description from seq where id = %s"
        cursor.execute(query_templ % seq_id)
        return cursor.fetchall()[0][0] or ""






if __name__ == "__main__":

    args = sys.argv
    if len(args) > 1:
        # got arguments, proceed
        args = sys.argv[1:]
        # [inputfile, outputdir, annotation_option]
        if os.path.exists(args[1]):
            print "Output directory already exists! Please name a new one."
            
        else:
            os.mkdir(args[1])
            print "Output directory created...", args[1]
            
            spd = spider3()
            spd.dispatcher(args)

    else:
        print """usage:
        python spider3.py inputfile outputdir [annotation_option]
        annotation_option should be one of [0, 1, 2]
        default is 0, for searching MetaFishNet sequences;
        1 for new annotation based on zebrafish and GO; 
        2 for annotation via remote NCBI blast only
            """
