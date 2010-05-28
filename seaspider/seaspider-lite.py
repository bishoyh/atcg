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
seaspider-lite.py
This is a trimmed version of SeaSpider for mapping new sequences to
MetaFishNet genes, returning the Ensembl gene IDs for hits.
It requires a local BLAST program, and MetaFishNet sequences.

This code is not optimized because it inherits the architecture of
the full version, which sits in Google Code SVN version control system,
http://code.google.com/p/atcg , with two additional features
    [1] annotate new fish sequences with gene ontology
        (requiring addtional database support);
    [2] search NCBI remotely for plain gene annotations.

use:
python spider3.py inputfile outputdir
"""

import sys, os, shelve
import time

from subprocess import *
from xml.etree import ElementTree # requires Python 2.5 +


#
# configuration. Modify the first two lines to your local BLAST path.
#
blast_path = '/home/yourpath/blast/bin/'
metafishnet_seqdb = '/home/yourpath/blast/mydb/metafishnet_seq.fasta'

path_blastall = blast_path + 'blastall'

GENUINE_HOMOLOGY = 33   ## only true if over 33 nt identical in alignment
BLAST_BATCH_SIZE = 200
TMPT_IN = 'blast_tmpt_input'
SPIDER_LOG = 'seaspider_log'
ANNOTATION_RESULT = 'seaspider_mfnmatch_result.txt'


#
# a BLAST parser
#
class singleparser:
    """
    parse BLAST xml result from string.
    handles multiple seq_query in one batch if newstyle.
    old style BLAST xml contains one seq_query per <BlastOutput>;
    while new style contains multiple.
    return { query_id : a list of hits, ... }, with each hit as a dict
    """
    def parse(self, s, check_homology=True, oldstyle=True):
        """
        parse BLAST xml result from string, wrapper function
        """
        tree = ElementTree.fromstring(s)
        return self.treeparse(tree, check_homology, oldstyle)
        
    def treeparse(self, tree, check_homology=True, oldstyle=True):
        """
        parsing a (XML) tree, return { query_id : hits, ... }
        """
        result_dict = {}
        iterationlist = tree.findall('BlastOutput_iterations/Iteration')
        
        if oldstyle:
            ## only one iteration
            queryid = tree.findtext('BlastOutput_query-def', '')
            result_dict[queryid] = self.parseiteration(iterationlist[0], check_homology)[1]
        else:
            for inode in iterationlist:
                (a, b) = self.parseiteration(inode, check_homology)
                result_dict[a] = b
        
        return result_dict
    
    def parseiteration(self, inode, check_homology=True):
        """
        parse an iteration node, return ('queryid', [hits])
        """
        hitlist = inode.findall('Iteration_hits/Hit')
        hits = []
        queryid = inode.findtext('Iteration_query-def', 'error_Iteration_query-def')
        for h in hitlist:
            try:
                parsed = self.parsehit(h, check_homology)
                if parsed:
                    hits.append(parsed)
            except:
                pass
                
        return (queryid, hits)

    def parsehit(self, hnode, check_homology=True):
        """
        parse a Hit node, return a dict of a number of keys,
        currently Hit_accession, Hit_def, and best E-value.
        ck_homology enforces a quality check of alignment: 
        over 33 identical nucleotides
        No error handling here; errors to be handled in self.parseiteration()
        """
        global GENUINE_HOMOLOGY
        hitdict, es = {}, []
        #es as a collector of E-values; check out validity of E-values first
        hsps = hnode.findall('Hit_hsps/Hsp')
        if check_homology:
            for h in hsps:
                if int(h.findtext('Hsp_identity')) > GENUINE_HOMOLOGY:
                    es.append(float(h.findtext('Hsp_evalue')))
        else:
            for h in hsps:
                es.append(float(h.findtext('Hsp_evalue')))
        if es:
            hitdict['Hit_accession'] = hnode.findtext('Hit_accession')
            hitdict['Hit_def'] = hnode.findtext('Hit_def', 'no def')
            hitdict['E_value'] = str(min(es))    #Best E-value
            # add more keys if you need them
        return hitdict

#
# persistent storage
#
class work_entry:
    """
    a piece of sequence
    class instance to be stored in shelved database while spider3 is working
    """
    def __init__(self, input_id='', sequence=''):
        
        self.input_id = input_id
        self.sequence = sequence
        self.homolog = ""#  accession id
        self.description = ""#  Hit_def
        self.evalue = ""

#
# main class
#
class spider3:
    """
    the annotation tool in SeaSpider
    """
    
    sp = singleparser()
    
    def dispatcher(self, args):
        """
        dispatcher subprogram according to args
        """
        workplace_db = shelve.open(os.path.join(args[1], 'workplacedb.dat'),
                                     writeback=True)
        workplace_db['records'] = []
        records = workplace_db['records']
        alist = self.read_fasta(args[0])
        for item in alist:
            entry = work_entry(input_id=item[0], sequence=item[1])
            records.append(entry)
        
        if len(args) == 2 or args[2] == '0':
            self.mfnmap(records, args[1])
        else:
            print "dispatcher error."

        workplace_db.close()

    def mfnmap(self, shelve_list, outputdir):
        """
        This option searches against the metabolic genes in MetaFishNet.
        Without GO analysis, mapping to Ensembl_id.
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
        
        for entry in worklist:
            try:
                entry.homolog = workdict[entry.input_id][0]['Hit_accession']
            except (KeyError, IndexError):
                pass
        
        out = open(os.path.join(outputdir, ANNOTATION_RESULT), 'w')
        for entry in shelve_list:
            out.write("\t".join( [entry.input_id, entry.homolog] ) + "\n")
        out.close()
    
    @staticmethod
    def read_fasta(infile):
        """
        read a fasta file, return a list
        of [(input_id, sequence), ...]
        """
        w = open(infile).read().split(">")
        alist = []
        for item in w:
            a = item.splitlines()
            if len(a) > 1:
                alist.append( (a[0].rstrip(), "".join(a[1:])) )
            elif item.strip():
                print "Not used: ", item
            else:
                pass
        print "Read input file...OK"
        return alist

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
        if blastxml:
            return self.sp.parse(blastxml, check_homology, False)
        else:
            print "Error in calling BLAST.\n", blasterr




if __name__ == "__main__":

    args = sys.argv
    if args:
        args = args[1:]
        if os.path.exists(args[1]):
            print "Output directory already exists! Please name a new one."
            
        else:
            os.mkdir(args[1])
            print "Output directory created...", args[1]
            
            spd = spider3()
            spd.dispatcher(args)

    else:
        print "USE: python seaspider-lite.py inputfile outputdir"
