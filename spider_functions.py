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
spider_functions.py

annotation tool for nt sequence libraries (e.g. EST)
either by Gene Ontology or homology on NCBI.

Part of SeaSpider,
which can be used as a standalone tool for annotations with minor tweaking.

Requires 
local blast installation (and my parser tool blastparser.py)
(blastall for local search, blastcl3 client for remote search on NCBI).
and locally installed GO-related databases

last updated: 04/20/2009 by Shuzhao Li
"""

import time, MySQLdb
from subprocess import *

from blastparser import *
from config_spider import *
from go_sleek_list import go_sleek

sp = singleparser()
mp = multiparser()


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
        else:
            print "Not used: ", item
    print "Read input file...OK"
    return alist



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
        
        # GO terms directly out of existing annotation on assicated homologs
        self.raw_go_terms = []
        
        # to obtain by self.go_analyze
        self.complete_go_terms = []
        self.end_go_terms = []  #children only
        self.go_sleek_terms = []  #possibly multiple entries
        self.go_sleek_desc = []
        self.is_metabolic = False

    def go_analyze(self, cursor):
        """
        process raw_go_terms to complete_go_terms and go_sleek
        by querying local GO database via cursor
        """
        
        # get all ancestors
        query_templ = "SELECT DISTINCT ancestor.acc FROM term \
            INNER JOIN graph_path ON (term.id=graph_path.term2_id)  \
            INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id) \
            WHERE term.acc = '%s'"
        ancestors = []
        for item in self.raw_go_terms:
            cursor.execute(query_templ % item)
            result = cursor.fetchall()
            for x in result:
                if x[0] and ":" in x[0] and x[0] != item:
                    ancestors.append(x[0])

        ancestors = list(set(ancestors))
        for item in self.raw_go_terms:
            if item not in ancestors:
                self.end_go_terms.append(item)
        
        self.complete_go_terms = list(set(self.raw_go_terms + ancestors))
        if "GO:0008152" in self.complete_go_terms:
            self.is_metabolic = True
        
        # the sleek/best GO terms
        for item in self.complete_go_terms:
            if item in go_sleek:
                self.go_sleek_terms.append(item)
        if not self.go_sleek_terms:
            self.go_sleek_terms = self.end_go_terms
        
        # fill in self.go_sleek_desc
        self.go_desc(cursor)

    def go_desc(self, cursor):
        """
        fill in self.go_sleek_desc
        """
        for acc in self.go_sleek_terms:
            query_templ = "select name from term where acc = '%s'"
            cursor.execute(query_templ % acc)
            try:
                self.go_sleek_desc.append( cursor.fetchall()[0][0] )
            except IndexError:
                print "Error in term query ", acc
                self.go_sleek_desc.append( "" )


