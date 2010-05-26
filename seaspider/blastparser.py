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
blastparser.py, part of SeaSpider,
parses BLAST xml results from string, and
handles results from single query or multiple queries
(blast with multiple entries is faster than doing separately).

Dealing with two xml formats from BLAST
blastall uses a new style, 
    where <Iteration_query-def> is defined; 
    each <Iteration> reports query on one sequence
blastcl3 uses an old style, 
    where each <BlastOutput> reports query on one sequence

## Example of data structure:
    output dict
    {queryid: hits, ...}
    where
    hits = [{'Hit_accession': 'CR650049',
             'Hit_def': 'Tetraodon nigroviridis full-length cDNA',
             'E_value': '1.03788e-79',
             ...}, {  }, ...]

last updated: 04/20/2009 by Shuzhao Li
"""


from xml.etree import ElementTree # requires Python 2.5 +
from wordlists import *

GENUINE_HOMOLOGY = 33   ## only true if over 33 nt identical in alignment


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



class multiparser(singleparser):
    """
    for parsing blastcl3 results, old style 
    where each seq_query is divided to a <BlastOutput>
    """
    def parse(self, s, check_homology=True, oldstyle=True):
        """
        This was written for parsing returns from blastcl3.
        blastcl3/NCBI returns multiple xml roots in one file (old style?)
        so cut it up to a list of xml strings, and process one by one.
        Return a list of tuples (queryid, hits).
        """
        result_dict = {}
        blocks = self.part_xml(s)
        #print "We've got %s blocks." %len(blocks)
        for b in blocks:
            node = ElementTree.fromstring(b)
            d = self.treeparse(node, True, True)
            result_dict.update(d)
        return result_dict
        
    def part_xml(self, s):
        """
        separate multiple entries to individual xml strings
        """
        blocks, b = [], 1
        while b > 0:
            a = s.find('<?xml version=')
            b = s[a+1:].find('<?xml version=')
            blocks.append(s[a:b+1])
            s = s[b+1:]
        # the last loop added '' to blocks
        blocks[-1] = s
        return blocks


#
# utility function(s)
#

def check_meaningful(s):
    """
    check a description str, if meaningful as annotation
    yet a rudimentary function, will improve in future
    """
    stamp = True
    for phrase in bad_phrases:
        if phrase in s:
            stamp = False
            break

    if stamp:
        alist = s.lower().split()
        for item in alist:
            if item in negative_list:
                stamp = False
                break
    
        for item in alist:
            if item in positive_list:
                stamp = True

    return stamp




if __name__ == '__main__':
    testxml = 'output2.xml' #a test
    sp = singleparser()
    w = open(testxml).read()
    print sp.parse(w)
