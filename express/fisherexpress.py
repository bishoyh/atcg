# -*- coding: utf-8 -*-
# under GNU GPL license as it contains third party code

"""
FisherExpress,
pathway enrichment analysis for MetaFishNet.
There are variants that work with MFN webserver, 
or combining gene and metabolites,
and with visualization features.

"""

from math import log, exp
#from dbmodels import *

TOTAL_ECS = 821
# this is enzymes in zebrafish, assuming same number for other fish
PATHWAY_FILE = 'fishpathway_1.9.6.txt'

from extgene2ec import ext2ec
from shmgene2ec import shmgene2ec
#     "C15_01_D02": ('2.7.4.3',),


def get_ecnums_by_gene(gene, species):
    try:
        if species == 'zebrafish':
            return ext2ec[gene]
        elif species == 'sheepshead minnow':
            return list(shmgene2ec[gene])
        else:
            return []
    except KeyError:
        return []

class FisherExactTest:
    """
    This implementation of Fisher's Exact Test was written by Aurélien Mazurie.
    Adapted from WordHoard project - http://wordhoard.northwestern.edu
    (edu.northwestern.at.utils.math.statistics.FishersExactTest)
    """
    def pvalue (self, k, n, C, G):
        um, lm = min(n, C), max(0, n + C - G)
        if (um == lm):
            return 1.0, 1.0, 1.0
        cutoff = self.__hypergeometric_probability(k, n, C, G)
        left_tail, right_tail, two_tailed = 0, 0, 0
        for i in range(lm, um + 1):
            p = self.__hypergeometric_probability(i, n, C, G)
            if (i <= k):
                left_tail += p
            if (i >= k):
                right_tail += p
            if (p <= cutoff):
                two_tailed += p
        left_tail = min(left_tail, 1)
        right_tail = min(right_tail, 1)
        two_tailed = min(two_tailed, 1)
        return left_tail, right_tail, two_tailed

    def enrichment (self, k, n, C, G):
        return (float(k) / C) / (float(n) / G)

    def evaluate (self, k, n, C, G):
        return self.enrichment(k, n, C, G), self.pvalue(k, n, C, G)

    def __hypergeometric_probability (self, i, n, C, G):
        return exp(
          self.__lncombination(C, i) +
          self.__lncombination(G - C, n - i) -
          self.__lncombination(G, n)
         )

    # Logarithm of the number of combinations of 'n' objects taken 'p' at a time
    def __lncombination (self, n, p):
        return \
          self.__lnfactorial(n) - \
          self.__lnfactorial(p) - \
          self.__lnfactorial(n - p)

    # Logarithm of n! with algorithmic approximation
    # Reference:
    #   Lanczos, C. 'A precision approximation of the gamma function',
    #   J. SIAM Numer. Anal., B, 1, 86-96, 1964."
    #   http://www.matforsk.no/ola/fisher.htm 
    def __lnfactorial (self, n):
        if (n <= 1):
            return 0
        else:
            return self.__lngamma(n + 1)

    def __lngamma (self, z):
        x = 0
        x += 0.1659470187408462e-06 / (z + 7)
        x += 0.9934937113930748e-05 / (z + 6)
        x -= 0.1385710331296526 / (z + 5)
        x += 12.50734324009056 / (z + 4)
        x -= 176.6150291498386 / (z + 3)
        x += 771.3234287757674 / (z + 2)
        x -= 1259.139216722289 / (z + 1)
        x += 676.5203681218835 / (z)
        x += 0.9999999999995183

        return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class Pathway:
    def linefeed(self, line):
        a = line.rstrip().split('\t')
        self.id = a[0]
        self.name = a[1].replace('"', '')
        self.rxns = a[2].split(';')
        self.enzymes = a[3].split(';')
        self.num_enzymes = int(a[4])
        self.cpds = a[5].split(';')
        self.cpd_num = int(a[6])
        


class FisherExpress:
    def __init__(self):
        self.FET = FisherExactTest()

    def get_pathways(self, PATHWAY_FILE):
        self.pathways = []
        w = open(PATHWAY_FILE).readlines()
        for line in w:
            if line.strip():
                P = Pathway()
                P.linefeed(line)
                self.pathways.append(P)
        print "Finished pathway read."

    def gather(self, inputlist, species='zebrafish'):
        """
        process input list
        """
        self.ecdict = {}
        self.eclist = []
        self.inputlist = set([x for x in inputlist if x])
        for gene in self.inputlist:
            ecs = get_ecnums_by_gene(gene, species)
            self.eclist += ecs
            for ec in ecs:
                if self.ecdict.has_key(ec):
                    self.ecdict[ec].append(gene)
                else:
                    self.ecdict[ec] = [gene]
        self.eclist = set(self.eclist)
        self.n = len(self.eclist)
        
    def ec_enrich_test(self, pathway_instance_list):
        """
        returns (p-value, num_overlap, overlap_features, pathway_instance)
        """
        result = []
        for P in pathway_instance_list:
            overlap_features = self.eclist.intersection( set(P.enzymes) )
            overlap_size = len( overlap_features )
            enrich_pvalue = self.FET.pvalue(overlap_size, self.n,
                                       P.num_enzymes, TOTAL_ECS)[1] 
                                       #ec_num as pathway size
            result.append( (enrich_pvalue, overlap_size, overlap_features, P) )
        return result
        
    def html_output(self, result):
        """
        to be fixed.
        """
        result.sort()
        counter = 0
        header = "<p>Input %d features, converted to %d enzymes. \
                 Number of total features is %d. </p>" %(self.n,
                 len(self.eclist), self.array_size)
        html = ''
        for r in result:
            counter += 1
            if r[0] < 0.05 or (r[0] < 0.5 and counter < 6):
                s = '<h3><a href="' + r[3].link + '">' + r[3].name + \
                '</a>, p-value=' + str(r[0]) + '</h3><p>overlap_size: ' +\
                str(r[1]) + ', pathway_size: ' + str(r[3].num_enzymes) +\
                '</p><p>Hits on this pathway:<br>'
                for ec in r[2]:
                    s += ec + ': ' + ' '.join(self.ecdict[ec]) + '<br>'
                s += '</p>'
                html += '<div class="analyzed_block">' + s + '</div>'
        return header + (html or 'No significant hit.')


    def write_mfn_result(self, result, outfile):
        # mfn_pathway   sig_ecs total_ecs   p-value
        # to add genes
        s = 'mfn_pathway\tselected_enzymes\tenzymes_in_pathway\tp-value\tECs\n'
        for r in result:
            s += '\t'.join([r[3].name, str(r[1]), str(r[3].num_enzymes), str(r[0]), 
                    ';'.join(r[2]) ]) + '\n'
            
        f = open(outfile, 'w')
        f.write(s)
        f.close()

    def f2f(self, infile, outfile, species='sheepshead minnow'):
        """
        wrapper function
        """
        self.get_pathways(PATHWAY_FILE)
        inputlist = open(infile).readlines()
        inputlist = [x.strip() for x in inputlist if x.strip()]
        self.gather(inputlist, species)
        result = self.ec_enrich_test(self.pathways)
        result.sort()
        self.write_mfn_result(result, outfile)
        

if __name__ == '__main__':
    
    import sys
    if len(sys.argv[1:]) < 3:
        print "USAGE: python fisherexpress.py infile outfile dre/shm"
    else:
        infile, outfile = sys.argv[1], sys.argv[2]
        dredict = {'dre': 'zebrafish', 'shm': 'sheepshead minnow'}
        FE = FisherExpress()
        FE.f2f(infile, outfile, species=dredict[sys.argv[3]])














