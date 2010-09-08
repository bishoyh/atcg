""" k2sb.py,
start with KGML, extract reactions/substrate/products,
then retrieve related info form KEGG web API, write to SBML

Use KEGG API; requires SOAPpy
Use SBML, via libsbml Python binding

sample KGML:
<pathway name="path:dre00010" org="dre" number="00010"
         title="Glycolysis / Gluconeogenesis"
         image="http://www.genome.jp/kegg/pathway/dre/dre00010.gif"
         link="http://www.genome.jp/dbget-bin/show_pathway?dre00010">
         
    <reaction name="rn:R01600" type="reversible">
        <substrate name="cpd:C00221"/>
        <product name="cpd:C01172"/>
    </reaction>

Developed with libsbml 3.3.2, now keeping up with ver 4.0.0.
Last modified 09/07/2010, Shuzhao Li
"""

from SOAPpy import WSDL
import libsbml
from xml.etree import ElementTree   # requires Python 2.5 +


LEGAL_SID = [chr(x) for x in range(48, 58)
            ] + [chr(x) for x in range(65, 91)
            ] + [chr(x) for x in range(97, 123)]

def sidify(s):
    """Note that SBML has strict requirements for the syntax of identifiers. 
    The following is summary of the definition of the SBML identifier type SId 
    (here expressed in an extended form of BNF notation):
      letter ::= 'a'..'z','A'..'Z'
      digit  ::= '0'..'9'
      idChar ::= letter | digit | '_'
      SId    ::= ( letter | '_' ) idChar*
    """
    global LEGAL_SID
    if s:
        news = 'sid_'
        for x in s:
            if x not in LEGAL_SID:
                news += '_'
            else:
                news += x
        s = news
    return s

class Reaction:
    """
    use EC and cpd ids instead of names.
    Use ecstr as a unit. I.e. treat an enzyme different from an enzyme complex.
    Leave the name translation and formatting to class vizpath.
    """
    def __init__(self):
        self.ID = ''
        self.enzymes = []
        self.reactants = []
        self.products = []
        self.cpds = []


class K2SB:
    """Converting KEGG XML (KGML) to SBML, with the help of KEGG web API"""
    def __init__(self):
        wsdl = 'http://soap.genome.jp/KEGG.wsdl'
        self.serv =  WSDL.Proxy(wsdl)
        self.docu = libsbml.SBMLDocument()
    
    def convert(self, kgml, sbml):
        """input kgml file, output sbml file"""
        w = open(kgml)
        tree = ElementTree.parse(w)
        note_dict = dict(tree.getroot().items())
        model = self.docu.createModel()
        model.setName(note_dict['name'])
        model.setAnnotation(self.d2note(note_dict))
        rxns = tree.findall('reaction')
        rlist = [self.parseRxn(x) for x in rxns]
        all_enzymes, all_cpds = [], []
        for r in rlist:
            all_enzymes += r.enzymes
            all_cpds += r.cpds
        all_enzymes, all_cpds = set(all_enzymes), set(all_cpds)
        
        for ec in all_enzymes:
            sp = model.createSpecies()
            sp.setId(sidify(ec))
            sp.setName(ec)
            sp.setSpeciesType("enzyme")
        
        for cpd in all_cpds:
            ss = model.createSpecies()
            ss.setId(sidify(cpd))
            ss.setName(cpd)
            ss.setSpeciesType("compound")
            
        for rxn in rlist:
            sbmrxn = model.createReaction()
            sbmrxn.setName(rxn.ID)
            for enz in rxn.enzymes:
                m = sbmrxn.createModifier()
                m.setId(sidify(enz))
            for sub in rxn.reactants:
                r = sbmrxn.createReactant()
                r.setId(sidify(sub))
            for prd in rxn.products:
                p = sbmrxn.createProduct()
                p.setId(sidify(prd))

        libsbml.writeSBML(self.docu, sbml)

    def d2note(self, adict):
        s = ''
        for x in adict.items():
            s += '\n\t' + '='.join(x)
        return s
    
    def parseRxn(self, rxn_element):
        rxn = Reaction()
        d = dict(rxn_element.items())
        rxn.ID = d['name']
        enzymes = self.serv.get_enzymes_by_reaction(d['name'])
        rxn.enzymes = [x for x in enzymes if x]
        
        reactants = rxn_element.findall('substrate')
        rxn.reactants = [dict(x.items())['name'] for x in reactants]
        products = rxn_element.findall('product')
        rxn.products = [dict(x.items())['name'] for x in products]
        
        rxn.cpds = rxn.reactants + rxn.products
        return rxn
    
