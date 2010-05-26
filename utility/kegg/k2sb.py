""" k2sb.py,
start with KGML, extract reactions/substrate/products,
then retrieve related info form KEGG web API,
write to SBML

Use KEGG API; requires SOAPpy
Use SBML, via libsbml Python binding

"""

from SOAPpy import WSDL
import libsbml
from xml.etree import ElementTree   # requires Python 2.5 +

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
        self.model = libsbml.Model(note_dict['name'])
        self.model.setAnnotation(self.d2note(note_dict))
        rxns = tree.findall('reaction')
        for rxn in rxns:
            self.model.addReaction(self.ParseRxn(rxn))
        self.docu.setModel(self.model)
        libsbml.writeSBML(self.docu, sbml)

    def d2note(self, adict):
        s = ''
        for x in adict.items():
            s += '\n\t' + '='.join(x)
        return s
    
    def ParseRxn(self, rxn_element):
        """parse a reaction from reaction tree_element, 
        return an instance of libsbml.Reaction"""
        d = dict(rxn_element.items())
        #[('type', 'reversible'), ('name', 'rn:R00710')]
        print d['name']
        # name, reversible
        rxn = libsbml.Reaction(d['name'])
        #, "", 0, d['type'] == 'reversible')
        # get all enzymes thru web API, as modifier
        enzymes = self.serv.get_enzymes_by_reaction(d['name'])
        for enz in enzymes:
            sp = libsbml.Species(enz)
            self.model.addSpecies(sp)
            rxn.addModifier(libsbml.ModifierSpeciesReference(enz))
        # get all substrates from tree_element, as reactants
        subs = rxn_element.findall('substrate')
        for sub in subs:
            sub_name = dict(sub.items())['name']
            ss = libsbml.Species(sub_name)
            self.model.addSpecies(ss)
            rxn.addReactant(libsbml.SpeciesReference(sub_name))
        # get all products from tree_element, as products
        prds = rxn_element.findall('product')
        for prd in prds:
            prd_name = dict(prd.items())['name']
            ss = libsbml.Species(prd_name)
            self.model.addSpecies(ss)
            rxn.addProduct(libsbml.SpeciesReference(prd_name))
        return rxn







""" sample KGML:
<pathway name="path:dre00010" org="dre" number="00010"
         title="Glycolysis / Gluconeogenesis"
         image="http://www.genome.jp/kegg/pathway/dre/dre00010.gif"
         link="http://www.genome.jp/dbget-bin/show_pathway?dre00010">
         
    <reaction name="rn:R01600" type="reversible">
        <substrate name="cpd:C00221"/>
        <product name="cpd:C01172"/>
    </reaction>
"""


