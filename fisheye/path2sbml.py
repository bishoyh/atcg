"""transform a .path file to SBML
part of MetaFishNet project.
For inferred reactions, the directionality may be unknown.
In these cases, all compounds were put in reactants; 
this scenario can be spotted by empty field of products.

Last modified 04/22/2009, Shuzhao Li
"""

import libsbml

class reaction:
    """
    use EC and cpd ids instead of names.
    Use ecstr as a unit. I.e. treat an enzyme different from an enzyme complex.
    Leave the name translation and formatting to class vizpath.
    """
    def __init__(self):
        self.ID = ''
        self.ecstr = ''
        self.reactants = []
        self.products = []
        self.cpds = []
        self.source = ''
        self.pathway = ''
        self.mark = '' #infish
        self.note = ''
        self.line = ''


class pathway:
    """
    Read a pathway from .path file in MetaFishNet project;
    write a SBML file.
    """
    def __init__(self):
        self.pathway = ''
        self.rxns = []
        self.eclist = []
        self.cmpds = []
        self.filename = ''

    def read_pathfile(self, infile):
        w = open(infile).readlines()
        # path file of predefined format
        hlist = []
        for line in w:
            a = line.split('\t')
            rxn = reaction()
            rxn.ID = a[0]
            rxn.ecstr = a[1]
            rxn.pathway = a[2].replace('"', '').strip()
            rxn.source = a[3]
            rxn.reactants = [x for x in a[4].strip().split(";") if x]
            rxn.products = [x for x in a[5].strip().split(";") if x]
            rxn.cpds =  rxn.reactants + rxn.products
            try:
                rxn.mark = a[6].strip()
                rxn.note = a[7].strip()
            except IndexError:
                pass
            hlist.append(rxn)
        self.rxns = hlist
        self.eclist = list(set([x.ecstr for x in self.rxns if x.ecstr]))
        self.tally_cpds()
        try:
            self.pathway = w[0].split('\t')[2].replace('"', '').strip()
        except IndexError:
            pass
        self.filename = infile.split(".")[-2]

    def tally_cpds(self):
        cmpds = [] 
        for r in self.rxns:
            cmpds += r.cpds
        self.cmpds = list(set(cmpds))

    def write_sbml(self, outfile = ''):
        if not outfile:
            outfile = self.filename + '.xml'
        docu = libsbml.SBMLDocument()
        model = libsbml.Model(self.pathway)
        model.setAnnotation(
                    "a pathway from MetaFishNet project")
        ectype = libsbml.SpeciesType("enzyme")
        cpdtype = libsbml.SpeciesType("compound")
        model.addSpeciesType(ectype)
        model.addSpeciesType(cpdtype)
        for ec in self.eclist:
            sp = libsbml.Species(ec)
            sp.setSpeciesType("enzyme")
            model.addSpecies(sp)
        for cpd in self.cmpds:
            ss = libsbml.Species(cpd)
            ss.setSpeciesType("compound")
            model.addSpecies(ss)
            
        for rxn in self.rxns:
            sbmrxn = libsbml.Reaction(rxn.ID)
            sbmrxn.addModifier(libsbml.ModifierSpeciesReference(rxn.ecstr))
            for sub in rxn.reactants:
                sbmrxn.addReactant(libsbml.SpeciesReference(sub))
            for prd in rxn.products:
                sbmrxn.addProduct(libsbml.SpeciesReference(prd))
            model.addReaction(sbmrxn)
            
        docu.setModel(model)
        libsbml.writeSBML(docu, outfile)




if __name__ == "__main__":
    import sys
    p = pathway()
    p.read_pathfile(sys.argv[1])
    p.write_sbml()
    




