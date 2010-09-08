"""transform a .path file to SBML
part of MetaFishNet project.
For inferred reactions, the directionality may be unknown.
In these cases, all compounds were put in reactants; 
this scenario can be spotted by empty field of products.

Developed with libsbml 3.3.2, now keeping up with ver 4.0.0.
initiation of libsbml.Model() and SpeciesType, update with libsbml
Last modified 06/08/2010, Shuzhao Li
"""

import libsbml


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
        self.id = ''
        self.pathway = ''
        self.rxns = []
        self.eclist = []
        self.cmpds = []
        self.filename = ''

    def read_pathfile(self, infile):
        self.id = infile.split('.')[0]
        w = open(infile).readlines()
        # path file of predefined format
        hlist = []
        for line in w:
            a = line.split('\t')
            rxn = reaction()
            rxn.ID = a[0]
            rxn.ecstr = a[1].strip()
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
        model = docu.createModel()
        model.setId(self.id)
        model.setName(self.pathway)
        model.setAnnotation(
                    "a pathway from MetaFishNet project")
        """
        ectype = libsbml.SpeciesType("enzyme")
        cpdtype = libsbml.SpeciesType("compound")
        model.addSpeciesType(ectype)
        model.addSpeciesType(cpdtype)
        """
        ectype = model.createSpeciesType()
        ectype.setId("enzyme")
        cpdtype = model.createSpeciesType()
        cpdtype.setId("compound")

        for ec in self.eclist:
            sp = model.createSpecies()
            sp.setId(sidify(ec))
            sp.setName(ec)
            sp.setSpeciesType("enzyme")
        for cpd in self.cmpds:
            ss = model.createSpecies()
            ss.setId(sidify(cpd))
            ss.setName(cpd)
            ss.setSpeciesType("compound")
            
        for rxn in self.rxns:
            sbmrxn = model.createReaction()
            sbmrxn.setId(rxn.ID)
            m = sbmrxn.createModifier()
            if rxn.ecstr:
                m.setId(sidify(rxn.ecstr))
            #m.setName(rxn.ecstr)
            for sub in rxn.reactants:
                r = sbmrxn.createReactant()
                r.setId(sidify(sub))
            for prd in rxn.products:
                p = sbmrxn.createProduct()
                p.setId(sidify(prd))
            
        #docu.setModel(model)
        libsbml.writeSBML(docu, outfile)




if __name__ == "__main__":
    import sys
    p = pathway()
    p.read_pathfile(sys.argv[1])
    p.write_sbml()
    




