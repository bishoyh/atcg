"""metabolicnet.py
metabolic network computing and visualization.
Identifiers mainly follow conventions in KEGG.
Their descriptions are supplemented by two external files,
dict_ec_def and dict_cpds_def.

Developed with libsbml 3.3.2, now keeping up with ver 4.0.0.
In using SBML files, names instead of ids are used,
because EC numbers can no longer be used as ids in new SBML specs.

"""

import pygraphviz as pgv
import networkx
# using dev verion 1, 1192; latest version conflict
import libsbml
import random
from numpy import mean

from dict_ec_def import *
from dict_cpds_def import *


# currency metabolites
# G11113 = C00008 = ADP
currency = ['C00001', 'C00080', 'C00007', 'C00006', 'C00005', 'C00003',
            'C00004', 'C00002', 'C00013', 'C00008', 'C00009', 'C00011', 
            'G11113', '',
            'H2O', 'H+', 'Oxygen', 'NADP+', 'NADPH', 'NAD+', 'NADH', 'ATP', 
            'Pyrophosphate', 'ADP', 'Orthophosphate', 'CO2',]

shapeslist = ['diamond', 'octagon', 'invtrapezium', 'triangle', 'box', 
        'trapezium', 'invtriangle', 'parallelogram', 'polygon','egg']

#
# a few utilities, to update into mnetwork()
#
def overlap(list1, list2):
    return set(list1).intersection(set(list2))
def quote(s):
    if s[0] == '"':
        return s
    else:
        return '"'+s+'"'
def count(cpd, edges):
    counter = 0
    for edge in edges:
        if cpd in edge:
            counter += 1
    return counter
def findhubs(inlist):
    newlist = []
    for x in inlist:
        if inlist.count(x) > 3 and x not in newlist:
            newlist.append(x)
    return newlist
def flatten(inlist):
    flattened = []
    for t in inlist:
        flattened.append(t[0])
        flattened.append(t[1])
    return flattened



class mnode:
    """
    A custom class of node, mapped in a mnetwork by internal dict.
    'concise' style shows EC numbers and compound numbers;
    'descriptive' style shows enzyme and compund names.
    """
    def __init__(self, id='', type='', vizstyle='concise'):
        self.id = id
        self.speciesType = type
        self.desc = ''
        self.label = id
        self.vizstyle = vizstyle

    def set_attr(self):
        # set attributes according to node type
        if self.speciesType == 'enzyme':
            if self.vizstyle == 'concise':
                self.beautified = self.beautifyec()
            else:
                self.beautified = self.beautifyec_desc()
            self.shape = 'ellipse'
            self.style = 'filled'
            self.color = '.7 .3 1.0' #by HSB system
        elif self.speciesType == 'compound':
            if self.vizstyle == 'concise':
                self.beautified = self.beautifycpd()
            else:
                self.beautified = self.beautifycpd_desc()
            self.shape = 'plaintext'
            self.style = ''
            self.color = ''
        else:
            print "Node type error!"

    def beautifyec(self):
        """
        break up long EC strings; concise style
        """
        ecs = self.id.split(";")
        return "\\n".join(ecs)

    def beautifyec_desc(self):
        """
        break up long EC strings; descriptive style
        """
        ecs = self.id.split(";")
        outlist = []
        for ec in ecs:
            try:
                outlist.append( self.long2short(dict_ec_def[ec]) )
            except KeyError:
                    outlist.append(ec)
        return "\\n".join(outlist)

    def beautifycpd(self):
        """
        adjust multiple compound to multiple lines; concise style
        """
        return "\\n".join(self.label.split(";"))
        
    def beautifycpd_desc(self):
        """
        adjust multiple compound to multiple lines; descriptive style
        """
        cpds = self.id.split(";")
        s = []
        for c in cpds:
            try:
                s.append( self.long2short(dict_cpds_def[c].split(";")[0]) )
            except KeyError:
                s.append( c )
        return "\\n".join(s)
        
    def long2short(self, s):
        """
        adjust long compound/enzyme names to multiple lines.
        Find '-' in [20:35], in not, break at [30]
        """
        if len(s) < 31:
            return s
        else:
            news = ''
            middle = s[20:35]
            while len(middle) > 14:
                pos = max(middle.find(" "), middle.find("-"))
                if pos < 0: #nothing found
                    news += s[:30] + '\\ \\n'
                    s = s[30:]
                else:
                    news += s[:pos+20] + '\\ \\n'
                    s = s[pos+20:]
                middle = s[20:35]
            news += s
            return news
        
    def dotformat(self):
        buf_str = '[label=' + quote(self.beautified) + ', shape=' + quote(self.shape)
        if self.style:
            buf_str += ', style=' + quote(self.style)
        if self.color:
            buf_str += ', color=' + quote(self.color)
        return buf_str + ']'

class mnetwork(networkx.DiGraph):
    """
    This class expands on networkx.DiGraph 
    for modeling a bipartite metabolic pathway/network.
    Currently, input is expected from a SBML file.
    Both node and edge are extended by internal dicts.
    Visualization is emphasized, via manipulation of dot str.
    Edges can be undirected - in which case, stored as (cpd, ec, 0)
    """
    def neighbors(self, node):
        """
        In current networkx ver 1.0.dev, 
        DiGraph.neighbors() returns the same result as DiGraph.successors().
        This function overwrites stock function.
        """
        return self.predecessors(node) + self.successors(node)
    
    def read(self, infile, vizstyle='concise'):
        """
        Only SBML is implemented now.
        nodedict extends node by mnode class.
        edgedict extends edge by a tubple identifier.
        """
        self.vizstyle = vizstyle
        self.nodedict = {}
        self.edgedict = {}
        self.note = [] #a list of str to record editing info
        self.read_sbml(infile)

    def markedge(self, edge):
        """
        This is a provisional hack to assign a tuple identifier on edge.
        As mnetwork contains both directed and undirected edges,
        the stock edge (a tuple) can not handle the situation. 
        A better way is probably to reimplement edge as a new class.
        """
        wdict = {0: "undirected", 1: "directed"}
        d = self.get_edge_data(edge[0], edge[1])
        edgex = list(edge)
        if d == 0:
            edgex.sort()
        self.edgedict[edge] = (edgex[0], edgex[1], wdict[d])
 
    def getbowtie(self, node):
        """
        return the (edge-node-edge) set for a node with degree 2;
        edges could be 'to', 'from', 'un', return sth like
        set[ ('to', ec1), ('from', ec2) ]
        """
        upstream, downstream = self.predecessors(node), self.successors(node)
        bowties = []
        for ec in downstream:
            if self.get_edge_data(node, ec) > 0: #directed
                bowties.append( ('to', ec) )
            else:
                bowties.append( ('un', ec) )
        for ec in upstream:
            bowties.append( ('from', ec) )
        return set(bowties)
 
    def read_sbml(self, infile):
        """
        edge may not be unique; self.edgedict[edge] is unique.
        """
        print "Working on ", infile
        reader = libsbml.SBMLReader()
        sbm = reader.readSBML(infile)
        model = sbm.getModel()
        self.name = model.getName()
        # transfer model attr to mnetwork
        nodes = model.getListOfSpecies()
        tagdict = {'':'', }
        for n in nodes:
            if n.getName():
                m = mnode(n.getName(), n.getSpeciesType(), vizstyle=self.vizstyle)
                m.set_attr()
                self.nodedict[n.getName()] = m
                self.add_node(n.getName())
                tagdict[n.getId()] = n.getName()
        rxns = model.getListOfReactions()
        for r in rxns:
            ec, edges = self.parse_rxn(r, tagdict)
            if ec:
                self.add_node(ec)
            self.add_edges_from(edges)
        for edge in self.edges():
            self.markedge(edge)
        
        if '' in self.nodes():
            self.remove_node('')
        print "Got model - ", self.name
        print "number of nodes:", self.__len__()
        
    def parse_rxn(self, rxn, tagdict):
        """
        get ec and edges from a libsbml.Reaction .
        Undirected edges are marked by a weight = 0; directed by 1.
        Use tagdict to convert SBML ids to names.
        If no EC is given, use rxn ID for visualization purpose.
        """
        edges = []
        ec = tagdict[rxn.getModifier(0).getId()]
        if not ec:
            ec = rxn.getId()
            m = mnode(ec, 'enzyme', vizstyle=self.vizstyle)
            m.set_attr()
            self.nodedict[ec] = m
            self.add_node(ec)
        reactants = rxn.getListOfReactants()
        products = rxn.getListOfProducts()
        if ec:
            if products:
                for p in products:
                    edges.append((ec, tagdict[p.getId()], 1))
                for r in reactants:
                    edges.append((tagdict[r.getId()], ec, 1))
            else:
                for r in reactants:
                    edges.append((tagdict[r.getId()], ec, 0))
        # not using direct edges btw cpds; use rxn.id in place of ec
        # Backup plan:
        # for rxns without an enzyme, return (cpd, cpd) edges.
        else:
            for r in reactants:
                for p in products:
                    edges.append((tagdict[r.getId()], tagdict[p.getId()], 1))
        return ec, edges
    
    def thincopy(self):
        """
        return a copy without currency metabolites.
        Add a checkpoint - if a currency metabolite has only < 4 edges,
        keep it.
        """
        newgraph = self.copy()
        newgraph.nodedict = self.nodedict
        newgraph.edgedict = self.edgedict
        newgraph.vizstyle = self.vizstyle
        for c in currency:
            try:
                if len(self.neighbors(c)) > 3:
                    newgraph.remove_node(c)
            except networkx.exception.NetworkXError:
                pass
        return newgraph
        
    def modularize(self):
        pass


    def concentrate_cpds(self):
        """
        merge excessive cpds attached to only one ecstr,
        to be reimplemented 
        I.e., remove these nodes and merge to a new node
        """
        def add_mnode(cpdlist):
            newnode = ";".join(cpdlist)
            #
            # may be problematic - watch out for vizstyle
            #
            newnodeclass = mnode(id=newnode, type="compound", vizstyle=self.vizstyle)
            newnodeclass.set_attr()
            self.nodedict[newnode] = newnodeclass
            self.add_node(newnode)
            return newnode
        def merge(ec):
            """
            Mergy by ec; for cpd nodes with degree = 1 or 2,
            and >3 similar edges. Remvoe these nodes;
            create a composite node instead and fill in edges
            """
            neighbors = self.neighbors(ec)
          
            cpds1 = [c for c in neighbors if self.degree(c) == 1]
            to_bin, from_bin, un_bin = [], [], []
            for c in cpds1:
                if self.has_edge(c, ec):
                    # can be directed or undirected
                    if self.get_edge_data(c, ec) > 0: #directed
                        to_bin.append(c)
                    else:
                        un_bin.append(c)
                else:
                    from_bin.append(c)
            if len(to_bin) > 3:
                newnode = add_mnode(to_bin)
                self.add_edge(newnode, ec, 1)
                self.markedge((newnode, ec))
                self.remove_nodes_from(to_bin)
            if len(from_bin) > 3:
                newnode = add_mnode(from_bin)
                self.add_edge(ec, newnode, 1)
                self.markedge((ec, newnode))
                self.remove_nodes_from(from_bin)
            if len(un_bin) > 3:
                newnode = add_mnode(un_bin)
                self.add_edge(newnode, ec, 0)
                self.markedge((newnode, ec))
                self.remove_nodes_from(un_bin)

            cpds2 = [c for c in neighbors if self.degree(c) == 2]
            cpds2bowties = [(c, self.getbowtie(c)) for c in cpds2]
            bowties = [x[1] for x in cpds2bowties]
            hubs = findhubs(bowties)
            for h in hubs:
                # set[ ('to', ec1), ('from', ec2) ]
                cpdlist = [x[0] for x in cpds2bowties if h == x[1]]
                newnode = add_mnode(cpdlist)
                for e in h:
                    if e[0] == 'un':
                        self.add_edge(newnode, e[1], 0)
                        self.markedge((newnode, e[1]))
                    elif e[0] == 'to':
                        self.add_edge(newnode, e[1], 1)
                        self.markedge((newnode, e[1]))
                    else: #'from'
                        self.add_edge(e[1], newnode, 1)
                        self.markedge((e[1], newnode))
                self.remove_nodes_from(cpdlist)
        ##
        candidate_ecs = [ec for ec in self.nodes() if
                         self.nodedict[ec].speciesType == 'enzyme'
                         and self.degree(ec) > 4]
        for ec in candidate_ecs:
            mergelist = merge(ec)

    def zoom(self, degreelimit=7):
        """
        breaking up the most connected cpds (such as currency metabolites);
        so that the overall visualization is less cluttered.
        """
        def move(node, shapenum):
            counter = 0
            edges = [x for x in self.edges() if node in x]
            for e in edges:
                counter += 1
                newnode = node + str(counter)
                #
                newnodeclass = mnode(id=newnode, type="compound", vizstyle=self.vizstyle)
                newnodeclass.label = node
                newnodeclass.beautified = newnodeclass.beautifycpd()
                newnodeclass.shape = shapeslist[shapenum]
                newnodeclass.style = 'filled'
                newnodeclass.color = 'yellow'
                self.nodedict[newnode] = newnodeclass
                self.add_node(newnode)
                #
                newedge = list(self.edgedict[e])
                newedge[newedge.index(node)] = newnode
                newedge = tuple(newedge)
                self.add_edge(newedge[0], newedge[1])
                self.edgedict[newedge[:2]] = newedge
            self.remove_node(node)

        zoomlist = self.zoomlist(degreelimit)
        shapenum = 0
        for node in zoomlist:
            move(node, shapenum)
            shapenum += 1


    def zoomlist(self, degreelimit=7):
        degrees = [(node, len(self.neighbors(node))) for node in self.nodes()
                   if self.nodedict[node].speciesType == 'compound']
        degrees.sort(cmp=lambda x,y: cmp(x[1], y[1]), reverse=True)
        m = mean([x[1] for x in degrees ])
        #return [x[0] for x in degrees if x[1] > degreelimit]
        # min(degreelimit, 4*m)][:8]
        return [x[0] for x in degrees[:3]]

    def sbml2png(self, infile, pngfile=""):
        self.read(infile)
        if self.__len__() > 13:
            self = self.thincopy()
            self.concentrate_cpds()
            #self.zoom()
        dotstr = self.write_dotstr()
        #self.write_dotfile(dotstr, infile[:-3]+'dot')
        self.draw_png(dotstr, pngfile)
        
    def draw_png(self, dotstr, pngfile=""):
        """
        concentrate = False, label_on = False, [...]
        """
        G=pgv.AGraph(dotstr)
        if not pngfile:
            pngfile = self.name + ".png"
        G.draw(pngfile, prog='dot')

    def write_dot_fromfile(self, infile, outfile, mark="MetaFishNet"):
        f = open(outfile, 'w')
        f.write(self.write_dotstr(infile, mark=mark))
        f.close()
    def write_dotfile(self, dotstr, outfile):
        f = open(outfile, 'w')
        f.write(dotstr)
        f.close()
        
    def write_dotstr(self, mark="MetaFishNet"):
        """
        create str for dot file
        """
        #self.concentrate_cpds()
        #if zoom_on:
        outstr = self.make_header(mark)
        outstr += self.write_nodes()
        wdict = {"undirected": "dir=none", "directed": ""}
        for edge in self.edges():
            outstr += self.write_dotline(edge, 
                                         wdict[self.edgedict[edge][2]])
        # will add modules?
        outstr += "    }\n"
        return outstr

    def make_header(self, mark="MetaFishNet"):
        # compress
        h = 'digraph G {\n     ratio = auto; size="36,60"; concentrate=true'
        h += ';\n fontsize=28; label="' + mark + ": " + self.name
        h += '";\n node [fontname=Helvetica, fontsize=18]; edge[style=bold];\n'
        return h

    def write_nodes(self):
        """
        write nodes into dot file, with proper formats.
        Leave out isolated nodes.
        """
        buf_str = ""
        for n in self.nodes():
            if n and self.neighbors(n):
                buf_str += '    "' + self.nodedict[
                              n].id +'" ' + self.nodedict[n].dotformat() + ';\n'
        return buf_str

    def write_dotline(self, (node1, node2), attr=""):
        """write a line in DOT file, node1 -> node2
        """
        return '    "' + node1 + '" -> "' + node2 + '" [' + attr + '];\n'


    def make_note(self):
        pass





#
#-------------------------------------------------------------------
#

if __name__ == "__main__":
    #
    # fish paths
    pass




