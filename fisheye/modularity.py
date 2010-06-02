"""
modularity.py
implementation of the modularity-finding algorithm from
Mark E.J. Newman (2006), PNAS 103:8577-8582
"""

from numpy import *
import re

testnetwork = 'karate/karate.gml'

class module:
    """
    A mudule is a list of nodes within a reference network.
    self.divide(nobj) returns delta_Q and [group1, group2].
    If the division is not successful, returned groups are empty.
    The reference network is passed in via argument nobj.
    """
    def __init__(self, nodes):
        """
        Keep track of nodes by their names (ids)
        because the indexes change from network matrix to subnet.
        node names -> node index via reference_network.crd
        """
        self.nodes = nodes
        self.num_nodes = len(nodes)
        # modularity matrix for this module/subnet, B^(g)
        self.modularity_matrix = []
        # s_i = 1 if node i belongs to group1, -1 if group2
        self.s = [0]*self.num_nodes
        self.delta_Q = 0
        self.to_be_divide = True

    def quick_divide(self, nobj):
        g1, g2 = self.eigen_divide(nobj)
        self.to_be_divide = False
        return self.delta_Q, [g1, g2]
        
        
    def slow_divide(self, nobj):
        """
        testing --
        Divide module first by eigenvector then fine-tuning.
        # testing
        #g1, g2 = ['1', '2', '9', '10', '15', '16', '21', '30', '31', '32', '33', '34'], ['3', '4', '5', '6', '7', '8', '19', '23', '24', '25', '26', '11', '12', '13', '14', '17', '18', '20', '22', '27', '28', '29']
        self.fetch_modularity_matrix(nobj)
        self.calculate_s(g1, g2)
        self.delta_Q = self.compute_delta_Q(nobj)
        """
        g1, g2 = self.eigen_divide(nobj)
        print "eigen_divide Q", self.delta_Q, '\n'
        if g1:
            [g1, g2] = self.fine_tune(g1, g2, nobj)
            print "after fine-tuning: ", self.delta_Q
        self.to_be_divide = False
        return self.delta_Q, [g1, g2]


    def eigen_divide(self, nobj):
        """
        divide module in two by leading eigenvector.
        This function reproduced the karate club data as given in 
        http://www-personal.umich.edu/~mejn/courses/2004/cscs535/problems5.pdf
        """
        self.fetch_modularity_matrix(nobj)
        [eigenvalues, eigenvectors] = linalg.eig(self.modularity_matrix)
        group1, group2 = [], []
        if eigenvalues.max() > 0:
            # vector corresponding to max eigenvalue, normalized
            lead_vector = eigenvectors[:, eigenvalues.argmax()]
            for ii in range(self.num_nodes):
                if lead_vector[ii] < 0:
                    group1.append(self.nodes[ii])
                else:
                    group2.append(self.nodes[ii])
            self.calculate_s(group1, group2)
            self.delta_Q = self.compute_delta_Q(nobj)
        return group1, group2

    def fine_tune(self, group1, group2, nobj):
        """
        check the gain from shuffling a single vertex,
        as by the Kernighan-Lin method. Return best groups.
        """
        #
        #
        # debugging ----------
        #
        #best_groups = [group1, group2]
        #best_dQ = self.delta_Q
        
        def shuffle(g1, g2, node):
            # have to copy g1, g2 instead of assign pointers
            gg1, gg2 = [x for x in g1], [x for x in g2]
            if node in gg1:
                gg2.append(node)
                gg1.remove(node)
            elif node in gg2:
                gg1.append(node)
                gg2.remove(node)
            else:
                raise Exception, "Node not in either group."
            return gg1, gg2
        def onepass(g1, g2, nlist):
            """
            shuffle one node btw g1 and g2,
            return the best node and state
            """
            nn, qq, gg = '', -1, [ [], [] ]
            while nlist:
                node = nlist.pop()
                ggg1, ggg2 = shuffle(g1, g2, node)
                self.calculate_s(ggg1, ggg2)
                tmp_Q = self.compute_delta_Q(nobj)
                if tmp_Q > qq:
                    nn, qq = node, tmp_Q
                    gg = [[x for x in ggg1], [x for x in ggg2]]
            return nn, qq, gg
        def oneround(g1, g2):
            bestqq, bestgg, nn = 0, [], ''
            nlist = g1 + g2 + [nn]
            gg = [[x for x in g1], [x for x in g2]]
            while nlist:
                nlist.remove(nn)
                passlist = [x for x in nlist]
                nn, qq, gg = onepass(gg[0], gg[1], passlist)
                #print "oneround .... nn=", nn, " ... gg[0] len", len(gg[0]), "... Q=", qq
                #print gg[0]
                if qq > bestqq:
                    bestqq, bestgg = qq, gg
            return bestqq, bestgg

        best_dQ, best_gg = oneround(group1, group2)
        counter = 0
        while best_dQ > self.delta_Q:
            self.delta_Q, best_groups = best_dQ, best_gg
            counter += 1
            print "round...", counter, "Q=", best_dQ
            best_dQ, best_gg = oneround(best_groups[0], best_groups[1])

        return best_groups


    def first_Q(self, nobj):
        """
        This is implementation of Eq.[2]. Not really used because 
        it is the same  as self.compute_delta_Q at the first network division.
        """
        s = array(self.s)
        return dot(dot(s.transpose(), nobj.modularity_matrix), s
                                            )/(4.0*nobj.num_edges)

    def compute_delta_Q(self, nobj):
        """
        compute gain in modularity, delta_Q, by Eq.[5].
        In Newman's paper, this step is suggested to be bypassed by Eq.[7].
        However, Eq.[5] can be written as (1/4m) * SUM[(s_i*s_j - 1)*B_(i,j)],
        and using B and B.sum() is more convinient in Numpy.
        I have not evaluated how my approach affects speed of the algorithm,
        though matrix operations in Numpy are respectable.
        """
        dQ = 0
        for ii in range(self.num_nodes):
            for jj in range(self.num_nodes):
                dQ += (self.s[ii]*self.s[jj] - 1) * self.modularity_matrix[ii, jj]
        return dQ/(4.0*nobj.num_edges)


    def calculate_s(self, g1, g2):
        for ii in range(self.num_nodes):
            if self.nodes[ii] in g1:
                self.s[ii] = 1
            else:
                self.s[ii] = -1

    def fetch_modularity_matrix(self, nobj):
        am = zeros( (self.num_nodes, self.num_nodes) )
        for ii in range(self.num_nodes):
            for jj in range(self.num_nodes):
                am[ii][jj] = nobj.modularity_matrix[nobj.crd[self.nodes[ii]],
                                                    nobj.crd[self.nodes[jj]]]
        self.modularity_matrix = am


class network:
    """
    A network is represented as a list of nodes and a list of edges,
    which can then be translated into an adjacency matrix.
    To stick to numpy consistency, all data are in "array" type.
    Please also note that Python indexing starts from 0.
    To use this class, self.nodes and self.edges should be obtained first,
    either from self.read_gml_network(infile) or elsewhere. 
    Then self.prime_network() fills in degrees and matrices etc.
    Finally, self.module_analyze() finds all modules in this network.
    """
    def __init__(self):
        """
        Network attributes, predefined for code clarity.
        Matrix in 'array' type (see class docstr above).
        node = vertex.
        An edge from a to b is written as (a, b).
        """
        self.desc = ''
        self.adjacency_matrix = []
        self.modularity_matrix = []
        self.nodes = []
        self.edges = []
        self.num_nodes = 0
        self.num_edges = 0
        self.degrees = []   # same sequence as in self.nodes
        self.modules = []   # module instances
        self.Q = 0          # modularity value, 0 for undivided network
        self.crd = {}       # coordinate for nodes


    def specsplit(self):
        """
        produce modules by eigenvector method alone
        """
        self.prime_network()
        m = module(self.nodes)
        self.modules = [m]
        divisible = [x for x in self.modules if x.to_be_divide]
        while divisible:
            for mod in divisible:
                #print "working on module with", mod.num_nodes
                delta_Q, newgroups = mod.quick_divide(self)
                # check if division successful
                if delta_Q > 0 and newgroups[0]:
                    self.Q += delta_Q
                    self.modules.remove(mod)
                    for g in newgroups:
                        self.modules.append(module(g))
                    
                # else mod is final
            divisible = [x for x in self.modules if x.to_be_divide]
        outstr = "Q=" + str(self.Q) + "\n"
        for m in self.modules:
            outstr += ",".join(m.nodes) + "\n"
        return outstr


    def module_analyze(self):
        """
        function to combine eigenvector based splitting and fine-tuning,
        still testing
        """
        self.prime_network()
        m = module(self.nodes)
        self.modules = [m]
        divisible = [x for x in self.modules if x.to_be_divide]
        """
        while divisible:
            for mod in divisible:
                print "working on module", mod.nodes
                delta_Q, newgroups = mod.divide(
                )
                self.Q += delta_Q
                self.modules.remove(mod)
                for g in newgroups:
                    self.modules += module(g)
            divisible = [x for x in self.modules if x.divisible]
        print self.Q, self.modules
        """
        ### testing
        
        print m.divide(self)
        



    #
    # house-keeping functions
    #

    def prime_network(self):
        """
        prepare network for analysis, starting from nodes and edges
        """
        if self.nodes and self.edges:
            self.num_edges = len(self.edges)
            self.num_nodes = len(self.nodes)
            self.make_node_index()
            self.make_adjacency_matrix()
            self.compute_degrees()
            self.make_modularity_matrix()
        else:
            raise Exception, "Network nodes or edges not properly defined!"

    def make_modularity_matrix(self):
        """
        Modularity matrix is defined in Eq.[3],
        B_(i,j) = A_(i,j) - (k_i*k_j)/2m.
        """
        am = zeros( (self.num_nodes, self.num_nodes) )
        for ii in range(self.num_nodes):
            for jj in range(self.num_nodes):
                am[ii, jj] = self.adjacency_matrix[ii, jj] - \
                         self.degrees[ii]*self.degrees[jj]/(2.0*self.num_edges)
        self.modularity_matrix = am
        
    def compute_degrees(self):
        """
        Compute degree for each node from adjacency matrix
        """
        if not self.degrees:
            for ii in range(self.num_nodes):
                self.degrees.append(self.adjacency_matrix[ii].sum())
        else:
            print "Degrees already exist!"

    def make_adjacency_matrix(self):
        """
        Make adjacency matrix from nodes and edges.
        Edges are treated as not directional, thus matrix is symmetrical.
        Multiple edges between two nodes are allowed.
        """
        am = zeros( (self.num_nodes, self.num_nodes) )
        for edge in self.edges:
            am[self.nodes.index(edge[0]), self.nodes.index(edge[1])] += 1
            # take out the next line will make a directional network
            am[self.nodes.index(edge[1]), self.nodes.index(edge[0])] += 1
        self.adjacency_matrix = am

    def make_node_index(self):
        """
        mapping node name/id to original index
        """
        for ii in range(self.num_nodes):
            self.crd[self.nodes[ii]] = ii

    #
    # The read_gml functions read network from GML format.
    # Alternatively, self.nodes and self.edges can be directly defined.
    #

    def read_gml_network(self, infile):
        (self.nodes, self.edges) = self.read_gml(infile)

    def read_gml(self, infile):
        """
        read a GML file 
        (http://www.infosun.fim.uni-passau.de/Graphlet/GML/gml-tr.html),
        such as the karate.gml used in this demo,
        and return lists of nodes and edges. Based on regular expressions.
        This is not a bullet-proof parser of GML. But if you have looked at 
        your file and will check the result, this function is likely to get 
        it right. Only ids are used; labels are ignored.
        An alternative Python GML parser can be found 
        in NetworkX (networkx.lanl.gov).
        """
        def scan_node(str_block):
            # input example: 'node   [     id 5   ]'
            a = str_block.split()
            #return the id after 'id'
            return a[a.index('id') + 1]
        def scan_edge(str_block):
            # input example: 'edge   [     source 4     target 2   ]'
            a = str_block.split()
            sid, tid = a[a.index('source') + 1], a[a.index('target') + 1]
            return (sid, tid)
        s = open(infile).read()
        # remove line breaks to ease re.findall
        s = " ".join(s.splitlines())
        nodes = re.findall('node.*?\\[.*?id.*?\\]', s)
        edges = re.findall('edge.*?\\[.*?source.*?\\]', s)
        print "Found %d nodes and %d edges." %(len(nodes), len(edges))
        nodes = [scan_node(x) for x in nodes]
        edges = [scan_edge(x) for x in edges]
        return (nodes, edges)





if __name__ == "__main__":
    print "\n\n\n\n\n\n"
    net = network()
    net.read_gml_network(testnetwork)
    #net.module_analyze()
    net.specsplit()
    
    #print m.adjacency_matrix.sum()
    #print m.modularity_matrix.sum()
    #m.eigen_divide()
    


