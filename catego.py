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
identify go_sleek as one of the categories:
    
    GO:0008150
    biological_process
    alt_id: GO:0000004
    alt_id: GO:0007582
    
    GO:0005575
    cellular component
    alt_id: GO:0008372
    
    GO:0003674
    molecular function
    alt_id: GO:0005554

"""

import MySQLdb
from go_sleek_list import go_sleek

biopro = ['GO:0008150', 'GO:0000004', 'GO:0007582']
celcom = ['GO:0005575', 'GO:0008372']
molfun = ['GO:0003674', 'GO:0005554']

def overlap(list1, list2):
    return set(list1).intersection(set(list2))

def go_analyze(cursor):
    
    godict = {}
    # get all ancestors
    query_templ = "SELECT DISTINCT ancestor.acc FROM term \
        INNER JOIN graph_path ON (term.id=graph_path.term2_id)  \
        INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id) \
        WHERE term.acc = '%s'"
    
    for item in go_sleek:
        godict[item] = []
        ancestors = []
        cursor.execute(query_templ % item)
        result = cursor.fetchall()
        for x in result:
            if x[0] and ":" in x[0] and x[0] != item:
                ancestors.append(x[0])

        if overlap(ancestors, biopro):
            godict[item].append('biological process')
        if overlap(ancestors, celcom):
            godict[item].append('cellular component')
        if overlap(ancestors, molfun):
            godict[item].append('molecular function')

    return godict

if __name__ == "__main__":
    LOCAL_GO_DB = 'go1123'
    dbc = MySQLdb.connect(user='guest', db=LOCAL_GO_DB)
    cursor = dbc.cursor()
    godict = go_analyze(cursor)
    
    f = open('cate_gosleek.txt', 'w')
    for key, value in godict.items():
        f.write(key + '\t' + ','.join(value) + '\n')
    f.close()









