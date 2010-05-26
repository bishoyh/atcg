""" batch_convert.py,
convert a batch of KEGG KGML to SBML
use k2sb

For zebrafish metabolic network, use path# < 1200
"""

from k2sb import *
import os

dre_path = 'drexml/'

files = os.listdir(dre_path)
c = K2SB()

for file in files:
    if file[:3]=='dre' and int(file.split('.')[0][-4:]) < 1200:
        #if 'sb_'+file not in files:
        print file
        c.convert(dre_path+file, dre_path+'sb_'+file)



