#FishEye is a visualization tool based on on Networkx and PyGraphviz.

Metabolic pathways are treated as bipartite networks. Many details of styling are
manipulated through mid-level markups. A number of optimizations are used to keep pathway graphs less cluttered. Networkx and Pygraphviz are described at: http://networkx.lanl.gov .


An example of FishEye:

![http://170.140.137.82/shuzhao/mfn1.9/fishthumbs/thumbnail_mfn1v9path161.png](http://170.140.137.82/shuzhao/mfn1.9/fishthumbs/thumbnail_mfn1v9path161.png)

The main input format is SBML. Main functions are in metabolicnet.py.
### Example of drawing a single pathway within Python console ###

```
>>> from metabolicnet import *
>>> mx = mnetwork()
>>> mx.read('mfnpath250.xml')
Got model -  MetaFishNet prov_path 30
number of nodes: 40
>>> newmx = mx.thincopy()
>>> newmx.concentrate_cpds()
>>> dotstr = newmx.write_dotstr()
>>> newmx.draw_png(dotstr)

```


### Example of using FishEye for making MetaFishNet graphs in batch ###

```

fishdir = 'fish/'
dotmark = 'MetaFishNet'

pathdir = 'pathways/'
sbmldir = 'sbml/'
dotdir = 'dots/'
pngdir = 'pngs/'
dotdir2 = 'extdots/'
pngdir2 = 'extpngs/'

import os, sys
sys.path.append('/home/habpi/projects/mfn/src/fisheye')

from path2sbml import *
from metabolicnet import *


#   convert paths to SBML
pathfiles = os.listdir(fishdir + pathdir)
for f in pathfiles:
    P = pathway()
    P.read_pathfile(fishdir + pathdir + f)
    P.write_sbml(fishdir + sbmldir + f[:-4] + 'xml')

#   visualize
sbfiles = os.listdir(fishdir + sbmldir)

for f in sbfiles:
    mx = mnetwork()
    mx.read(fishdir + sbmldir + f, 'concise')
    if mx.__len__() > 13:
        mx = mx.thincopy()
        mx.concentrate_cpds()
        if len(mx.edges()) > 360:
            mx.zoom()

    dotstr = mx.write_dotstr(dotmark)
    out = open(fishdir + dotdir + f[:-3] + 'dot', 'w')
    out.write(dotstr)
    out.close()
    
    mx.draw_png(dotstr, fishdir + pngdir + f[:-3] + 'png')


for f in sbfiles:
    mx = mnetwork()
    mx.read(fishdir + sbmldir + f, 'descriptive')
    
    if mx.__len__() > 13:
        mx = mx.thincopy()
        mx.concentrate_cpds()
        if len(mx.edges()) > 360:
            mx.zoom()
    
    dotstr = mx.write_dotstr(dotmark)
    out = open(fishdir + dotdir2 + f[:-3] + 'dot', 'w')
    out.write(dotstr)
    out.close()
    
    mx.draw_png(dotstr, fishdir + pngdir2 + f[:-3] + 'png')

```