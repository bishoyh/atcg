#how to programmatically extract network model from KEGG

# Extract network model from KEGG database #

KEGG is a great resource. However, it is difficult to automatically extract a network model from KEGG data because they tend to mix visual elements with true network structure. My program does it by combining KGML files and database API.

See also http://atcg.appspot.com/3

Code available in "Source", trunk/utility/kegg/. It demonstrates 1) how to fetch KGML files by ftp 2) extract model and write SBML 3) in a batch mode.

Another program, KEGG2SBML (http://sbml.org/Software/KEGG2SBML) takes a similar strategy. Its newer version was not available by the time I wrote my code. Still, my code may be more accessible: it's under 100 lines of Python versus 3000 lines of Perl.



### update ###

As KEGG has newer organization of their FTP site and some of the pathways, the retrieval point in the scripts should be changed accordingly to e.g.
` ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/organisms/dre/ `

At least for zebrafish, all KGML files under number 1100 are conventional metabolic pathways. File dre01100.xml is a total map and file dre04070.xml is for a signaling pathway.