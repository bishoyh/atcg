#mummichog, modular and pathway analysis of metabolomics
# mummichog: predicting network activity for high throughput metabolomics #

FOR VERSION 1.0.x

mummichog is a Python program for analyzing data from high throughput, untargeted metabolomics. It leverages the organization of metabolic networks to predict functional activity directly from feature tables, bypassing metabolite identification. The features include

  1. computing significantly enriched metabolic pathways;
  1. identifying significant modules in the metabolic network;
  1. visualization of top networks in web browser
  1. visualization that also plugs into Cytoscape
  1. tentative annotations
  1. metabolic models for different species through plugins

For older versions, please refer to
https://code.google.com/p/atcg/wiki/MummichogPreviousVersion

## Version 1 is here ##
As Google disallows package download here, newer versions will be linked to our own server, at
<a href='http://ps270548.dreamhost.com/init/downloads/mummichog-1.0.5.zip'>mummichog 1.0.5</a>.
Major updates include

  1. improved network visualization in browser
  1. added support for single input file, where threshold of input significant list is automatically determined
  1. output files in Cytoscape 3 and Excel formats



---

## Install ##

  * Download latest version of mummichog.
  * Mummichog is a Python 2.x package, that requires a few libraries (numpy and scipy, networkx, xlsxwriter). While these dependency libraries are common and easy to install, we recommend to install Anaconda (a scientific Python distribution, http://continuum.io/downloads), which will fulfill all dependencies.

  * You can unpack mummichog to where you like. No installation is necessary for mummichog. Go to the test/ directory to test mummichog on your computer:
```
python ../mummichog/main.py -f testdata.txt -o myoutput
```
This will produce a bunch of messages in the console and a new directory containing analysis reports.


---

## Input format ##

  * A tab-delimited text file is used as input. One feature per line. Any rows starting with '#' will be skipped. Each line has to contain the following tab-delimited numeric fields:

> m/z, retention time, p-value, statistic score

Example of input file:
```
mz      rtime   p-value t-score
186.0185697     463.4557469     0.000149751400132       3.82298721075
279.1773473     90.0607054      0.000399613326314       3.56621327128
344.1330624     124.1798842     0.000998323061251       -3.31200757387
215.9641894     132.6805079     0.00105418285794        -3.2963838682
177.0323244     77.31339162     0.00121065359218        3.25639199833
296.0973768     135.3519579     0.00171645907855        -3.15367022007
527.3784209     593.6222016     0.00176815004959        -3.14481340469
```

  * Older versions use a two-file input format, one file for significant list and the other for reference list. This format is still supported, and can be used to specify significant list if one wishes to avoid automated selection in the one-file format.



---

## Usage ##

mummichog can be copied as a standalone program, no installation is necessary if the library dependence is satisfied. Each run will generate a time-stamped directory containing reports and graphs. The program is called either in relative path or absolute path. It is strongly recommended that you run mummichog from your data directory not the program directory. When evoked without argument, mummichog displays options:

```
skylee@rhino:~/test$ python ../mummichog/main.py

    Usage example:
    python main.py -f mydata.txt -o myoutput
    
        -i, --input: file of selected significant m/z features, one per line
        -r, --reference: file containing all m/z features
        -f, --infile: single file as input, 
              containing all features with tab-delimited columns
              m/z, retention time, p-value, statistic score
        
        -n, --network: network model to use (default human_mfn), 
              [human, human_mfn, mouse, fly, yeast]
        
        -o, --output: output file identification string (default 'mcgresult')
        -k, --workdir: directory for all data files.
              Default is current directory.
        
        -m, --mode: analytical mode of mass spec, [positive, negative, dpj].
              Default is dpj, a short version of positive.
        -u, --instrument: [5, 10, 25, FTMS, ORBITRAP].
              Any integer is treated as ppm. Default is 10. 
              Instrument specific functions may be implemented.
              
        -p, --permutation: number of permutation to estimate null distributions.
              Default is 100.
        -z,   --force_primary_ion: M+H[+] (M-H[-] for negative mode) must be 
              present for a predicted metabolite, [True, False].
              Default is False.
              
        -d, --modeling: modeling permutation data, [no, gamma].
              Default is gamma.
        -e, --evidence: cutoff score for metabolite to be in activity network.
              Default is 3.

```

An example run will yield screen messages below:

```
skylee@rhino:~/test$ python ../mummichog/main.py -f testdata.txt -o rhino9t -e 6

[logo]

mummichog version 0.9.2 r44:244a00cb7bbc 

Started @ Tue Jun 11 09:51:20 2013

Loading metabolic network MFN_1.10.2...
cpds with MW: 2016
Got 397 significant features from 9116 references

Pathway Analysis...
query_set_size = 210 compounds
total_feature_num = 1445 compounds
Resampling, 100 permutations to estimate background ...
 1 2 3 4... 100
Pathway background is estimated on 11900 random pathway values

Modular Analysis, using 100 permutations ...
 1 2 3 4... 99 100
Null distribution is estimated on 1227 random modules
User data yield 14 network modules

Got ActivityNetwork of 46 metabolites.

Annotation was written to test/1370958677.93.rhino9t/csv/_tentative_featurematch_rhino9t.csv

Pathway analysis report was written to 
test/1370958677.93.rhino9t/csv/mcg_pathwayanalysis_rhino9t (.csv and .xlsx)

Modular analysis report was written to
test/1370958677.93.rhino9t/csv/mcg_modularanalysis_rhino9t (.csv and .xlsx)

Inspected network report was written to
test/1370958677.93.rhino9t/csv/InspectedNodes_ActivityNetwork (.csv and .xlsx)

Exporting top modules to test/1370958677.93.rhino9t/sif/...

HTML report was written to
test/1370958677.93.rhino9t/result.html

Finished @ Tue Jun 11 09:52:23 2013

```




---

## Output ##

A default run will generate a directory structure like this

```
1370958677.93.rhino9t/
    result.html
    mummichog.log
    tsv/
        InspectedNodes_ActivityNetwork.csv
        mcg_pathwayanalysis_rhino9t.csv
        mcg_pathwayanalysis_rhino9t.xlsx
        mcg_modularanalysis_rhino9t.csv
        mcg_modularanalysis_rhino9t.xlsx
        _tentative_featurematch_rhino9t.csv
        _tentative_featurematch_rhino9t.xlsx
    sif/
        ActivityNetwork.txt
        module_1.txt
        ...
    web/
        ...
```

A summary report is given in "result.html", which can be viewed in modern web browsers (excluding IE 8 or older). Internet connection is required to link a visualization library.

Full tables, results from annotation, pathway analysis and network module analysis, are given under the tsv/ directory. Results are provided in two formats: tab-delimited (.tsv) and Excel (.xlsx).


## Visualization ##

**Web browser based visual**
The "result.html" visualizes the activity network and up to top 5 modules, through javascript based technologies, while software development continues.

**Cytoscape and .txt**
Network/pathway can be described in .txt files. One can use Cytoscape (cytoscape.org) to visualize the .txt files. Cytoscape is a powerful tool to work on network graphs in a friendly graphic interface. Please refer to Cytoscape's guides for details.

**Graphviz and DOT (deprecated)**
Older versions of mummichog used Graphviz/DOT based visualization. This is deprecated now, but still in the source code so that curious developers may tweak it.



---

## Q & A ##

**I ran the program. Now what?**

Open the "result.html" in the output directory using a web browser for a summary report. More data are under the tsv/ directory.


**What's the difference between Pathway analysis and Modular analysis?**

Pathways are predefined units with human knowledge (hence possibly more sensitive). Analysis via network structure (modular) is less biased. A module can be within a pathway or in between several pathways. These two approaches are rather complementary.


**How are m/z features mapped to metabolites?**

Mummichog considers about 20 derivatives, e.g. M+H[+] and M(C13)+2H[2+], collected from literature and expert polling.
Matching stringency can be a single value, e.g. 10 ppm, or a function that is determined on the instrument.
Comparing an input m/z to data in our model will give 0 to several tentative matches to metabolites.
(The likelihood gets greater when pathway and module information is considered, the underlying idea of mummichog.)
There is also an option "-z", which enforces the presence of M+H[+] or M-H[-] in predicted metabolites.


**How are the p-values computed in the Pathway analysis?**

The relationship between m/z features and compounds is many-to-many.
The pathway enrichment test is done in compound space not in m/z
space.

For example in the test above, the list of selected significant m/z features converts to 117 compounds; the reference list converts to 1414 compounds (background universe). If we have 2 hits on a pathway of 7 metabolites (present in the reference; the pathway may be larger), Fisher exact test computes on a contingency table from these numbers, and produces a right-tail p-value = 0.108. If we had 3 hits on the pathway, the p-value = 0.015, indicating a more significant enrichment.

If and only if the number of matched compounds is greater than
the number of corresponding m/z features in a specific pathway, the
latter is used in Fisher exact test in order not to inflate the
significance. This is a practical workaround as we don't have a good
estimation of the underlying statistical structure.

Fisher exact test, however, only gives the enrichment significance without taking into account of the probability of mapping user's input m/z features to pathways. We amend this by resampling from the reference list in a similar way as Berriz et al 2003 (Bioinformatics 19:2502). The p-values from resampling are modeled as a Gamma distribution, and the EASE scores (Hosack et al. Genome Biology 2003, 4:[R70](https://code.google.com/p/atcg/source/detail?r=70)) of real data are converted to adjusted p-values on the CDF (cumulative distribution function) of the Gamma model. The EASE score is a variant of Fisher exact test by removing one hit in the pathway analysis. This will penalize pathways with fewer hits, thus fits the rationale of mummichog, "strength from groups".


**How are the p-values computed in the Modular analysis?**

From the example input list of significant m/z features, we get 117 tentatively matched compounds. These compounds are looked up in our human metabolic network model, and modules are identified by their connectivity. An activity score is calculated per module, based on both enrichment and connectivity.

Similarly, we can pick up random m/z features from our reference list, and repeat the whole process to calculate scores for these random modules. Repeated sampling generates a null distribution, indicating the background activity of our data. Module p-values are calculated on this random background, again modeled as a Gamma distribution as for pathway analysis. Since permutation is done on the reference data from the same experiment, analytical bias is largely avoided.


**How many permutations should I do with the -p argument?**

For Modular analysis, the screen message will report the number of random modules used to estimate null distribution. This number is much greater than -p argument, because each permutation generates many subnetworks (modules).
Similarly for Pathway analysis, each permutation will go through ~120 pathways.

The default number of 100 for -p is usually enough to get a good estimation. Of course, the more permutations, the more accurate your p-value gets (and the longer the program runs).


**Why are there pathway names in my modular analysis?**

They are there just to show what pathways the metabolites in a module belong to. This does NOT mean a statistical significance of any of these pathways, which is tested in the pathway analysis.

**What are the metabolic models behind mummichog?**

There is a human reference model integrated from KEGG, UCSD BiGG (Duarte et al. 2007) and Edinburgh Model (Ma et al. 2007). The integration process was described in
  * `Li S, Pozhitkov A, Ryan RA, Manning CS, Brown-Peterson N, Brouwer M. Constructing a fish metabolic network model. Genome Biology. 2010; 11(11):R115`.
  * Pathways in this model can be browsed at http://metafishnet.appspot.com/hbrowse .

Metabolic models from BioCyc (biocyc.org - human, mouse, fly and yeast) are now included. Support for other species will be added in the future.


**How does mummichog relate to metabolomics databases?**

The current metabolic models in the scientific community are based on the substantial biochemical research over the past century, and lately on genome sequences. Among these, KEGG is a leading force. Metabolomics has developed on the advanced analytical instruments, where much of the data are not yet connected to the preexisting metabolic models. This is the case for Metlin, HMDB and MMCD.

Mummichog tries to make the connection between these two disparate fields. However, a significant effort will be needed to incorporate more chemical data into metabolic models. Mummichog's top role is functional analysis. If one's primary goal is to match a spectral feature, the metabolomics databases will be the right place to go.

**What's the difference between 'dpj' and 'positive' in MS operation mode?**

Associated with each mode is a table to compute isotopic and adduct forms. 'positive' uses the generic table for positive mode, while 'dpj' uses a subset of 'positive' to better fit the protocol in our lab.

**How are currency metabolites dealt with?**

Currency metabolites are ubiquitous and not informative in analyzing networks structures.
The default setting is to exclude them from default model, which has little effect on the prediction results.
The list of currency metabolites in mummichog is ['H2O', 'H+', 'Oxygen', 'NADP+', 'NADPH', 'NAD+', 'NADH', 'ATP', 'Pyrophosphate', 'ADP', 'Orthophosphate', 'CO2'].


**How do I update to a new version?**

You can uppack the newer version and run it side by side with the older version. There's no conflict. Mummichog always records version numbers and parameters in the output.


**Is mummichog open source?**

Getting there but not quite yet. Source code is released, but we have to verify the packaged metabolic data (not the code per se) before we can truly open source the program.




---


Authors: Shuzhao Li, Andrei Todor