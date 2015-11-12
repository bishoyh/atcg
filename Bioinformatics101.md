#How to start with bioinformatics data analysis

# Introduction #

The alternative title would be "how to move from Excel to command line". Bioinformatics is an important tool in many fields (of course there's a trendy term "data science"). Each field has a specific set of tools and work flows, which will not be covered here. This is meant to be a practical primer to get people started, then they can become experts on their own.


# Tools to use #

  * **UNIX Shell.** Shell, the common command line interface on UNIX-like systems (including Mac OS and Linux), is how one gets intimate with the computer. Graphic interface is great, but can be limiting because it is the road some developers paved for you. Shell is the AWD Jeep that enables you go off-road. Yes, there's Windows DOS. But scientific computing is traditionally centered on UNIX (among many other reasons).

  * **R**. R is a computing environment, heavy on statistics. A bioinformatician will have to use R at some point because some tools are only available there. But proficiency in R can be the only thing a scientist needs to do all data analysis.

  * **Python.** Python is the all-purpose, enterprise grade computer language that has been replacing Perl in the bioinformatics world.

  * **Data visualization.** One can get high-quality figures by using Python, R, or even MS Excel. I use MeV (tm4.org) and Gitools (gitools.org) for quick heat maps and a few things. Cytoscape (cytoscape.org) is the most popular network visualization tool. One can always use a good image editor (e.g. Adobe Photoshop or GIMP).


# Setting up a working computer #

R and Python run on all three major operating systems (OS): Windows, Mac OS and Linux. Well, almost. In reality, I had a lot more troubles with some R packages on Linux than on Windows. Installing some Python libraries is so much easier on Linux than on Windows. Although one can live with a single OS, a bioinformatician has to know Linux, even if only because most servers run on Linux or something similar. Mac OS has a very polished user interface, and is compatible with most Linux programs. The downside is that data exchange can be cumbersome (e.g. writing to NTFS can be a pain), and a Mac computer with high CPU/storage gets quite pricey.

Recommendation in our lab is to run Linux VM on Windows.
A VM (virtual machine) is a self-contained guest OS. The host OS here is MS Windows, for two reasons:
  * Some mass spectrometry software is Windows only. We can't get around that.
  * Compatibility is necessary to work with the rest of world.

I use Virtualbox (virtualbox.org) on my Windows laptop computer, and run Linux Mint 17 Mate as guest VM. The Linux VM is my main work environment.

A few rules of good practice:
  * Separation of data from tools;
  * Regular data back up at different locations;
  * Keeping notes of analyses (e.g. using README or log files).


# Data analysis using Python #

By installing the Anaconda distribution (http://continuum.io/downloads), one gets all common Python tools for scientific computing. Among them,
  * Scipy has many functionalities, including scipy.stats library.
  * Scikit-learn is a great collection of machine learning tools.
  * MatPlotLib started as a clone of Matlab style plotting tool, but has evolved into a powerful library with well polished APIs.
  * pandas is data analysis library that gets much of R flavors.

> ... to show example codes...



# Data analysis using R #

Install R and Bioconductor, and a few more depending on your field.



# The rest is on your own... almost #
  * One can probably become a master of common computational techniques by finishing this course (highly recommended):
https://software-carpentry.org/lessons.html

  * Get a "support group" if you can. There are also online places one can ask questions: http://stackoverflow.com; https://www.biostars.org.
