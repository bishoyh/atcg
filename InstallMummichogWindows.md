# mummichog installation on MS Windows

## Installation of mummichog on MS Windows ##

### Overview ###

Please refer to the main manual if you look for the installation on Linux/Mac OS. This page describes the installation on MS Windows, which will require some level of computer skills. You will need Python 2.x, numpy, scipy and networkx. Their specific versions may vary.


### Step-by-step ###

1) Download Python: python-2.7.2.msi from http://python.org. (This guide uses 32-bit packages throughout.) Install Python by double-clicking the installer and following its instructions.

2) Download and install numpy: numpy-1.6.1-win32-superpack-python2.7.exe from http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/ . Install numpy by double-clicking the installer and following its instructions.

3) Download and install scipy: scipy-0.10.1-win32-superpack-python2.7.exe from http://sourceforge.net/projects/scipy/files/scipy/0.10.1/ . Install scipy by double-clicking the installer and following its instructions.

4) Now verify your Python installation. Open a command line window, type "python" and hit "ENTER". This should bring up the Python interactive mode. If "python" is not found, add it to your system path by doing the following:
right click "My Computer" or via control panel, goto System Properties\Advanced\Environment Variable (see image below), edit "Path" to add "C:\Python27", depending on your installed location. Save change and reboot the computer.

![http://atcg.googlecode.com/files/set_win_path2.png](http://atcg.googlecode.com/files/set_win_path2.png)

5) Download networkx: networkx-1.5.zip from http://networkx.lanl.gov/download/networkx/ . Unzip the file. In a command line window, change to the directory of your unzipped networkx, run
```
python setup.py install
```
This should install networkx to your Python libraries.

**Skip to Step 14) for mummichog version 0.9+, which no longer requires Graphviz/pygraphviz. Steps 6-13) are ONLY for version 0.8 and older.**


---


Next, it will take some effort to install the visualization components, Graphviz and pygraphviz. You will need a C compiler, and make sure pygraphviz finds Graphviz. Administrator's privilege is required for at least one of the packages (Graphviz).

6) Download a MinGW installer: mingw-get-inst-20110802.exe from http://sourceforge.net/projects/mingw/files/Installer/mingw-get-inst/ . Please NOTE that I used version 20110802, because the newer versions have compatibility issues (in gcc-4.6). Install MinGW by double-clicking the installer and following the instructions. You only need the C compiler for this current purpose.

7) Download Graphviz: graphviz-2.28.0.msi from http://www.graphviz.org/Download_windows.php . Install Graphviz by double-clicking the installer and following its instructions. It may require administrator's privilege.

8) Now verify the programs are in your system path. Find the "Path" entry in the same way as described in step 4. It's a long string - you should have segments like "C:\Python27;C:\MinGW\bin;C:\Program Files\Graphviz 2.28\bin", depending on the installed locations. If not, edit "Path" and add them. Save change and reboot the computer.

9) Here begins the fun part with pygraphviz (a discussion on pygraphviz on Windows can be seen at http://stackoverflow.com/questions/2798858/installing-pygraphviz-on-windows-python-2-6).
Download pygraphviz-1.1.tar.gz from http://networkx.lanl.gov/download/pygraphviz/ .
Unpack the pygraphviz source package, go to its directory. Edit the setup.py file to define (**depending on your Graphviz location**):
```
library_path=r"C:\Program Files\Graphviz 2.28\lib\release\lib"
include_path=r"C:\Program Files\Graphviz 2.28\include\graphviz"
```
Save the file.

10) In a command line window, go to your unpacked pygraphviz directory and run
```
python setup.py build -c mingw32 install
```
An error may occur like this
```
  ...........
  File "C:\Python27\lib\distutils\ccompiler.py", line 1121, in gen_lib_options
    opt = compiler.runtime_library_dir_option(dir)
  File "C:\Python27\lib\distutils\unixccompiler.py", line 285, in runtime_librar
y_dir_option
    compiler = os.path.basename(sysconfig.get_config_var("CC"))
  File "C:\Python27\lib\ntpath.py", line 198, in basename
    return split(p)[1]
  File "C:\Python27\lib\ntpath.py", line 170, in split
    d, p = splitdrive(p)
  File "C:\Python27\lib\ntpath.py", line 125, in splitdrive
    if p[1:2] == ':':
TypeError: 'NoneType' object is not subscriptable
```

11) The error came from a bug in Python distutils library. Edit your "C:\Python27\lib\distutils\unixccompiler.py" file (line 285, in runtime\_library\_dir\_option), change
```
    compiler = os.path.basename(sysconfig.get_config_var("CC"))
```
to
```
    compiler = 'gcc'
```

12) Rerun `python setup.py build -c mingw32 install` to compile and install pygraphviz. (You should restore your "C:\Python27\lib\distutils\unixccompiler.py" file afterward!)

13) If you have got this far, congratulations! Now let's fix a bug in pygraphviz. Find your pygraphviz installation (it should be something like C:\Python27\Lib\site-packages\pygraphviz\", and edit file "agraph.py". In function "_get\_prog", around line 1237, change
```
        return runprog
```
to
```
        return '"' + runprog + '"'
```
This is to double quote the Graphviz command, so that it will be called correctly even if space characters are in its path._


---


14) Download the latest version of mummichog from http://code.google.com/p/atcg/downloads/list . Unpack it to wherever you like. In a command line window, go to the unpacked mummichog directory and test your installation:
```
cd test
python ..\mummichog\main.py -i test.sig -r test.ref -o mytest
```


### Sample installation files ###
Here are the list of files I tested on WinXP and Windows 7 (in the order of installation):
```
python-2.7.2.msi
numpy-1.6.1-win32-superpack-python2.7.exe
scipy-0.10.1-win32-superpack-python2.7.exe
networkx-1.5.zip
mingw-get-inst-20110802.exe  
graphviz-2.28.0.msi          
[add PATH; reboot computer]
pygraphviz-1.1.tar.gz
mummichog-0.5.16.zip (please use the latest version)
```




Going back to mummichog manual:
http://code.google.com/p/atcg/wiki/mummichog_for_metabolomics