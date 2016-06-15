#!/bin/bash
# Installation Guide for VisTrails on Arch Linux

# Seems like VisTrails don't work with python 3 (got it working on python 2.7), the error message with python 3 is 'File "vistrails/run.py", line 146 except SystemExit, e:'. The most complete list of dependencies seem to be: http://www.vistrails.org/index.php/Building_From_Source

# The python2 package manager (pip2) was installed to simplify the installation of the python library 'usagestats' (and its dependency 'requests'). The VisTrails didn't run without this python library.

pacman -S python2-pip # will install python2
pip2 install requests usagestats

# After this you can run vistrails as root and let it download the dependencies (I don't tested if this work), or install the dependencies by pacman (below).

pacman -S python2-pyqt4

# With this, vistrails can open the GUI, and will prompt you if you want to install additional packages. It will use pip to do that, so there probably no difference between installing them by hand, or letting vistrails do the work. But for me, there was a problem when I tried executing vistrails as root (dialog was completely gray). So I installed the following packages by pip2.

pip2 install scipy iphyton certifi scikit-learn file_archive xlrd sqlalchemy matplotlib tej backports.ssl_match_hostname

# VTK is needed for the examples on the vistrails guide (https://sourceforge.net/projects/vistrails/files/vistrails/v2.2.4/VisTrails.pdf/download). Fell free to remove this line if you will not follow the examples on the vistrails guide.
#pacman -S vtk
# NOTE: was not possible to get vtk working even installing the vtk or vtk6 package. Will search for more information.

