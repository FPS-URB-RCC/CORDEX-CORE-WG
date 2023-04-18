# PyNCplot
# Python to compute diagnositcs coded in Fortran from netCDF files.
# L. Fita, Centro de Investigaciones del Mar y la Atmosfera (CIMA), Argentina
# From L. Fita work in different places: CCRC (Australia), LMD (France)
# More information at: http://www.xn--llusfb-5va.cat/python/PyNCplot
#
# pyNCplot and its component nc_var.py comes with ABSOLUTELY NO WARRANTY. 
# This work is licendes under a Creative Commons 
#   Attribution-ShareAlike 4.0 International License (http://creativecommons.org/licenses/by-sa/4.0)
#
## Python script to indepdententy test different components 
## L. Fita, CIMA. 
####### ####### ##### ##### #### ### ## #

from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import os
import re
import subprocess as sub
import numpy.ma as ma

# Importing generic tools file 'generic_tools.py' and the others
import gen_tools as gen

# Importing Fortran modules and the others
import module_ForDef as fdef
import module_ForBas as fbas
import module_ForGen as fgen
import module_ForSci as fsci

diagexp = """Diagnostic to compute as [diagn],[value1],[value2],...
  [laplace],[varn],[dimn1],[dimn2]: laplacian of [varn] along dimensions named [dimn1] and [dimn2]
  """
availdiag = ['[laplace],[varn],[dimn1],[dimn2]']

runsliceexp = """',' separated values of the slice as [dimn]|[slicev]
  [slicev]: chunks to process from dimension [dimn]"""

sliceexp = """',' separated values of the running slice to perform the analysis as [dimn]|[slicev]
  [value] = -1, all the values
  [value] = -9, last values
  [value] = int, a single value
  [value] = [beg, end, frq], from a beginning to an end with a given frequency"""

parser = OptionParser()
parser.add_option("-d", "--Diagnostic", dest="diag", help=diagexp, 
  metavar="LABEL")
parser.add_option("-f", "--Folder", dest="fold", help="Folder with data", 
  metavar="LABEL")
parser.add_option("-g", "--Debug", dest="debug", help="debug mode", 
  metavar="VALUE")
parser.add_option("-H", "--HeaderMiddleTile", dest="hmt", 
  help="',' separated sections of [Header],[Middle],[Taile] of the names of the files as [fold]/[Header]*[Middle]*[Tail]", 
  metavar="LABELS")
parser.add_option("-r", "--RunSlice", dest="runslice", help=sliceexp, metavar="VALUES")
parser.add_option("-s", "--Slice", dest="slice", help=sliceexp, metavar="VALUES")

(opts, args) = parser.parse_args()

mainn = 'Slicing_computations.py'

#######    #######
## MAIN
    #######
    
# Diagnostic
if opts.diag is None:
    print (errormsg)
    print ('  ' + mainn + ": no diagnostic provided !!")
    print ("     a gicvn diagnostic must be provided as -d [diag]")
    quit(-1)
else:
    diag = opts.diag.split(',')
    
# Folder
if opts.fold is None:
    print (errormsg)
    print ('  ' + mainn + ": no folder provided !!")
    print ("     a given folder must be provided as -f [folder]")
    quit(-1)

# Debug
if opts.debug is None:
    print (warnmsg)
    print ('  ' + mainn + ": no debug value is given, assuming no debugging !!")
    debug = False
else:
    debug = gen.Str_Bool(opts.debug)
    
# HMT
if opts.hmt is None:
    print (errormsg)
    print ('  ' + mainn + ": no header, middle and tail part for files to process " +\
      "provided !!")
    print ("     strings must be provided as -H [Header],[Middle],[Tail]")
    quit(-1)
else
    epxa = '[Header],[Middle],[Tail]'
    gen.check_arguments(mainn, opts.hmt, expa, ',')
    header = opts.hmt.split(',')[0]
    middle = opts.hmt.split(',')[1]
    tail = opts.hmt.split(',')[2]

# runslice
if opts.runslice is None:
    print (warnmsg)
    print ('  ' + mainn + ": no runnig slice to avoid memory issues is given, " +    \
      "assuming to process all the data at once!!")
    runslice = None
else
    runslice = opts.runslice + ''

# slice
if opts.slice is None:
    print (warnmsg)
    print ('  ' + mainn + ": no slice to data is given processing all data")
    slicev = None
else
    slicev = opts.slice + ''

# Getting the files
files = gen.files_folder_HMT(folder=opts.fold, head=header, middle=middle, tail=tail,\
  rmfolder=False)

# Getting the slice
varslice = gen.Str_DicSlice(opts.slice, ',', '|')

# Including the running dimension
runslice = stringS_dictvar(opts.runslice, ',', '|')

for dimn in runslice.keys():
    rdicv = runslice[dimn]
    if dimn in varslice:
        dicv = varslice[dimn]
        if type(dicv) == type(1):
            if dicv == -1:
                varslice[dimn] = [0,dicv,rdicv]
    else:
        varslice[dimn] = [0, -1, rdicv]

if diag[0] == 'laplace':
    


else:
    print (errormsg)
    print ('  ' + mainn + ": diagnostic '" + diag[0] + "' not ready !!")
    print ('    available ones:', availdiag)
    quit(-1)

