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
import module_ForDiag as fdiag
import module_ForGen as fgen
import module_ForSci as fsci

diagexp = """Diagnostic to compute as
  [laplace],[dimn1],[dimn2]: laplacian along dimensions named [dimn1] and [dimn2]
  """

runsliceexp = """',' separated values of the slice as [dimn]|[slicev]
  [value] = -1, all the values
  [value] = -9, last values
  [value] = int, a single value
  [value] = [beg, end, frq], from a beginning to an end with a given frequency"""

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
parser.add_option("-HMT", "--HeaderMiddleTile", dest="hmt", 
  help="',' separated sections of [Header],[Middle],[Taile] of the names of the files as [fold]/[Header]*[Middle]*[Tail]", 
  metavar="LABELS")
parser.add_option("-r", "--RunSlice", dest="runslice", help=sliceexp, metavar="VALUES")
parser.add_option("-s", "--Slice", dest="slice", help=sliceexp, metavar="VALUES")

(opts, args) = parser.parse_args()

mainn = 'Slicing_computations.py'

#######    #######
## MAIN
    #######

# Debug
if opts.debug is None:
    print (warnmsg)
    print ('  ' + mainn + ": no debug value is given, assuming no debugging !!")
    debug = False
else:
    debug = Str_Bool(opts.debug)

files_folder_HMT(folder='.',head='',middle='',tail='', rmfolder=True)
