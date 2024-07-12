## e.g. # python3 all_tamon-ver.py
# L. Fita, UBA, CIMA, IFAECI, C. A. Buenos Aires, Argentina
### Computing a scalar index from which summarize the impact of cities in 
#     atmospheric circulation. 
# Using monthly spatial means around the city as reference, the urbdyn index computes 
#   the spatial mean of anomalies of the 2x2 grid-box around the city center up to 
#   850 hPa of 3-dimensional data.
#   It is assumed that if the city does not impact the atmospheric cirulation, around 
#   the city, the spatial anomaly of 3-dimensional variables has a zero mean 

import subprocess as sub
import numpy as np
import os
import re
from netCDF4 import Dataset as NetCDFFile
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from optparse import OptionParser
import sys
#from cStringIO import StringIO
import numpy.ma as ma
import datetime as dt
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.shapereader as shpreader

import nc_var_tools as ncvar
import diag_tools as diag
import drawing_tools as drw
import generic_tools as gen

fscratch = True

ifold = '/home/lluis.fita/process/mod/CORDEX-CORE/urban'

secs = ['map', 'SN', 'WE']

# interested levels
iilevs = ['925', '850']

kfig = 'png'

colb = 'hot_r'

tpotmn = 265. 
tpotmx = 320.

tpotn = 265. 
tpotx = 320.
pmin = 840.

tpotan = -2.
tpotax = 2.
colab = 'seismic'

lonlatbox = None

#######    #######
## MAIN
    #######
p0 = 1000.
RCp = 0.286

parser = OptionParser()
parser.add_option("-v", "--Verbose", dest="verbose", help="Need verbose", metavar="LABEL")

(opts, args) = parser.parse_args()

####### ###### ##### #### ### ## #

debug = gen.Str_Bool(opts.verbose)

Mn='_ymonmean'
files = gen.files_folder_HMT(folder=ifold, head='ta',middle=Mn, tail='nc',           \
  rmfolder=False)
Nfiles = len(files)

firstcity = True
ifile = 0
levs = []
cities = {}
lcities = []
lcitygrs = []
alldoms = []
allgcms = []
allrcms = []
# Getting list of available cities
for filen in files:

    fn = os.path.basedir(filen)
    dom = filen.split('/')[0]
    cityn = filen.split('/')[1]
    if not gen.searchInlist(alldoms, dom): alldoms.append(dom)
    if not gen.searchInlist(lcities, cityn): lcities.append(cityn)
    
    gcm = fn.split('_')[1]
    rcm = fn.split('_')[2]

    citygr = cityn + '_' + gcm + '_' + rcm
    tcitygr = (cityn, gcm, rcm)
    lcitygrs.append(tcitygr)

    if not gen.searchInlist(allgcms, gcm): allgcms.append(gcm)
    if not gen.searchInlist(allrcms, rcm): allrcms.append(rcm)
    
    if firstcity:
        iff =  dom + '/' + cityn
        
        Mn=gcm + '_' + rcm + '_ymonmean'
        cityfiles = gen.files_folder_HMT(folder=iff, head='ta',middle=Mn, tail='nc', \
          rmfolder=True)

        for cfilen in cityfiles:
            lev = re.sub("[^0-9]", "", filen.split('/')[2].split('_')[0])
            levs.append(int(lev))
        levs.sort()
        Nlevs = len(levs)
        
        prevcitygr
        firstcty = False

        cities[(cityn, gcm, rcm)] = [dom]

    else:
        if citygr != prevcitygr:
            cities[(cityn, gcm, rcm)] = [dom]
            prevcitygr = citygr

Ncities = len(lcities)
print ('amount of files:', Nfiles, 'from', Ncities, 'cities found: ', lcities)

if debug:
    print ('  city values [dom] [gcm] [rcm]_______')
    for cityn  in lcities:
        print (cities[cityn])
        
presv = levs[Nlevs-1]

# Getting anomalies
for citygr in citygrs:
    cityn = citygr[0]
    gcm = ciytgr[1]
    rcm = citygr[2]

    dicv = cities[citygr]
    
    dom = dicv[0]
    # Horizontal resolutions are different
    if dom == 'EUR-11':
        Nx = 4
    else:
        Nx = 2

    ilev = 0
    for lev in iilevs:       
        filen = dom +'/'+ cityn +'/ta'+ lev + '_' + gcm + '_' + rcm + '_ymonmean.nc'

        onc = NetCDFFile(filen, 'r')

        ncvar.check_varInfile('main', filen, onc, 'ta' + str(lev)+'cyclemean')
        ovar = onc.variables['ta' + str(lev) + 'cyclemean']
        if ilev == 0:
            dx = ovar.shape[2]
            dy = ovar.shape[1]
            values = np.full((12,2,dy,dx), gen.fillValueF) 
            olon = onc.variables['lon']
            olat = onc.variables['lat']

            lon = olon[:]
            lat = olat[:]

            rcmn = gcm + ' ' + rcm

            fxfoldn=ifxfold +'/'+ dom +'/'+ gcm + '/historical/r1i1p1/' + rcm
                
            loncity = lon[dy//2,dx//2]
            latcity = lat[dy//2,dx//2]

            # Getting orography
            orogfiles = gen.files_folder_HMT(folder=fxfoldn, head='*/fx/orog/',      \
              middle='orog', tail='nc', rmfolder=False)
            orogfn = orogfiles[0][1:len(orogfiles[0])]
            oncorog = NetCDFFile(orogfn, 'r')
            oorog = oncorog.variables['orog']
            orog = oorog[:]
            oorlon = oncorog.variables['lon']
            oorlat = oncorog.variables['lat']
            orlon = oorlon[:]
            orlat = oorlat[:]
            oncorog.close()

            [iorog, jorog], dist = gen.mindist(orlon, orlat, loncity, latcity)

            iox = iorog - Nx
            eox = iorog + Nx + 1
            ioy = jorog - Nx
            eoy = jorog + Nx + 1
            orogv = orog[ioy:eoy,iox:eox]

        ta = ovar[:] 
        values[:,ilev,:,:] = ta[:]*(p0/lev)**(RCp)
        onc.close()

        ilev = ilev + 1
    # end of levels
   
    # Computing anomalies
    valanom = values[:].mean(axis=(2,3))

    # FROM: https://en.wikipedia.org/wiki/Pressure_altitude
    fac1 = 1013.25
    fac2 = 44307.694
    fac3 = 5.25530
    topo = fac1*(1-orogv/fac2)**fac3

    avals = np.full((12,2,dy,dx), gen.fillValueR)
    # Masking outside topography
    for it in range(12):
        for iz in range(2):
            mat = values[it,iz,:,:] - valanom
            avals[it,iz,:,:] = ma.array(mat, mask=topo < presv)
        
    allamean = avals.mean(axis(1,2,3))
    levamean = avals.mean(axis(2,3))
    
    dicv = dicv + [values] + [avals] + [allamean] + [levamean]
    
    cities[citygr] = dicv

# Plotting

# Direct all rounded time-series
Nrow = 1
Ncol = 1
 
ofign = 'urbdyn_tahmon_all'
ofignS = ofign  + '.png'
if fscratch: sub.call('rm ' + ofignS, shell=True)
if not os.path.isfile(ofignS):

    minv = gen.fillValueR
    maxv = -gen.fillValueR

    fig, axmat = plt.subplots(Nrow,Ncol)
 
    ifig = 1
    ax = plt.subplot(Nrow,Ncol,ifig)
    for citygr in citygrs:
        dicv = cities[citygr]
        
        allmean = dicv[3]
        an = allmean.min()
        ax = allmean.max()
        
        if citygr == citygrs[Nfiles-1]:
            il = ax.plot(allmean[0:11], allmean[1:12], '-x', color='black')
            for it in range(11):
                ax.annotate(shortmon[it], xy=(allmean[it], allmean[it+1]),           \
                  color='red')
        else:
            il = ax.plot(allmean[0:11], allmean[1:12], '-', color='gray')
    
        if an < minv: minv = an
        if ax > maxv: maxv = ax
    
    ax.set_xlabel('month (it)')
    ax.set_xlabel('month (it+1)')
     
    ax.set_title('Anual cycle of urbdyn from CORDEX-CORE')

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)


