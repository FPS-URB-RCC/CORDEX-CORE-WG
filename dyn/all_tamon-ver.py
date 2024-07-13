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

ifxfold = '/datos/MOD/CORDEX/CMIP5'
ifold = '/home/lluis.fita/process/mod/CORDEX-CORE/urban'

# Amount of grid points around the city (double for EUR-11)
Nx22 = 2
Nx11 = 4

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

# Cities with non-ASCII characters
nonASCII = {
  'Goi-nia': "Goiânia", 'XI-an': "XI'an", 'S-o_Paulo':  "São_Paulo"
}

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
        
# netcDF creation
filen = 'all_tamon-ver.nc'
if fscratch: sub.call('rm ' + filen, shell=True)
if not os.path.isfile(filen):

    Mn='_ymonmean'
    ifoldS = './*/*'
    files = gen.files_folder_HMT(folder=ifoldS, head='ta',middle=Mn, tail='nc',      \
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
    Ncit22 = 0
    Ncit11 = 0
    # Getting list of available cities
    for filen in files:

        fn = os.path.basename(filen)
        if fn.count('None') != 0: continue
        dom = filen.split('/')[1]
        cityn = filen.split('/')[2]
    
        if dom == 'EUR-11': Ncit11 = Ncit11 + 1
        else:Ncit22 = Ncit22 + 1
    
        if not gen.searchInlist(alldoms, dom): alldoms.append(dom)
        if not gen.searchInlist(lcities, cityn): lcities.append(cityn)
    
        gcm = fn.split('_')[1]
        rcm = fn.split('_')[2]
        if debug: print ('    ', dom, cityn, gcm, rcm)

        citygr = cityn + '_' + gcm + '_' + rcm
        tcitygr = (cityn, gcm, rcm)
        lcitygrs.append(tcitygr)

        if not gen.searchInlist(allgcms, gcm): allgcms.append(gcm)
        if not gen.searchInlist(allrcms, rcm): allrcms.append(rcm)
    
        if firstcity:
            iff = './' + dom + '/' + cityn
        
            Mn=gcm + '_' + rcm + '_ymonmean'
            cityfiles = gen.files_folder_HMT(folder=iff, head='ta',middle=Mn,        \
              tail='nc', rmfolder=True)

            for cfilen in cityfiles:
                lev = re.sub("[^0-9]", "", cfilen.split('/')[3].split('_')[0])
                levs.append(int(lev))
            levs.sort()
            Nlevs = len(levs)
        
            prevcitygr = tcitygr
            firstcity = False

            cities[tcitygr] = [dom]

        else:
            if tcitygr != prevcitygr:
                cities[tcitygr] = [dom]
                prevcitygr = tcitygr

    Ncities = len(lcities)
    Ndom = len(alldoms)
    Ngcm = len(allgcms)
    Nrcm = len(allrcms)
    print ('amount of files:', Nfiles, 'from', Ncities, 'cities found: ', lcities)

    if debug:
        print ('  city values [cityn] [gcm] [rcm] [dom]_______')
        for citygr  in lcitygrs:
            print (citygr, cities[citygr])

    # netCDF
    ## 
    onewnc = NetCDFFile(filen, 'w')

    # dimensions
    newdim = onewnc.createDimension('lon22', 2*Nx22+1)
    newdim = onewnc.createDimension('lon11', 2*Nx11+1)
    newdim = onewnc.createDimension('lat22', 2*Nx22+1)
    newdim = onewnc.createDimension('lat11', 2*Nx11+1)
    newdim = onewnc.createDimension('cit22', Ncit22)
    newdim = onewnc.createDimension('cit11', Ncit11)
    newdim = onewnc.createDimension('city', Ncities)    
    newdim = onewnc.createDimension('mon', 12)
    newdim = onewnc.createDimension('pres', 2)
    newdim = onewnc.createDimension('drg', 4)
    newdim = onewnc.createDimension('StrLength', 64)   
    newdim = onewnc.createDimension('dom', Ndom)
    newdim = onewnc.createDimension('gcm', Ngcm)
    newdim = onewnc.createDimension('rcm', Nrcm)
    
    # Variable-dimension
    newvar = onewnc.createVariable('lon22', 'f', ('lon22'))
    newvar[:] = np.arange(-Nx22*0.22,(Nx22+2)*0.22)
    newattr=ncvar.basicvardef(newvar, 'longitude', 'Longitude relative city center', \
      'degrees_east')

    newvar = onewnc.createVariable('lon11', 'f', ('lon11'))
    newvar[:] = np.arange(-Nx11*0.11,(Nx11+2)*0.11)
    newattr=ncvar.basicvardef(newvar, 'longitude', 'Longitude relative city center', \
      'degrees_east')

    newvar = onewnc.createVariable('lon11', 'f', ('lon11'))
    newvar[:] = np.arange(-Nx22*0.22,(Nx22+2)*0.22)
    newattr=ncvar.basicvardef(newvar, 'longitude', 'Longitude relative city center', \
      'degrees_east')

    newvar = onewnc.createVariable('lon11', 'f', ('lon11'))
    newvar[:] = np.arange(-Nx11*0.11,(Nx11+2)*0.11)
    newattr=ncvar.basicvardef(newvar, 'longitude', 'Longitude relative city center', \
      'degrees_east')

    newvar = onewnc.createVariable('mon', 'i', ('mon'))
    newvar[:] = np.arange(12)
    newattr = ncvar.basicvardef(newvar, 'mon', 'month of the year', '-')

    newvar = onewnc.createVariable('pres', 'f', ('pres'))
    newvar[:] = [92500., 85000.]
    newattr = ncvar.basicvardef(newvar, 'pres', 'air pressure', 'Pa')
    newvar.setncattr('presmax', 92500.)

    # Variables
    newvarcityn = onewnc.createVariable('cityn', 'c', ('city', 'StrLength'))
    newattr=ncvar.basicvardef(newvarcityn11, 'city_name', 'Name of the city', '-')

    newvarcityn11 = onewnc.createVariable('cityn11', 'c', ('cit11', 'StrLength'))
    newattr=ncvar.basicvardef(newvarcityn11, 'city_name', 'Name of the city', '-')
    
    newvarcityn22 = onewnc.createVariable('cityn22', 'c', ('cit22', 'StrLength'))
    newattr=ncvar.basicvardef(newvarcityn11, 'city_name', 'Name of the city', '-')

    newvardom = onewnc.createVariable('dom', 'c', ('dom', 'StrLength'))
    newattr=ncvar.basicvardef(newvardom, 'city_name', 'Name of the city', '-')

    newvargcm = onewnc.createVariable('gcm', 'c', ('gcm', 'StrLength'))
    newattr=ncvar.basicvardef(newvargcm,'gcm','Name of the global climate model','-')

    newvarrcm = onewnc.createVariable('rcm', 'c', ('rcm', 'StrLength'))
    newattr=ncvar.basicvardef(newvarrcm,'rcm','Name of the regional climate model',  \
      '-')

    newcitydrg = onewnc.createVariable('citydrg', 'i', ('drg','city'))
    newattr=ncvar.basicvardef(newcitydrg, 'citydrg','domain gcm rcm of the city','K')

    dimn11 = ('cit11', 'mon', 'pres', 'lat11', 'lon11')
    newvarvals11 = onewnc.createVariable('tah11', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarvals11, 'tah', 'Potential temperature', 'K')

    dimn22 = ('cit22', 'mon', 'pres', 'lat22', 'lon22')
    newvarvals22 = onewnc.createVariable('tah22', 'f', dimn22)
    newattr=ncvar.basicvardef(newvarvals22, 'tah', 'Potential temperature', 'K')

    dimn11 = ('cit11', 'mon', 'pres')
    newvarvalsm11 = onewnc.createVariable('tah11mean', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarvalsm11, 'tahmean', 'Mean potential temperature',\
      'K')

    dimn22 = ('cit22', 'mon', 'pres')
    newvarvalsm22 = onewnc.createVariable('tah22mean', 'f', dimn22)
    newattr=ncvar.basicvardef(newvarvalsm22, 'tahmean', 'Mean potential temperature',\
      'K')

    dimn11 = ('cit11', 'mon', 'pres', 'lat11', 'lon11')
    newvarvalsa11 = onewnc.createVariable('tah11anom', 'f', dimn11,                  \
      fill_Value=gen.fillValueR)
    newattr=ncvar.basicvardef(newvarvalsa11, 'tahanom', 'Anomaly potential ' +       \
      'temperature', 'K')

    dimn22 = ('cit22', 'mon', 'pres', 'lat22', 'lon22')
    newvarvalsa22 = onewnc.createVariable('tah22anom', 'f', dimn22,                  \
      fill_Value=gen.fillValueR)
    newattr=ncvar.basicvardef(newvarvalsa22, 'tahanom', 'Anomaly potential ' +       \
      'temperature', 'K')

    dimn11 = ('cit11', 'lat11', 'lon11')
    newvarorog11 = onewnc.createVariable('orog11', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarorog11, 'orog', 'orography', 'm')

    dimn22 = ('cit22', 'lat22', 'lon22')
    newvarorog22 = onewnc.createVariable('orog22', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarorog22, 'orog', 'orography', 'm')
    
    ncvar.add_global_PyNCplot(onewnc, 'main', 'all_tamon-ver.py', '1.0', True)

    onewnc.sync()

    # Maximum pressure value
    presv = levs[Nlevs-1]
    
    # Getting anomalies
    icit = 0
    icit11 = 0
    icit22 = 0
    for citygr in lcitygrs:
        cityn = citygr[0]
        gcm = citygr[1]
        rcm = citygr[2]

        dicv = cities[citygr]
    
        dom = dicv[0]
        newvarcityn[icit,0:len(cityn)] = cityn[:]
        newvardom[icit,0:len(dom)] = dom[:]
        newvargcm[icit,0:len(gcm)] = gcm[:]
        newvarrcm[icit,0:len(rcm)] = rcm[:]
    
        newcitydrg[2,icit] = alldoms.index(dom)
        newcitydrg[3,icit] = allgcms.index(gcm)
        newcitydrg[4,icit] = allrcms.index(rcm)
    
        # Horizontal resolutions are different
        if dom == 'EUR-11':
            Nx = Nx11
            newvarcityn11[icit11,0:len(cityn)] = cityn[:]
            newcitydrg[0,icit] = icit11
            newcitydrg[1,icit] = -9
            icit11 = icit11 + 1
        else:
            Nx = Nx22
            newvarcityn22[icit22,0:len(cityn)] = cityn[:]
            newcitydrg[0,icit] = -9
            newcitydrg[1,icit] = icit22
            icit22 = icit22 + 1

        ilev = 0
        for lev in iilevs:       
            filen = dom +'/'+ cityn +'/ta'+ lev +'_'+ gcm +'_'+ rcm + '_ymonmean.nc'
            onc = NetCDFFile(filen, 'r')

            ncvar.check_varInfile('main', filen, onc, 'ta' + str(lev)+'cyclemean')
            ovar = onc.variables['ta' + str(lev) + 'cyclemean']
            if ilev == 0:
                dx = ovar.shape[2]
                dy = ovar.shape[1]
                values = np.full((12,2,2*Nx+1,2*Nx+1), gen.fillValueF) 
                olon = onc.variables['lon']
                olat = onc.variables['lat']

                lon = olon[:]
                lat = olat[:]

                rcmn = gcm + ' ' + rcm

                fxfoldn=ifxfold +'/'+ dom +'/'+ gcm + '/historical/r1i1p1/' + rcm
                
                loncity = lon[dy//2,dx//2]
                latcity = lat[dy//2,dx//2]
                ix = dx//2 - Nx
                ex = dx//2 + Nx + 1
                iy = dy//2 - Nx
                ey = dy//2 + Nx + 1

                # Getting orography
                orogfiles = gen.files_folder_HMT(folder=fxfoldn, head='*/fx/orog/',  \
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

            ta = ovar[:,iy:ey,ix:ex] 
            levv = np.float(lev)
            values[:,ilev,:,:] = ta[:]*(p0/levv)**(RCp)
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

        avals = np.full((12,2,2*Nx+1,2*Nx+1), gen.fillValueR)
        # Masking outside topography
        for it in range(12):
            for iz in range(2):
                mat = values[it,iz,:,:] - valanom[it,iz]
                avals[it,iz,:,:] = ma.array(mat, mask=topo < presv)
            
        allamean = avals.mean(axis=(1,2,3))
        levamean = avals.mean(axis=(2,3))
    
        dicv = dicv + [values] + [avals] + [allamean] + [levamean]
    
        cities[citygr] = dicv

        if dom == 'EUR-11':
            newvarvals11[icit11-1,:,:,:,:] = values
            newvarvalsa11[icit11-1,:,:,:,:] = avals
            newvarvalsm11[icit11-1,:,:] = valanom
            newvarorog11[icit11-1,:,:] = orogv
        else:
            newvarvals22[icit22-1,:,:,:,:] = values
            newvarvalsa22[icit22-1,:,:,:,:] = avals
            newvarvalsm22[icit22-1,:,:] = valanom
            newvarorog22[icit22-1,:,:] = orogv
    
        icity = icity + 1
        onewnc.sync()

    Ncitygr = len(lcitygrs)
    onewnc.sync()
    onewnc.close()
    print ("Successfull cration of file '" + filen + "' !!")
    
else
    onewnc = NetCDFFile(filen, 'r')
    
    Ncitygr = len(onewnc.dimensions['city'])
    
    onewcitydrg = onewnc.variables['citydrg']
    onewvarcityn = onewnc.Variables['cityn']
    onewvarcityn11 = onewnc.Variables['cityn11']
    onewvarcityn22 = onewnc.Variables['cityn22']
    onewvardom = onewnc.Variables['dom']
    onewvargcm = onewnc.Variables['gcm']
    onewvarrcm = onewnc.Variables['rcm']
    onewvarvals11 = onewnc.Variables['tah11']
    onewvarvals22 = onewnc.Variables['tah22']
    onewvarvalsm11 = onewnc.Variables['tah11mean']
    onewvarvalsm22 = onewnc.Variables['tah22mean']
    onewvarvalsa11 = onewnc.Variables['tah11anom']
    onewvarvalsa22 = onewnc.Variables['tah22anom']
    onewvarorog11 = onewnc.Variables['orog11']
    onewvarorog22 = onewnc.Variables['orog22']
    
    newcitydrg = onewcitydrg[:]
    newvarcityn = onewvarcityn[:]
    newvarcityn11 = onewvarcityn11[:]
    newvarcityn22 = onewvarcityn22[:]
    newvardom = onewvardom[:]
    newvargcm = onewvargcm[:]
    newvarrcm = onewvarrcm[:]
    newvarvalues11 = onewvarvals11[:]
    newvarvalues22 = onewvarvals22[:]
    newvarvaluesm11 = onewvarvalsm11[:]
    newvarvaluesm22 = onewvarvalsm22[:]
    newvarvaluesa11 = onewvarvalsa11[:]
    newvarvaluesa22 = onewvarvalsa22[:]
    newvarorog11 = onewvarorog11[:]
    newvarorog22 = onewvarorog22[:]
    
    # Maximum pressure value
    opres = oewnc.variables['pres']
    presv = opres.presmax

    onewnc.close()

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
    for icit in range(Ncitygr):
    
        drg = newcitydrg[:,icit]
        if drg[0] != -9:
            meanv = newvarvaluesm11[drg[0],:,:,:,:]
        else:
            meanv = newvarvaluesm22[drg[1],:,:,:,:]
            
        allmean = meanv.mean(axis=(1,2,3))
        xv = allmean[:]
        yv = allmean[1:12]+[allmean[0]
        
        ann = allmean.min()
        anx = allmean.max()
        if icit == Ncitygr:
        
            il = ax.plot(xv, yv, '-x', color='black')
            for it in range(12):
                ax.annotate(gen.shortmon[it], xy=(xv[it], yv[it]), color='red')
        else:
            il = ax.plot(allmean[0:11], allmean[1:12], '-', color='gray')
    
        if ann < minv: minv = ann
        if anx > maxv: maxv = anx
    
    xtrm = np.max([np.abs(ann), anx])
    xtrm = xtrm*1.05
    
    ax.set_xlim(-xtrm, xtrm)
    ax.set_ylim(-xtrm, xtrm)

    ax.set_xlabel('month (it)')
    ax.set_xlabel('month (it+1)')
     
    ax.set_title('Anual cycle of urbdyn from CORDEX-CORE')

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)


