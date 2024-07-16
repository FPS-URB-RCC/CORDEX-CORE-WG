## e.g. # python3 all_tamon-ver.py -v yes
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

fscratch = False
gscratch = False

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
  'Goi-nia': "Goiânia", 'Xi-an': "Xi'an", 'S-o_Paulo':  "São_Paulo"
}

# Variable to compute
varn = 'tah'

# Amount of extreme cities to plot 
Ncityxtrm = 7

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
allfilen = 'all_tamon-ver.nc'
if fscratch: sub.call('rm ' + allfilen, shell=True)
if not os.path.isfile(allfilen):

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
            rlevs = list(levs)
            rlevs.sort(reverse=True)
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
    onewnc = NetCDFFile(allfilen, 'w')

    # dimensions
    newdim = onewnc.createDimension('lon22', 2*Nx22+1)
    newdim = onewnc.createDimension('lon11', 2*Nx11+1)
    newdim = onewnc.createDimension('lat22', 2*Nx22+1)
    newdim = onewnc.createDimension('lat11', 2*Nx11+1)
    newdim = onewnc.createDimension('cit22', Ncit22)
    newdim = onewnc.createDimension('cit11', Ncit11)
    newdim = onewnc.createDimension('city', Nfiles)    
    newdim = onewnc.createDimension('mon', 12)
    newdim = onewnc.createDimension('pressure', Nlevs)
    newdim = onewnc.createDimension('pres', 2)
    newdim = onewnc.createDimension('drg', 5)
    newdim = onewnc.createDimension('StrLength', 64)   
    newdim = onewnc.createDimension('dom', Ndom)
    newdim = onewnc.createDimension('gcm', Ngcm)
    newdim = onewnc.createDimension('rcm', Nrcm)
    onewnc.sync()
    
    # Variable-dimension
    newvar = onewnc.createVariable('lon22', 'f', ('lon22'))
    newvar[:] = np.arange(-Nx22*0.22,(Nx22+1)*0.22,0.22)
    newattr=ncvar.basicvardef(newvar, 'longitude', 'Longitude relative city center', \
      'degrees_east')

    newvar = onewnc.createVariable('lon11', 'f', ('lon11'))
    newvar[:] = np.arange(-Nx11*0.11,(Nx11+1)*0.11,0.11)
    newattr=ncvar.basicvardef(newvar, 'longitude', 'Longitude relative city center', \
      'degrees_east')

    newvar = onewnc.createVariable('lat22', 'f', ('lat22'))
    newvar[:] = np.arange(-Nx22*0.22,(Nx22+1)*0.22,0.22)
    newattr=ncvar.basicvardef(newvar, 'latitude', 'Latitude relative city center',   \
      'degrees_north')

    newvar = onewnc.createVariable('lat11', 'f', ('lat11'))
    newvar[:] = np.arange(-Nx11*0.11,(Nx11+1)*0.11,0.11)
    newattr=ncvar.basicvardef(newvar, 'latitude', 'Latitude relative city center',   \
      'degrees_north')

    newvar = onewnc.createVariable('mon', 'i', ('mon'))
    newvar[:] = np.arange(12)
    newattr = ncvar.basicvardef(newvar, 'mon', 'month of the year', '-')

    newvar = onewnc.createVariable('pressure', 'f', ('pressure'))
    newvar[:] = rlevs
    newattr = ncvar.basicvardef(newvar, 'pressure', 'air pressure', 'Pa')
    newvar.setncattr('presmax', 92500.)

    newvar = onewnc.createVariable('pres', 'f', ('pres'))
    newvar[:] = [92500., 85000.]
    newattr = ncvar.basicvardef(newvar, 'pres', 'air pressure', 'Pa')
    newvar.setncattr('presmax', 92500.)

    # Variables
    newvarcityn = onewnc.createVariable('cityn', 'c', ('city', 'StrLength'))
    newattr=ncvar.basicvardef(newvarcityn, 'city_name', 'Name of the city', '-')

    newvarcityn11 = onewnc.createVariable('cityn11', 'c', ('cit11', 'StrLength'))
    newattr=ncvar.basicvardef(newvarcityn11, 'city_name', 'Name of the city', '-')
    
    newvarcityn22 = onewnc.createVariable('cityn22', 'c', ('cit22', 'StrLength'))
    newattr=ncvar.basicvardef(newvarcityn22, 'city_name', 'Name of the city', '-')

    newvardom = onewnc.createVariable('dom', 'c', ('dom', 'StrLength'))
    newattr=ncvar.basicvardef(newvardom, 'city_name', 'Name of the city', '-')

    newvargcm = onewnc.createVariable('gcm', 'c', ('gcm', 'StrLength'))
    newattr=ncvar.basicvardef(newvargcm,'gcm','Name of the global climate model','-')

    newvarrcm = onewnc.createVariable('rcm', 'c', ('rcm', 'StrLength'))
    newattr=ncvar.basicvardef(newvarrcm,'rcm','Name of the regional climate model',  \
      '-')

    newcitydrg = onewnc.createVariable('citydrg', 'i', ('drg','city'))
    newattr=ncvar.basicvardef(newcitydrg, 'citydrg','domain gcm rcm of the city','K')

    newvarcitylon = onewnc.createVariable('citylon', 'f', ('city'))
    newattr=ncvar.basicvardef(newvarcitylon,'citylon','Longitude of the city center',\
      'degrees_east')
    
    newvarcitylat = onewnc.createVariable('citylat', 'f', ('city'))
    newattr=ncvar.basicvardef(newvarcitylon,'citylat','Latitude of the city center', \
      'degrees_north')

    dimn11 = ('cit11', 'mon', 'pressure', 'lat11', 'lon11')
    newvarvals11 = onewnc.createVariable(varn + '11', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarvals11, varn, 'Potential temperature', 'K')

    dimn22 = ('cit22', 'mon', 'pressure', 'lat22', 'lon22')
    newvarvals22 = onewnc.createVariable(varn + '22', 'f', dimn22)
    newattr=ncvar.basicvardef(newvarvals22, varn, 'Potential temperature', 'K')

    dimn11 = ('cit11', 'mon', 'pres')
    newvarvalsm11 = onewnc.createVariable(varn + '11mean', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarvalsm11, varn + 'mean', 'Mean potential temperature',\
      'K')

    dimn22 = ('cit22', 'mon', 'pres')
    newvarvalsm22 = onewnc.createVariable(varn + '22mean', 'f', dimn22)
    newattr=ncvar.basicvardef(newvarvalsm22, varn + 'mean', 'Mean potential temperature',\
      'K')

    dimn11 = ('cit11', 'mon', 'pres', 'lat11', 'lon11')
    newvarvalsa11 = onewnc.createVariable(varn + '11anom', 'f', dimn11,              \
      fill_value=gen.fillValueR)
    newattr=ncvar.basicvardef(newvarvalsa11, varn + 'anom', 'Anomaly potential ' +       \
      'temperature', 'K')

    dimn22 = ('cit22', 'mon', 'pres', 'lat22', 'lon22')
    newvarvalsa22 = onewnc.createVariable(varn + '22anom', 'f', dimn22,                  \
      fill_value=gen.fillValueR)
    newattr=ncvar.basicvardef(newvarvalsa22, varn + 'anom', 'Anomaly potential ' +       \
      'temperature', 'K')

    dimn11 = ('cit11', 'lat11', 'lon11')
    newvarorog11 = onewnc.createVariable('orog11', 'f', dimn11)
    newattr=ncvar.basicvardef(newvarorog11, 'orog', 'orography', 'm')

    dimn22 = ('cit22', 'lat22', 'lon22')
    newvarorog22 = onewnc.createVariable('orog22', 'f', dimn22)
    newattr=ncvar.basicvardef(newvarorog22, 'orog', 'orography', 'm')
    
    ncvar.add_global_PyNCplot(onewnc, 'main', 'all_tamon-ver.py', '1.0', True)

    onewnc.sync()

    newvals = ncvar.writing_str_nc(newvardom, alldoms, 64)
    newvals = ncvar.writing_str_nc(newvargcm, allgcms, 64)
    newvals = ncvar.writing_str_nc(newvarrcm, allrcms, 64)

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

        print (icit, ':', icit11, icit22, '|', cityn, dom, gcm, rcm)
        ilev = 0
        for levv in rlevs:       
            lev = str(levv)
            filen = dom +'/'+ cityn +'/ta'+ lev +'_'+ gcm +'_'+ rcm + '_ymonmean.nc'
            onc = NetCDFFile(filen, 'r')

            ncvar.check_varInfile('main', filen, onc, 'ta' + str(lev)+'cyclemean')
            ovar = onc.variables['ta' + str(lev) + 'cyclemean']
            if ilev == 0:
                dx = ovar.shape[2]
                dy = ovar.shape[1]
                values = np.full((12,Nlevs,2*Nx+1,2*Nx+1), gen.fillValueF) 
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

                newvarcitylon[icit] = loncity
                newvarcitylat[icit] = latcity

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
        valanom = values[:,0:2,:,:].mean(axis=(2,3))

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
    
        icit = icit + 1
        onewnc.sync()

    Ncitygr = len(lcitygrs)
    onewnc.sync()
    #onewnc.close()
    print ("Successfull cration of file '" + allfilen + "' !!")
    
else:
    onewnc = NetCDFFile(allfilen, 'r')
    
    Ncitygr = len(onewnc.dimensions['city'])
    
    onewcitydrg = onewnc.variables['citydrg']
    onewvarcityn = onewnc.variables['cityn']
    onewvarcityn11 = onewnc.variables['cityn11']
    onewvarcityn22 = onewnc.variables['cityn22']
    onewvardom = onewnc.variables['dom']
    onewvargcm = onewnc.variables['gcm']
    onewvarrcm = onewnc.variables['rcm']
    onewvarvals11 = onewnc.variables[varn + '11']
    onewvarvals22 = onewnc.variables[varn + '22']
    onewvarvalsm11 = onewnc.variables[varn + '11mean']
    onewvarvalsm22 = onewnc.variables[varn + '22mean']
    onewvarvalsa11 = onewnc.variables[varn + '11anom']
    onewvarvalsa22 = onewnc.variables[varn + '22anom']
    onewvarorog11 = onewnc.variables['orog11']
    onewvarorog22 = onewnc.variables['orog22']
    varu = onewvarvals11.units
    
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
    opres = onewnc.variables['pres']
    pres = opres[:]
    presv = opres.presmax

    onewnc.close()

    # List of cities
    lcities = []
    for icit in range(Ncitygr):
        citynS = gen.byte_String(newvarcityn[icit,:])
        if not gen.searchInlist(lcities, citynS): lcities.append(citynS)

# Plotting

# Direct individual rounded time-series
Nrow = 1
Ncol = 2
 
for icit in range(Ncitygr):
    citynS = gen.byte_String(newvarcityn[icit,:])
    
    drg = newcitydrg[:,icit]
    if drg.mask[0]: continue
    if drg[0] != -9:
        meanv = newvarvaluesa11[drg[0],:,:,:,:]
        orog = newvarorog11[drg[0],:,:]
    else:
        meanv = newvarvaluesa22[drg[1],:,:,:,:]
        orog = newvarorog22[drg[1],:,:]

    # Masking topography points
    for it in range(meanv.shape[0]):
        for iz in range(meanv.shape[1]):
            meanv[it,iz,:,:] = np.where(orog > 50., gen.fillValueR, meanv[it,iz,:,:])
    meanv = ma.masked_equal(meanv, gen.fillValueR)
            
    domS = gen.byte_String(newvardom[drg[2],:])
    gcmS = gen.byte_String(newvargcm[drg[3],:])
    rcmS = gen.byte_String(newvarrcm[drg[4],:])
    if debug: print (icit,':', citynS, domS, gcmS, rcmS)
   
    ofign = domS + '/'  + citynS + '/urbdyn_' + varn + 'mon_' + gcmS + '_' + rcmS
    if citynS in nonASCII:
        citygS = nonASCII[citynS]
    else:
        citygS = citynS + ''
    citygS = citygS.replace('_',' ')
    
    ofignS = ofign  + '.png'
    if gscratch: sub.call('rm ' + ofignS, shell=True)
    if not os.path.isfile(ofignS):
        print ("  plotting '" + ofignS + "' ...")

        minv = gen.fillValueR
        maxv = -gen.fillValueR

        fig, axmat = plt.subplots(Nrow,Ncol)
 
        # Circular plot
        ifig = 1
        ax = plt.subplot(Nrow,Ncol,ifig)

        allmean = meanv.sum(axis=(1,2,3))
        xv = list(allmean[:]) + [allmean[0]]
        yv = list(allmean[1:12])+list(allmean[0:2])
        
        ann = allmean.min()
        anx = allmean.max()
#        for iz in range(meanv.shape[1]):
#            for iy in range(meanv.shape[2]):
#                for ix in range(meanv.shape[3]):
#                    iixv = list(meanv[:,iz,iy,ix]) + [meanv[0,iz,iy,ix]]
#                    iiyv = list(meanv[1:12,iz,iy,ix]) + list(meanv[0:2,iz,iy,ix])
#                    il = ax.plot(iixv, iiyv, '-x', color='gray')
#                    ann = meanv[:,iz,iy,ix].min()
#                    anx = meanv[:,iz,iy,ix].max()
#                    if ann < minv: minv = ann
#                    if anx > maxv: maxv = anx
        il = ax.plot(xv, yv, '-x', color='black')
        for it in range(12):
            ax.annotate(gen.shortmon[it], xy=(xv[it], yv[it]), color='k')
        if ann < minv: minv = ann
        if anx > maxv: maxv = anx

        xtrm = np.max([np.abs(minv), maxv])
        xtrm = xtrm*1.15
        ilxx = ax.plot([-xtrm,xtrm], [0,0], '-', color='black', linewidth=1.)
        ilyy = ax.plot([0,0,], [-xtrm,xtrm], '-', color='black', linewidth=1.)
        print ('circular ann:', ann, 'anx:', anx, 'xtrm', xtrm)
    
        ax.set_xlim(-xtrm, xtrm)
        ax.set_ylim(-xtrm, xtrm)

        ax.set_xlabel('month (it)')
        ax.set_ylabel('month (it+1)')
        ax.grid()

        ax.set_title('urbdynmean anual cycle', fontsize=8)
        minv = gen.fillValueR
        maxv = -gen.fillValueR

        # Linear plot
        ifig = 2
        ax = plt.subplot(Nrow,Ncol,ifig)
            
        for iz in range(meanv.shape[1]):
            if iz == 0:
                color = '#AA0000'
            else:
                color = '#0000AA'
            for iy in range(meanv.shape[2]):
                for ix in range(meanv.shape[3]):
                    ann = meanv[:,iz,iy,ix].min()
                    anx = meanv[:,iz,iy,ix].max()
                    if ann < minv: minv = ann
                    if anx > maxv: maxv = anx
                    if iy == 0 and ix == 0:
                        presS = str(int(pres[iz]/100.)) + ' hPa'
                        il = ax.plot(range(12), meanv[:,iz,iy,ix], '-x', color=color,\
                          linewidth=0.5, label=presS)
                    else:
                        il = ax.plot(range(12), meanv[:,iz,iy,ix], '-x', color=color,\
                          linewidth=0.5)
        il = ax.plot(range(12), allmean, '-x', color='black', label='sum')
    
        xtrm = np.max([np.abs(minv), maxv])

        ax2 = ax.twinx()
        ax.set_ylim(-xtrm, xtrm)
        ax2.set_ylim(-xtrm, xtrm)

        ax.set_xticks(arange(12))
        ax.set_xticklabels(gen.shortmon,fontsize=6)
        ytickspos = ax.get_yticks()
        Nytcks = len(ytickspos)
        ax.set_yticklabels(['']*Nytcks)

        ax.set_xlabel('month (it)')
        ax.set_ylabel('anomaly (' + gen.units_lunits(varu) + ')')

        ax.legend()
        ax.grid()

        ax.set_title('urbdyn(i,j,k) [' + str(meanv.shape[2]) + 'x' +                 \
          str(meanv.shape[2]) +  '] anual cycle', fontsize=8)
     
        #fig.suptitle(citygS + ' ' + gcmS + ' ' + rcmS + ' anual cycle of urbdyn ' +  \
        #  'from CORDEX-CORE',fontsize=10)
        fig.suptitle(citygS + ' ' + gcmS + ' ' + rcmS + ' anual cycle of ' + varn +  \
          ' urbdyn', fontsize=11)
     
        #ax.set_title('Anual cycle of urbdyn from CORDEX-CORE')

        drw.output_kind(kfig, ofign, True)
        if debug: sub.call('display ' + ofignS + ' &', shell=True)

        #quit()

# Direct individual rounded time-series by gcm-rcm couples
Nrow = 1
Ncol = 2

for cityn in lcities:
    icityvalues = []
    for icit in range(Ncitygr):
        citynS = gen.byte_String(newvarcityn[icit,:])
        if citynS != cityn: continue
        icityvalues.append(icit)

    ivvc = 0
    Ncityvalues = len(icityvalues)
    icit = icityvalues[ivvc]
    drg = newcitydrg[:,icit]
    if drg.mask[0]: continue
                
    domS = gen.byte_String(newvardom[drg[2],:])
    gcmS = gen.byte_String(newvargcm[drg[3],:])
    rcmS = gen.byte_String(newvarrcm[drg[4],:])

    citynS = cityn + ''
    ofign = domS + '/'  + citynS + '/urbdyn_' + varn + 'mon_all-gcm-rcm'
    ofignS = ofign  + '.png'
    if gscratch: sub.call('rm ' + ofignS, shell=True)
    if not os.path.isfile(ofignS):
        print ("  plotting '" + ofignS + "' ...")

        minv = gen.fillValueR
        maxv = -gen.fillValueR
        
        allmeanv = np.full((12,Ncityvalues), gen.fillValueR)

        fig, axmat = plt.subplots(Nrow,Ncol)
     
        # Circular plot
        ifig = 1
        ax = plt.subplot(Nrow,Ncol,ifig)
        for ivvc in range(Ncityvalues):
            icit = icityvalues[ivvc]
            drg = newcitydrg[:,icit]
            if drg.mask[0]: continue
            if drg[0] != -9:
                meanv = newvarvaluesa11[drg[0],:,:,:,:]
                orog = newvarorog11[drg[0],:,:]
            else:
                meanv = newvarvaluesa22[drg[1],:,:,:,:]
                orog = newvarorog22[drg[0],:,:]

            # Masking topography points
            for it in range(meanv.shape[0]):
                for iz in range(meanv.shape[1]):
                    meanv[it,iz,:,:] = np.where(orog > 50., gen.fillValueR,          \
                      meanv[it,iz,:,:])
            meanv = ma.masked_equal(meanv, gen.fillValueR)
                
            domS = gen.byte_String(newvardom[drg[2],:])
            gcmS = gen.byte_String(newvargcm[drg[3],:])
            rcmS = gen.byte_String(newvarrcm[drg[4],:])
            if debug: print ('    ', drg, ':', gcmS, rcmS)
       
            if citynS in nonASCII:
                citygS = nonASCII[citynS]
            else:
                citygS = citynS + ''
            citygS = citygS.replace('_',' ')
    
            allmean = meanv.sum(axis=(1,2,3))
            xv = list(allmean[:]) + [allmean[0]]
            yv = list(allmean[1:12])+list(allmean[0:2])
            allmeanv[:,ivvc] = allmean
            
            ann = allmean.min()    
            anx = allmean.max()
            colv = drw.colorsauto[drg[3]+1]
            markv = drw.pointkindsauto[drg[4]+1]
            labS = gcmS + ' ' + rcmS
            if debug: print ('        ' + labS + ':',colv, markv, '<>', ann, anx)
            il = ax.plot(xv, yv, '-', color=colv, marker=markv, label=labS)
            if ann < minv: minv = ann
            if anx > maxv: maxv = anx

            if ivvc == 2: break

        allmeanv = ma.masked_equal(allmeanv, gen.fillValueR)
        allmv = allmeanv.mean(axis=1)
        xv = list(allmv[:]) + [allmv[0]]
        yv = list(allmv[1:12])+list(allmv[0:2])

        il = ax.plot(xv, yv, '-', color='k', marker='*', linewidth=1., label='mean')
        for it in range(12):
            ax.annotate(gen.shortmon[it], xy=(xv[it], yv[it]), color='k')
    
        xtrm = np.max([np.abs(minv), maxv])
        xtrm = xtrm*1.05
        ilxx = ax.plot([-xtrm,xtrm], [0,0], '-', color='black', linewidth=1.)
        ilyy = ax.plot([0,0,], [-xtrm,xtrm], '-', color='black', linewidth=1.)
        if debug: print ('circular ann:', ann, 'anx:', anx, 'xtrm', xtrm)
        
        ax.set_xlim(-xtrm, xtrm)
        ax.set_ylim(-xtrm, xtrm)
  
        ax.set_xlabel('month (it)')
        ax.set_ylabel('month (it+1)')
        ax.grid()
        ax.legend(ncol=1, fontsize=6)
    
        ax.set_title('urbdynmean anual cycle', fontsize=8)
        minv = gen.fillValueR
        maxv = -gen.fillValueR
    
        # Linear plot
        ifig = 2
        ax = plt.subplot(Nrow,Ncol,ifig)
                
        for ivvc in range(Ncityvalues):
            icit = icityvalues[ivvc]
            drg = newcitydrg[:,icit]
            if drg.mask[0]: continue
            if drg[0] != -9:
                meanv = newvarvaluesa11[drg[0],:,:,:,:]
                orog = newvarorog11[drg[0],:,:]
            else:
                meanv = newvarvaluesa22[drg[1],:,:,:,:]
                orog = newvarorog22[drg[0],:,:]

            # Masking topography points
            for it in range(meanv.shape[0]):
                for iz in range(meanv.shape[1]):
                    meanv[it,iz,:,:] = np.where(orog > 50., gen.fillValueR,          \
                      meanv[it,iz,:,:])
            meanv = ma.masked_equal(meanv, gen.fillValueR)

            allmean = meanv.sum(axis=(1,2,3))
            ann = allmean.min()
            anx = allmean.max()

            colv = drw.colorsauto[drg[3]+1]
            markv = drw.pointkindsauto[drg[4]+1]
            labS = gcmS + ' ' + rcmS
            il = ax.plot(range(12), allmean, '-', color=colv, marker=markv,          \
              linewidth=0.5, label=labS)
            if ann < minv: minv = ann
            if anx > maxv: maxv = anx
        
        il = ax.plot(range(12), allmv, '-', color='k', marker='*', linewidth=1.,     \
          label=labS)
        xtrm = np.max([np.abs(minv), maxv])
        xtrm = xtrm*1.05
        il = ax.plot([0,11], [0,0], '-', color='k', linewidth=0.75)
    
        ax2 = ax.twinx()
        ax.set_ylim(-xtrm, xtrm)
        ax2.set_ylim(-xtrm, xtrm)
    
        ax.set_xticks(arange(12))
        ax.set_xticklabels(gen.shortmon,fontsize=6)
        ytickspos = ax.get_yticks()
        Nytcks = len(ytickspos)
        ax.set_yticklabels(['']*Nytcks)
    
        ax.set_xlabel('month (it)')
        ax.set_ylabel('anomaly (' + gen.units_lunits(varu) + ')')
    
        ax.grid()
        #ax.legend()
    
        ax.set_title('urbdyn(i,j,k) [' + str(meanv.shape[2]) + 'x' +                 \
          str(meanv.shape[2]) +  '] anual cycle', fontsize=8)
         
        #fig.suptitle(citygS + ' ' + gcmS + ' ' + rcmS + ' anual cycle of urbdyn ' +  \
        #  'from CORDEX-CORE',fontsize=10)
        fig.suptitle(citygS + ' anual cycle of ' + varn +  \
          ' urbdyn', fontsize=11)
         
        #ax.set_title('Anual cycle of urbdyn from CORDEX-CORE')
    
        drw.output_kind(kfig, ofign, True)
        if debug: sub.call('display ' + ofignS + ' &', shell=True)
    
        #quit()

# Direct all rounded time-series
Nrow = 1
Ncol = 2
 
ofign = 'urbdyn_' + varn + 'mon_all'
ofignS = ofign  + '.png'
cityxtrms = gen.order_cols()
if gscratch: sub.call('rm ' + ofignS, shell=True)
if not os.path.isfile(ofignS):
    print ("  plotting '" + ofignS + "' ...")

    minv = gen.fillValueR
    maxv = -gen.fillValueR
    xtrv = -gen.fillValueR

    fig, axmat = plt.subplots(Nrow,Ncol)
 
    ifig = 1
    # Circular
    ax = plt.subplot(Nrow,Ncol,ifig)
    for icit in range(Ncitygr):
        citynS = gen.byte_String(newvarcityn[icit,:])
    
        drg = newcitydrg[:,icit]
        if drg.mask[0]: continue
        if drg[0] != -9:
            meanv = newvarvaluesa11[drg[0],:,:,:,:]
            orog = newvarorog11[drg[0],:,:]
        else:
            meanv = newvarvaluesa22[drg[1],:,:,:,:]
            orog = newvarorog22[drg[0],:,:]

        # Masking topography points
        if np.any(orog.std() > 50.): continue
# This decompensates sums
#        for it in range(meanv.shape[0]):
#            for iz in range(meanv.shape[1]):
#                meanv[it,iz,:,:] = np.where(orog > 50., gen.fillValueR,              \
#                  meanv[it,iz,:,:])
#        meanv = ma.masked_equal(meanv, gen.fillValueR)

        gcmS = gen.byte_String(newvargcm[drg[3],:])
        rcmS = gen.byte_String(newvarrcm[drg[4],:])
        cityndrg = citynS + '_' + gcmS + '_' + rcmS
            
        allmean = meanv.sum(axis=(1,2,3))
        xv = list(allmean[:]) + [allmean[0]]
        yv = list(allmean[1:12])+list(allmean[0:2])

        ann = allmean.min()
        anx = allmean.max()
        if icit == Ncitygr-1:
            il = ax.plot(xv, yv, '-x', color='black')
            for it in range(12):
                ax.annotate(gen.shortmon[it], xy=(xv[it], yv[it]), color='red')
        else:
            il = ax.plot(xv, yv, '-x', color='gray')
    
        if ann < minv: 
            minv = ann
            if debug: print ('  new minium for ' + citynS, minv)
        if anx > maxv: 
            maxv = anx
            if debug: print ('  new maxium for ' + citynS, maxv)
        # Absolute xtreme
        xtr = np.max([np.abs(ann),anx])
        if icit == 0:
            cityxtrms.first([xtr,icit,cityndrg],'reverse')
        else:
            cityxtrms.fill([xtr,icit,cityndrg])

        #if icit == 0: break

    if debug:
        print ('    first', 2*Ncityxtrm, 'xtreme cities _______')
        dicvdrg = cityxtrms.cols[2]
        dicvval = cityxtrms.cols[0]
        dicvicit = cityxtrms.cols[1]
        for iic in range(Ncityxtrm*2):
            print ('    ', iic,dicvdrg[iic], dicvval[iic], dicvicit[iic])
    
    # re-Plotting first 2 cities from cityxtrms
    dicv = cityxtrms.cols[1]
    cityndrg = cityxtrms.cols[2]
    plotted = []
    iivv = 0
    for icit in dicv:
        citynS = gen.byte_String(newvarcityn[icit,:])
        if gen.searchInlist(plotted, cityndrg[iivv]): continue
        plotted.append(cityndrg[iivv])
    
        if citynS in nonASCII:
            citygS = nonASCII[citynS]
        else:
            citygS = citynS + ''
        citygS = citygS.replace('_',' ')
    
        drg = newcitydrg[:,icit]
        if drg.mask[0]: continue
        if drg[0] != -9:
            meanv = newvarvaluesa11[drg[0],:,:,:,:]
        else:
            meanv = newvarvaluesa22[drg[1],:,:,:,:]
            
        allmean = meanv.sum(axis=(1,2,3))
        xv = list(allmean[:]) + [allmean[0]]
        yv = list(allmean[1:12])+list(allmean[0:2])

        colv = drw.colorsauto[iivv]
        il = ax.plot(xv, yv, '-x', color=colv, label=citygS)
        for it in range(12):
            ax.annotate(gen.shortmon[it], xy=(xv[it], yv[it]), color=colv, fontsize=6)
        iivv = iivv + 1
        if iivv > Ncityxtrm: break

    xtrm = np.max([np.abs(minv), maxv])
    xtrm = xtrm*1.05
    if debug: print ('ann:', ann, 'anx:', anx, 'xtrm', xtrm)
    
    ax.set_xlim(-xtrm, xtrm)
    ax.set_ylim(-xtrm, xtrm)

    ax.set_xlabel('month (it)')
    ax.set_ylabel('month (it+1)')
    ax.grid()
    ax.legend(ncol=2, fontsize=6)

    ax.set_title('anual cycle')
    minv = gen.fillValueR
    maxv = -gen.fillValueR

    ifig = ifig + 1
    # Linear
    ax = plt.subplot(Nrow,Ncol,ifig)
    for icit in range(Ncitygr):
        citynS = gen.byte_String(newvarcityn[icit,:])
    
        drg = newcitydrg[:,icit]
        if drg.mask[0]: continue
        if drg[0] != -9:
            meanv = newvarvaluesa11[drg[0],:,:,:,:]
            orog = newvarorog11[drg[0],:,:]
        else:
            meanv = newvarvaluesa22[drg[1],:,:,:,:]
            orog = newvarorog22[drg[0],:,:]

        # Masking topography points
        if np.any(orog.std() > 50.): continue
            
        allmean = meanv.sum(axis=(1,2,3))

        ann = allmean.min()
        anx = allmean.max()
        il = ax.plot(range(12), allmean, '-x', color='gray')   
    
    # re-Plotting first 2 cities from cityxtrms
    dicv = cityxtrms.cols[1]
    cityndrg = cityxtrms.cols[2]
    plotted = []
    iivv = 0
    for icit in dicv:
        citynS = gen.byte_String(newvarcityn[icit,:])
        if gen.searchInlist(plotted, cityndrg[iivv]): continue
        plotted.append(cityndrg[iivv])
    
        if citynS in nonASCII:
            citygS = nonASCII[citynS]
        else:
            citygS = citynS + ''
        citygS = citygS.replace('_',' ')
    
        drg = newcitydrg[:,icit]
        if drg.mask[0]: continue
        if drg[0] != -9:
            meanv = newvarvaluesa11[drg[0],:,:,:,:]
        else:
            meanv = newvarvaluesa22[drg[1],:,:,:,:]

        allmean = meanv.sum(axis=(1,2,3))

        colv = drw.colorsauto[iivv]
        il = ax.plot(range(12), allmean, '-x', color=colv, label=citygS)
        iivv = iivv + 1
        if iivv > Ncityxtrm: break

    ax2 = ax.twinx()
    ax.set_ylim(-xtrm, xtrm)
    ax2.set_ylim(-xtrm, xtrm)
    
    ax.set_xticks(arange(12))
    ax.set_xticklabels(gen.shortmon,fontsize=6)
    ytickspos = ax.get_yticks()
    Nytcks = len(ytickspos)
    ax.set_yticklabels(['']*Nytcks)
    
    ax.set_xlabel('month (it)')
    ax.set_ylabel('anomaly (' + gen.units_lunits(varu) + ')')
    
    ax.grid()
    #ax.legend()
    
    ax.set_title('urbdyn(i,j,k) anual cycle', fontsize=8)

    fig.suptitle('Anual cycle of urbdyn ' + varn)

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)

