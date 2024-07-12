## e.g. # python3 tamon_pressure_ver.py -c Buenos_Aires -d SAM-22 -g NCC-NorESM1-M -r ICTP-RegCM4-7 -v no
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

secs = ['SN', 'WE']

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
parser.add_option("-c", "--City", dest="city", help="city name", metavar="LABEL")
parser.add_option("-d", "--Domain", dest="dom", help="domain", metavar="LABEL")
parser.add_option("-g", "--GCM", dest="gcm", help="GCM", metavar="LABEL")
parser.add_option("-r", "--RCM", dest="rcm", help="RCM", metavar="LABEL")
parser.add_option("-v", "--Verbose", dest="verbose", help="Need verbose", metavar="LABEL")

(opts, args) = parser.parse_args()

####### ###### ##### #### ### ## #

debug = gen.Str_Bool(opts.verbose)

if opts.dom == 'EUR-11':
    Nx = 6
else:
    Nx = 3

ifold = opts.dom + '/' + opts.city

Mn=opts.gcm + '_' + opts.rcm + '_ymonmean'
files = gen.files_folder_HMT(folder=ifold, head='ta',middle=Mn, tail='nc', rmfolder=True)

levs = []
for filen in files:
    lev = re.sub("[^0-9]", "", filen.split('/')[2].split('_')[0])
    levs.append(int(lev))
levs.sort()
Nlevs = len(levs)

ilev = 0
for lev in levs:
    for filen in files:
        if filen.count(str(lev)) > 0:
            onc = NetCDFFile(filen, 'r')

            ncvar.check_varInfile('main', filen, onc, 'ta' + str(lev)+'cyclemean')
            ovar = onc.variables['ta' + str(lev) + 'cyclemean']
            if ilev == 0:
                dx = ovar.shape[2]
                dy = ovar.shape[1]
                values = np.full((12,Nlevs,dy,dx), gen.fillValueF) 
                olon = onc.variables['lon']
                olat = onc.variables['lat']

                lon = olon[:]
                lat = olat[:]

                rcmn = opts.gcm + ' ' + opts.rcm

                fxfoldn=ifxfold +'/'+ opts.dom +'/'+ opts.gcm + '/historical/r1i1p1/'+\
                  opts.rcm
                
                loncity = lon[dy//2,dx//2]
                latcity = lat[dy//2,dx//2]
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
                # Getting sftlf
                #sftlffiles = gen.files_folder_HMT(folder=fxfoldn, head='*/fx/sftlf/',\
                #  middle='sftlf', tail='nc', rmfolder=False)
                #sftlffn = sftlffiles[0][1:len(sftlffiles[0])]
                #oncsftlf = NetCDFFile(sftlffn, 'r')
                #osftlf = oncsftlf.variables['sftlf']
                #sftlf = osftlf[:]
                #oncsftlf.close()

            ta = ovar[:] 
            values[:,ilev,:,:] = ta[:]*(p0/lev)**(RCp)
            onc.close()

            ilev = ilev + 1

# Anomalies
ix = dx//2 - Nx
ex = dx//2 + Nx + 1
iy = dy//2 - Nx
ey = dy//2 + Nx + 1
valanom = values[:,:,iy:ey,ix:ey].mean(axis=(2,3))

# Plotting
proj = drw.DefineMap('PlateCarree', [0.,'None'])
datacrs = drw.DefineMap('PlateCarree', [0.,'None'])
Mervals, Parvals, Coastvals, Countryvals, Stsvals, Rivervals, varprojvals,           \
 mapv, lonlatL, kinddatasource = drw.Projection_values('auto', 'auto')

Nrow = 3
Ncol = 4

Ncol, Nrow, NOpanels, dmer, dpar = drw.NcolNrow_figure(Nrow*Ncol, kind='colfix,4')
figsizev, txth = drw.figureMap_size(lonlatbox, figwidth=8, Ncol=Ncol, Nrow=Nrow,     \
  titpercen=1., dpi=200)

# Direct 2D maps
ofign = opts.dom + '/' + opts.city + '/thaymon_' + opts.gcm + '_' + opts.rcm +       \
  '_map'
ofignS = ofign  + '.png'
if fscratch: sub.call('rm ' + ofignS, shell=True)
if not os.path.isfile(ofignS):

    iox = iorog - Nx
    eox = iorog + Nx + 1
    ioy = jorog - Nx
    eoy = jorog + Nx + 1

    iz = Nlevs-1
    presv = levs[iz]

    vals = values[:,iz,iy:ey,ix:ex] + 0.
    orogv = orog[ioy:eoy,iox:eox]
    #sftlfv = sftlf[ioy:eoy,iox:eox]

    x = lon[iy:ey,ix:ex]
    y = lat[iy:ey,ix:ex]
    
    # FROM: https://en.wikipedia.org/wiki/Pressure_altitude
    fac1 = 1013.25
    fac2 = 44307.694
    fac3 = 5.25530
    topo = fac1*(1-orogv/fac2)**fac3

    # Masking outside topography
    for it in range(12):
        vals[it,:,:] = ma.array(vals[it,:,:], mask=topo < presv)
    topo = ma.masked_greater(topo, presv)

    print  ('  map range:', vals.min(), ',', vals.max())
    
    fig,axmat = plt.subplots(Nrow,Ncol)

    ifig=1
    for imon in range(12):
        ax = plt.subplot(Nrow,Ncol,ifig,projection=proj)

        # Shaded
        im = drw.plot_pcolormesh(ax, x, y, vals[imon,:,:], datacrs, colb, 'auto',   \
          tpotn, tpotx) 

        # Contour 0.25 K
        il=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotn, tpotx, 'fixc',    \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',0.25,'%g',True,1,0,0.25,'solid'], labplot=None)

        # Contour 1. K
        il1=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotn, tpotx, 'fixc',   \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',1.0,'%g',True,6,1,0.5,'solid'], labplot=None)

        # Topography
        lucolors = drw.usgscolors
        lucols = []
        cblabs = []
        Ncolors = len(list(lucolors.keys()))
        Ncc = 0
        for iv in range(Ncolors):
            vvals = lucolors[iv+1]
            # Without repeated values
            if not gen.searchInlist(cblabs,vvals[0]):
                cblabs.append(vvals[0])
                lucols.append(vvals[2])
                Ncc = Ncc + 1
        ccs = np.arange(0.,Ncc*1.,Ncc*1/(Ncc*1.+1.))+1.5
        Ncolors = len(lucols)
        bounds = range(Ncolors)
        cmap = mpl.colors.ListedColormap(lucols)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        topo1 = np.where(topo < presv, 32., -9.)
        topo1 = ma.masked_equal(topo1, -9.)

        tl = plt.pcolormesh(x, y, topo1, transform=datacrs, cmap=cmap, vmin=1,       \
          vmax=Ncolors)

        # sftlfs: this is not land cover !
        #topo1 = topo+levs[Nlevs-1]-(levs[Nlevs-1]-pmin)*0.025
        #for ii in range(Nx*2):
        #    xx = [x1D[ii],x1D[ii+1]]
        #    y1 = [topo[ii],topo[ii+1]]
        #    y2 = [topo1[ii],topo1[ii+1]]
        #    print ('  Lluis', ii, ':', xx, y1, y2, int(sftlfv[ii]))
        #    col = drw.colorsauto[int(sftlfv[ii])]
        #    sl = ax.fill_between(xx, y1, y2, color=col) 
        
        # City center
        ip = ax.annotate('x', xy=(x[Nx,Nx],y[Nx,Nx]), color='#0000AA', fontsize=10)

        if dmer[ifig-1] == 'bottom':
            dm = True
        else:
            dm = False
        if dpar[ifig-1] == 'left':
            dp = True
        else:
            dp = False

        drw.plot_mapdrawing(ax, x, y, datacrs, LonLatbox=lonlatbox, botlab=dm,       \
          leftlab=dp, kind='NaturalEarth,10m')
        #drw.plot_mapdrawing(ax, x, y, datacrs, LonLatbox=lonlatL, botlab=dm,         \
        #  leftlab=dp, mervals=Mervals, parvals=Parvals, coastvals=Coastvals,         \
        #  countryvals=Countryvals, stsvals=Stsvals, rivervals=Rivervals, kind=kinddatasource)

        ax.set_title(gen.shortmon[imon])

        ifig = ifig + 1

    drw.add_colorbar(fig, im, '', orientation='vertical', cbtkfmt='%g',              \
      cbtkrot=0., cbtksize=8, posize=[0.91, 0.30, 0.03, 0.40], labpos=[0.98,0.5],    \
      labfigcoords='figure fraction', labrot=90, labsize=10, cbticks='auto',         \
      cblabs='auto')
    ax.annotate('tha ($K$)', xy=(0.92,0.25), xycoords='figure fraction', fontsize=10)

    labS = rcmn + ' mon tahmean ' + opts.city.replace('_',' ')
    fig.suptitle(labS, fontsize=11)

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)

# Anom 2D maps
ofign = opts.dom + '/' + opts.city + '/thaymon-anom_' + opts.gcm + '_' + opts.rcm +  \
  '_map'
ofignS = ofign + '.png'
if fscratch: sub.call('rm ' + ofignS, shell=True)
if not os.path.isfile(ofignS):

    ix = dx//2 - Nx
    ex = dx//2 + Nx + 1
    iy = dy//2 - Nx
    ey = dy//2 + Nx + 1
    iox = iorog - Nx
    eox = iorog + Nx + 1
    ioy = jorog - Nx
    eoy = jorog + Nx + 1
    
    iz = Nlevs-1
    presv = levs[iz]

    vals = values[:,iz,iy:ey,ix:ex] + 0.
    orogv = orog[ioy:eoy,iox:eox]
    #sftlfv = sftlf[ioy:eoy,iox:eox]

    for it in range(12):
        vals[it,:,:] = values[it,iz,iy:ey,ix:ex] - valanom[it,iz]

    x = lon[iy:ey,ix:ex]
    y = lat[iy:ey,ix:ex]

    # FROM: https://en.wikipedia.org/wiki/Pressure_altitude
    fac1 = 1013.25
    fac2 = 44307.694
    fac3 = 5.25530
    topo = fac1*(1-orogv/fac2)**fac3

    # Masking outside topography
    for it in range(12):
        vals[it,:,:] = ma.array(vals[it,:,:], mask=topo < presv)

    print  ('  anom map range:', vals.min(), ',', vals.max())
    
    fig,axmat = plt.subplots(Nrow,Ncol)

    ifig=1
    for imon in range(12):
        ax = plt.subplot(Nrow,Ncol,ifig,projection=proj)

        # Shaded
        im = drw.plot_pcolormesh(ax, x, y, vals[imon,:,:], datacrs, colab, 'auto',    \
          tpotan, tpotax)

        # Contour 0.25 K
        il=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotan, tpotax, 'fixc',  \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',0.25,'%g',True,1,0,0.25,'solid'], labplot=None)

        # Contour 1. K
        il1=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotan, tpotax, 'fixc', \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',1.0,'%g',True,6,1,0.5,'solid'], labplot=None)

        # Topography
        lucolors = drw.usgscolors
        lucols = []
        cblabs = []
        Ncolors = len(list(lucolors.keys()))
        Ncc = 0
        for iv in range(Ncolors):
            vvals = lucolors[iv+1]
            # Without repeated values
            if not gen.searchInlist(cblabs,vvals[0]):
                cblabs.append(vvals[0])
                lucols.append(vvals[2])
                Ncc = Ncc + 1
        ccs = np.arange(0.,Ncc*1.,Ncc*1/(Ncc*1.+1.))+1.5
        Ncolors = len(lucols)
        bounds = range(Ncolors)
        cmap = mpl.colors.ListedColormap(lucols)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        topo1 = np.where(topo < presv, 32., -9.)
        topo1 = ma.masked_equal(topo1, -9.)

        tl = plt.pcolormesh(x, y, topo1, transform=datacrs, cmap=cmap, vmin=1,       \
          vmax=Ncolors)

        # sftlfs: this is not land cover !
        #topo1 = topo+levs[Nlevs-1]-(levs[Nlevs-1]-pmin)*0.025
        #for ii in range(Nx*2):
        #    xx = [x1D[ii],x1D[ii+1]]
        #    y1 = [topo[ii],topo[ii+1]]
        #    y2 = [topo1[ii],topo1[ii+1]]
        #    print ('  Lluis', ii, ':', xx, y1, y2, int(sftlfv[ii]))
        #    col = drw.colorsauto[int(sftlfv[ii])]
        #    sl = ax.fill_between(xx, y1, y2, color=col) 
        
        # City center
        ip = ax.annotate('x', xy=(x[Nx,Nx],y[Nx,Nx]), color='#0000AA', fontsize=10)

        if dmer[ifig-1] == 'bottom':
            dm = True
        else:
            dm = False
        if dpar[ifig-1] == 'left':
            dp = True
        else:
            dp = False

        drw.plot_mapdrawing(ax, x, y, datacrs, LonLatbox=lonlatbox, botlab=dm,       \
          leftlab=dp, kind='NaturalEarth,10m')
        #drw.plot_mapdrawing(ax, x, y, datacrs, LonLatbox=lonlatL, botlab=dm,         \
        #  leftlab=dp, mervals=Mervals, parvals=Parvals, coastvals=Coastvals,         \
        #  countryvals=Countryvals, stsvals=Stsvals, rivervals=Rivervals, kind=kinddatasource)

        ax.set_title(gen.shortmon[imon])

        ifig = ifig + 1

    drw.add_colorbar(fig, im, '', orientation='vertical', cbtkfmt='%g',              \
      cbtkrot=0., cbtksize=8, posize=[0.91, 0.30, 0.03, 0.40], labpos=[0.98,0.5],    \
      labfigcoords='figure fraction', labrot=90, labsize=10, cbticks='auto',         \
      cblabs='auto')
    ax.annotate('tha ($K$)', xy=(0.92,0.25), xycoords='figure fraction', fontsize=10)

    labS = rcmn + ' mon anom. tahmean ' + opts.city.replace('_',' ')
    fig.suptitle(labS, fontsize=11)

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)

# Direct vert cros sec. figures
for secn in secs:
    ofign = opts.dom + '/' + opts.city + '/thaymon_' + opts.gcm + '_' + opts.rcm +   \
      '_' + secn
    ofignS = ofign + '.png'
    if fscratch: sub.call('rm ' + ofignS, shell=True)
    if os.path.isfile(ofignS): continue

    if secn == 'SN':
        ix = dy//2 - Nx
        ex = dy//2 + Nx + 1 
        iox = jorog - Nx
        eox = jorog + Nx + 1 
        vals = values[:,:,ix:ex,dx//2] + 0.
        orogv = orog[iox:eox,iorog]
        #sftlfv = sftlf[iox:eox,iorog]

        x1D = lat[ix:ex,dx//2]
        y1D = levs
        x, y = np.meshgrid(x1D,y1D)

        xaxisS = 'latitude'

    else:
        ix = dx//2 - Nx
        ex = dx//2 + Nx + 1
        iox = iorog - Nx
        eox = iorog + Nx + 1 
        vals = values[:,:,dy//2,ix:ex] + 0.
        orogv = orog[jorog,iox:eox]
        #sftlfv = sftlf[jorog,iox:eox]

        x1D = lon[dy//2,ix:ex]
        y1D = levs
        x, y = np.meshgrid(x1D,y1D)

        xaxisS = 'longitude'

    # FROM: https://en.wikipedia.org/wiki/Pressure_altitude
    fac1 = 1013.25
    fac2 = 44307.694
    fac3 = 5.25530
    topo = fac1*(1-orogv/fac2)**fac3

#    # Masking outside topography
#    for it in range(12):
#        for iz in range(vals.shape[1]):
#            vals[it,iz,:] = ma.array(vals[it,iz,:], mask=topo < presv)

    print  ('  sec', secn, 'range:', vals.min(), ',', vals.max())
    
    fig,axmat = plt.subplots(Nrow,Ncol)

    ifig=1
    for imon in range(12):
        ax = plt.subplot(Nrow,Ncol,ifig)

        # Shaded
        im = drw.plot_pcolormeshNOmap(ax, x, y, vals[imon,:,:], colb, 'auto',        \
          tpotn, tpotx) 

        # Contour 0.25 K
        il=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotn, tpotx, 'fixc',    \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',0.25,'%g',True,1,0,0.25,'solid'], labplot=None)

        # Contour 1. K
        il1=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotn, tpotx, 'fixc',   \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',1.0,'%g',True,6,1,0.5,'solid'], labplot=None)
        # Topography
        basetopo = [levs[Nlevs-1]]*(Nx*2+1)
        tl = ax.fill_between(x1D, basetopo, topo, color='#B2B2B2', step='post') 

        # sftlfs: this is not land cover !
        #topo1 = topo+levs[Nlevs-1]-(levs[Nlevs-1]-pmin)*0.025
        #for ii in range(Nx*2):
        #    xx = [x1D[ii],x1D[ii+1]]
        #    y1 = [topo[ii],topo[ii+1]]
        #    y2 = [topo1[ii],topo1[ii+1]]
        #    print ('  Lluis', ii, ':', xx, y1, y2, int(sftlfv[ii]))
        #    col = drw.colorsauto[int(sftlfv[ii])]
        #    sl = ax.fill_between(xx, y1, y2, color=col) 
        
        # City center
        plinen=levs[Nlevs-1]-(levs[Nlevs-1]-pmin)*0.10
        ip = ax.plot([x1D[Nx],x1D[Nx]],[levs[Nlevs-1],plinen], '-', color='#0000AA', linewidth=2)

#        ax.set_ylim(levs[Nlevs-1], levs[0])
        ax.set_ylim(levs[Nlevs-1], pmin)

        if dmer[ifig-1] != 'bottom':
            Nticks = len(ax.get_xticks())
            ax.set_xticklabels(['']*Nticks)
        else:
            ax.set_xlabel(xaxisS)
        if dpar[ifig-1] != 'left':
            Nticks = len(ax.get_yticks())
            ax.set_yticklabels(['']*Nticks)
        else:
            ax.set_ylabel('pressure ($hPa$)')

        ax.set_title(gen.shortmon[imon])

        ifig = ifig + 1

    drw.add_colorbar(fig, im, '', orientation='vertical', cbtkfmt='%g',              \
      cbtkrot=0., cbtksize=8, posize=[0.91, 0.30, 0.03, 0.40], labpos=[0.98,0.5],    \
      labfigcoords='figure fraction', labrot=90, labsize=10, cbticks='auto',         \
      cblabs='auto')
    ax.annotate('tha ($K$)', xy=(0.92,0.25), xycoords='figure fraction', fontsize=10)

    labS = secn + ' vert. for ' + rcmn + ' mon tahmean ' + opts.city.replace('_',' ')
    fig.suptitle(labS, fontsize=11)

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)

for secn in secs:
    ofign = opts.dom + '/' + opts.city + '/thaymon-anom_' + opts.gcm + '_' +         \
      opts.rcm + '_' + secn
    ofignS = ofign + '.png'
    if fscratch: sub.call('rm ' + ofignS, shell=True)
    if os.path.isfile(ofignS): continue

    if secn == 'SN':
        ix = dy//2 - Nx
        ex = dy//2 + Nx + 1 
        iox = jorog - Nx
        eox = jorog + Nx + 1 
        vals = values[:,:,ix:ex,dx//2] + 0.
        for it in range(12):
            for iz in range(valanom.shape[1]):
                vals[it,iz,:] = values[it,iz,ix:ex,dx//2] - valanom[it,iz]
        orogv = orog[iox:eox,iorog]
        #sftlfv = sftlf[iox:eox,iorog]

        x1D = lat[ix:ex,dx//2]
        y1D = levs
        x, y = np.meshgrid(x1D,y1D)

        xaxisS = 'latitude'

    else:
        ix = dx//2 - Nx
        ex = dx//2 + Nx + 1
        iox = iorog - Nx
        eox = iorog + Nx + 1 
        vals = values[:,:,dy//2,ix:ex] + 0.
        for it in range(12):
            for iz in range(valanom.shape[1]):
                vals[it,iz,:] = values[it,iz,dy//2,ix:ex] - valanom[it,iz]
        orogv = orog[jorog,iox:eox]
        #sftlfv = sftlf[jorog,iox:eox]

        x1D = lon[dy//2,ix:ex]
        y1D = levs
        x, y = np.meshgrid(x1D,y1D)

        xaxisS = 'longitude'

    # FROM: https://en.wikipedia.org/wiki/Pressure_altitude
    fac1 = 1013.25
    fac2 = 44307.694
    fac3 = 5.25530
    topo = fac1*(1-orogv/fac2)**fac3

#    # Masking outside topography
#    for it in range(12):
#        for iz in range(vals.shape[0]):
#            vals[it,iz,:] = ma.array(vals[it,iz,:], mask=topo < presv)

    print  ('  sec anom', secn, 'range:', vals.min(), ',', vals.max())
    
    fig,axmat = plt.subplots(Nrow,Ncol)

    ifig=1
    for imon in range(12):
        ax = plt.subplot(Nrow,Ncol,ifig)

        # Shaded
        im = drw.plot_pcolormeshNOmap(ax, x, y, vals[imon,:,:], colab, 'auto',       \
          tpotan, tpotax)

        # Contour 0.25 K
        il=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotan, tpotax, 'fixc',  \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',0.25,'%g',True,1,0,0.25,'solid'], labplot=None)

        # Contour 1. K
        il1=drw.plot_contour(ax, vals[imon,:,:], x, y, None, tpotan, tpotax, 'fixc', \
          None, varlabpos=[0.01,0.05], varlabcoords='figure fraction',               \
          parameters=['#000000',1.0,'%g',True,6,1,0.5,'solid'], labplot=None)
        # Topography
        basetopo = [levs[Nlevs-1]]*(Nx*2+1)
        tl = ax.fill_between(x1D, basetopo, topo, color='#B2B2B2', step='post') 

        # sftlfs: this is not land cover !
        #topo1 = topo+levs[Nlevs-1]-(levs[Nlevs-1]-pmin)*0.025
        #for ii in range(Nx*2):
        #    xx = [x1D[ii],x1D[ii+1]]
        #    y1 = [topo[ii],topo[ii+1]]
        #    y2 = [topo1[ii],topo1[ii+1]]
        #    print ('  Lluis', ii, ':', xx, y1, y2, int(sftlfv[ii]))
        #    col = drw.colorsauto[int(sftlfv[ii])]
        #    sl = ax.fill_between(xx, y1, y2, color=col) 
        
        # City center
        plinen=levs[Nlevs-1]-(levs[Nlevs-1]-pmin)*0.10
        ip = ax.plot([x1D[Nx],x1D[Nx]],[levs[Nlevs-1],plinen], '-', color='#0000AA', \
          linewidth=2)

#        ax.set_ylim(levs[Nlevs-1], levs[0])
        ax.set_ylim(levs[Nlevs-1], pmin)

        if dmer[ifig-1] != 'bottom':
            Nticks = len(ax.get_xticks())
            ax.set_xticklabels(['']*Nticks)
        else:
            ax.set_xlabel(xaxisS)
        if dpar[ifig-1] != 'left':
            Nticks = len(ax.get_yticks())
            ax.set_yticklabels(['']*Nticks)
        else:
            ax.set_ylabel('pressure ($hPa$)')

        ax.set_title(gen.shortmon[imon])

        ifig = ifig + 1

    drw.add_colorbar(fig, im, '', orientation='vertical', cbtkfmt='%g',              \
      cbtkrot=0., cbtksize=8, posize=[0.91, 0.30, 0.03, 0.40], labpos=[0.98,0.5],    \
      labfigcoords='figure fraction', labrot=90, labsize=10, cbticks='auto',         \
      cblabs='auto')
    ax.annotate('tha ($K$)', xy=(0.92,0.25), xycoords='figure fraction',             \
      fontsize=10)

    labS = secn + ' vert. for ' + rcmn + ' mon anom. tahmean ' +                     \
      opts.city.replace('_',' ')
    fig.suptitle(labS, fontsize=11)

    drw.output_kind(kfig, ofign, True)
    if debug: sub.call('display ' + ofignS + ' &', shell=True)

