import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table
import astropy.units as u
import astropy.coordinates as ac
from scipy.optimize import leastsq

import plot_util as pp

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/source_class/hetdex_v2/'
#path = './hetdex_v2/'

# Gaus catalogue
lofargcat = Table.read(path+'mosaic.pybdsm.gaul.fits')

# Source catalogue
lofarcat = Table.read(path+'mosaic.pybdsm.srl.fits')


maskDC0 = lofarcat['DC_Maj'] == 0

maskS = lofarcat['S_Code'] == 'S'
maskM = lofarcat['S_Code'] == 'M'
maskC = lofarcat['S_Code'] == 'C'
Nsources = len(lofarcat)
NS = np.sum(maskS)
NM = np.sum(maskM)
NC = np.sum(maskC)

print '{n:d} sources'.format(n=Nsources)
print '{n:d} S'.format(n=NS)
print '{n:d} M'.format(n=NM)
print '{n:d} C'.format(n=NC)


# get nearest neighbour for all sources
c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,c,nthneighbor=2)



# source classes - parameters & masks

size_large = 15.           # in arcsec

# compact isolated S sources
size1 = 5.           # in arcsec
separation1 = 30.    # in arcsec
mask1 = (lofarcat['DC_Maj']*3600. < size1) & (f_nn_sep2d > separation1*u.arcsec) & (maskS)
# LR

# extended isolated S sources
size2 = 15.          # arcsec -- how large can this be
mask2 = (lofarcat['DC_Maj']*3600. >= size1) & (lofarcat['DC_Maj']*3600. < size2) & (f_nn_sep2d > separation1*u.arcsec) & (maskS) 
# LR


# isolated M sources
separation3 = 30.    # in arcsec
mask3 = (f_nn_sep2d > separation3*u.arcsec) & (~maskS) & (lofarcat['DC_Maj']*3600. <= size_large)
#

# not isolated
mask4 = (f_nn_sep2d <= separation3*u.arcsec) & (lofarcat['DC_Maj']*3600. <= size_large)
# LR/VC ?? to investigate further

# large extended
mask5 = (lofarcat['DC_Maj']*3600. >= size_large)
#VC

print '# Source classes #'
print '{n:d} S compact isolated (s<5", NN>30")'.format(n=np.sum(mask1))
print '{n:d} S extended isolated (5<s<15", NN>30")'.format(n=np.sum(mask2))
print '{n:d} M isolated (s<15", NN>30")'.format(n=np.sum(mask3))
print '{n:d} not isolated (NN<30")'.format(n=np.sum(mask4))
print '{n:d} large (s>15")'.format(n=np.sum(mask5))

# all visual sources




### diagnostic plots ###

# plot size distribution
f, ax = pp.paper_single_ax()
ax.hist(3600.*lofarcat['DC_Maj'][~maskDC0], range=(0,80), bins=100, histtype='step', label='All')
ax.hist(3600.*lofarcat['DC_Maj'][~maskDC0&maskS], range=(0,80), bins=100, histtype='step', label='S')
ax.hist(3600.*lofarcat['DC_Maj'][~maskDC0&maskM], range=(0,80), bins=100, histtype='step', label='M')
ax.hist(3600.*lofarcat['DC_Maj'][~maskDC0&maskC], range=(0,80), bins=100, histtype='step', label='C')
ax.set_xlabel('Deconvolved Major Axis [arcsec]')
ax.set_ylabel('N')
ax.legend()


# plot nearest neighbour distribution
f,ax = pp.paper_single_ax()
ax.hist(f_nn_sep2d.to('arcsec').value, bins=100, histtype='step', label='All')
ax.hist(f_nn_sep2d.to('arcsec').value[maskS], bins=100, histtype='step', label='S')
ax.set_xlabel('Nearest source [arcsec]')
ax.set_ylabel('N')
ax.legend()


# 2D histogram : size-nearest neighbour distance
# for 'S' sources
f,ax = pp.paper_single_ax()
X =  f_nn_sep2d.to('arcsec').value[~maskDC0&maskS]
Y = 3600.*lofarcat['DC_Maj'][~maskDC0&maskS]
H, xe, ye =  np.histogram2d( X, Y, bins=(100,100), normed=True)
H2 = H.T
xc = (xe[1:] +xe[:-1] )/2.
yc = (ye[1:] +ye[:-1] )/2.
c = ax.contour(xc, yc, H2, [0.5])
xind = np.sum(X>xe[:,np.newaxis],axis=0)-1
yind = np.sum(Y>ye[:,np.newaxis],axis=0)-1
Hval = H2[yind,xind]
c = ax.scatter(X, Y,c=Hval,s=10, edgecolor='none',zorder=1)
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
ax.hlines(size1,x1,x2,colors='k',linestyle='dashed')
ax.vlines(separation1,y1,y2,colors='k',linestyle='dashed')
ax.set_xlabel('NN separation [arcsec]')
ax.set_ylabel('DCmaj [arcsec]')
ax.contour(xc, yc, H2)


# and 'M' sources
f,ax = pp.paper_single_ax()
X =  f_nn_sep2d.to('arcsec').value[~maskDC0&maskM]
Y = 3600.*lofarcat['DC_Maj'][~maskDC0&maskM]
H, xe, ye =  np.histogram2d( X, Y, bins=(100,100), normed=True)
H2 = H.T
xc = (xe[1:] +xe[:-1] )/2.
yc = (ye[1:] +ye[:-1] )/2.
c = ax.contour(xc, yc, H2, [0.5])
xind = np.sum(X>xe[:,np.newaxis],axis=0)-1
yind = np.sum(Y>ye[:,np.newaxis],axis=0)-1
Hval = H2[yind,xind]
c = ax.scatter(X, Y,c=Hval,s=10, edgecolor='none',zorder=1)
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
ax.vlines(separation1,y1,y2,colors='k',linestyle='dashed')
ax.set_xlabel('NN separation [arcsec]')
ax.set_ylabel('DCmaj [arcsec]')
ax.contour(xc, yc, H2)

