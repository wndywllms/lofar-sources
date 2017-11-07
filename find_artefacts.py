#!/usr/bin/python

'''
get_intersecting_sources
find which sources in a catalogue intersect, enclose or are enclosed by another source...
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from matplotlib import patches
from astropy.table import Table, Column
import astropy.coordinates as ac
import astropy.units as u
import os


path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.presort.fits'


artefactlistfile = 'gg_artefact_case1_3-fixed.fits'



lofarcat = Table.read(lofarcat_file_srt)

artefactlist = Table.read(artefactlistfile)

# for now, no artefacts
artefact = np.zeros(len(lofarcat),dtype=bool)
if 'artefact' not in lofarcat.colnames:
    lofarcat.add_column(Column(artefact,'artefact'))
for n in artefactlist['Source_Name']:
    ni = np.where(lofarcat['Source_Name']==n)[0][0]
    lofarcat['artefact'][ni] = True    

#artefact = np.zeros(len(lofarcat),dtype=bool)
##bright compact sources
#selind1 = np.where((lofarcat['Maj'] < 8.) & (lofarcat['Total_flux'] > 100.))[0]
#c = ac.SkyCoord(lofarcat['RA'][selind1], lofarcat['DEC'][selind1], unit="deg")
## faint large roundish sources
#selind = np.where((lofarcat['Maj'] > 30.) & (lofarcat['Min']/lofarcat['Maj'] > 0.5) & (lofarcat['Total_flux'] < 10.))[0]
#for ii in selind:
    #l = lofarcat[ii]
    #cl = ac.SkyCoord([l['RA']], [l['DEC']], unit="deg")
    #idx1, idxself, sep, _ =  cl.search_around_sky(c, l['Maj']*u.arcsec)
    #if len(idx1) > 0:
        #artefact[ii] = 1

#artefact = np.zeros(len(lofarcat),dtype=bool)
#selind1 = np.where((lofarcat['Maj'] < 8.) & (lofarcat['Total_flux'] > 100.))[0]
#c = ac.SkyCoord(lofarcat['RA'][selind1], lofarcat['DEC'][selind1], unit="deg")



#artefact = np.zeros(len(lofarcat),dtype=bool)
#selind1 = np.where((lofarcat['Maj'] < 8.) & (lofarcat['Total_flux'] > 100.) & (lofarcat['NN5_sep'] < 60.))[0]
#artefact[selind] =1
    

## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)