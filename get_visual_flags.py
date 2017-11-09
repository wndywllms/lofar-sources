#!/usr/bin/python

'''
get_visual_flags
add flags to catalogue based on visual inspection of subclasses of sources
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from astropy.table import Table, Column, join
import astropy.coordinates as ac
import astropy.units as u
import os

#################################################################################

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.presort.fits'



lofarcat = Table.read(lofarcat_file_srt)

nhuge_2masx_vc_cat_file = 'fullsample/sample_all_src_clean_large_faint_nhuge_2masx-vflag.fits'
nhuge_2masx_vc_cat = Table.read(nhuge_2masx_vc_cat_file)

lofarcat.sort('Source_Name')
tt=join(lofarcat, nhuge_2masx_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','nhuge_2masx_flag')


lofarcat.add_column(tt['nhuge_2masx_flag'])

#################################################################################



## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)