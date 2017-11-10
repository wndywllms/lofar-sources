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

### nhuge_2masx

nhuge_2masx_vc_cat_file = 'fullsample/sample_all_src_clean_large_faint_nhuge_2masx-vflag.fits'
nhuge_2masx_vc_cat = Table.read(nhuge_2masx_vc_cat_file)

if 'nhuge_2masx_flag' in lofarcat.colnames:
    lofarcat.remove_column('nhuge_2masx_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, nhuge_2masx_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','nhuge_2masx_flag')


lofarcat.add_column(tt['nhuge_2masx_flag'])

### clustered


clustered_vc_cat_file = 'fullsample/sample_all_src_clean_small_nisol_clustered-vflag.fits'
clustered_vc_cat = Table.read(clustered_vc_cat_file)

if 'clustered_flag' in lofarcat.colnames:
    lofarcat.remove_column('clustered_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, clustered_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','clustered_flag')


lofarcat.add_column(tt['clustered_flag'])



### large faint clustered


Lclustered_vc_cat_file = 'testsample_large/sample_all_src_clean_large_faint_nhuge_n2masx_nisol_clustered-vflag.fits'
Lclustered_vc_cat = Table.read(Lclustered_vc_cat_file)

if 'Lclustered_flag' in lofarcat.colnames:
    lofarcat.remove_column('Lclustered_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, Lclustered_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','Lclustered_flag')


lofarcat.add_column(tt['Lclustered_flag'])




#################################################################################



## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)