import os 
from astropy.table import Table, Column, join, vstack
import numpy as np
'''
ID_flag
1 - ML
2 - 2MASX
3 - LGZ
30 - LGZ pending
4 - no id possible
5 - TBC
'''



'''
From the LOFAR catalogue:
Source_Name
RA
E_RA
E_RA_tot
DEC
E_DEC
E_DEC_tot
Peak_flux
E_Peak_flux
E_Peak_flux_tot
Total_flux
E_Total_flux
E_Total_flux_tot
Maj
E_Maj
Min
E_Min
PA
E_PA
Isl_rms
S_Code
Mosaic_ID
Isl_id

Append:
ID_flag
ID_wisename
ID_psname
ID_2masxname
ID_ra
ID_dec
ML_LR
LGZ_flags?

'''






if __name__=='__main__':

    ### Required INPUTS
    # lofar source catalogue, gaussian catalogue and ML catalogues for each


    path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'

    lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.gaus.fixed.fits'
    lofarcat_orig_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.fits'

    # PS ML - matches for sources and gaussians
    psmlcat_file = path+'lofar_pw.fixed.fits'
    psmlgcat_file = path+'lofar_gaus_pw.fixed.fits'

    # sorted output from flowchart
    lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.sorted.fits'

    # LGZ output
    #lgz_compcat_file = os.path.join(path,'LGZ_v0/HETDEX-LGZ-comps-v0.5.fits')
    lgz_cat_file = os.path.join(path,'LGZ_v0/HETDEX-LGZ-comps-v0.5-filtered.fits') # ! name !
    lgz_remove_file = os.path.join(path,'LGZ_v0/remove.txt')

    merge_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v0.1.fits')    

    lofarcat0 = Table.read(lofarcat_orig_file)
    lofarcat_sorted = Table.read(lofarcat_file_srt)
    lgz_remove = [l.rstrip() for l in open(lgz_remove_file,'r').readlines()]
    
    psmlcat = Table.read(psmlcat_file)
    
    #lgz_compcat = Table.read(lgz_compcat_file)
    lgz_cat = Table.read(lgz_cat_file)


    
    ## remove sources associated/flagged by LGZ v0
    # ideally this would just remove the components in the LGZ comp catalogue  - but using legacy catalogues mean that these don't directly map onto the new sources
    # martin has produced remove.txt to do this.


    lgz_select = np.ones(len(lofarcat_sorted), dtype=bool)
    #for si,s in enumerate(lofarcat_sorted['Source_Name']):
        #if s in lgz_remove:
            #lgz_select[si] = False
    lgz_remove = np.unique(lgz_remove)
    tlgz_remove = Table([Column(lgz_remove,'Source_Name'), Column(np.ones(len(lgz_remove)),'LGZ_remove')])
    lofarcat_sorted.sort('Source_Name')
    tc = join(lofarcat_sorted, tlgz_remove, join_type='left')
    tc.sort('Source_Name')
    lgz_select = (tc['LGZ_remove']!=1)


    print 'Removing {n:d} sources associated in LGZv0'.format(n=np.sum(~lgz_select))
    lofarcat_sorted = lofarcat_sorted[lgz_select]



    ## remove artefacts
    # all the artefacts identified and visually confirmed in the flowchart process
    print 'Throwing away {n:d} artefacts'.format(n=np.sum(lofarcat_sorted['Artefact_flag'] == 1))
    lofarcat_sorted = lofarcat_sorted[lofarcat_sorted['Artefact_flag'] == 0]
    
    print 'left with {n:d} sources'.format(n=len(lofarcat_sorted))


    ### add some needed columns
    
    lofarcat_sorted.add_column(Column(np.zeros(len(lofarcat_sorted),dtype='S60'),'ID_name'))
    #lofarcat_sorted.add_column(Column(['None']*len(lofarcat_sorted),'ID_name'))
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'ID_ra'))
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'ID_dec'))


    # handle TBC

    # handle 2MASX sources
    ## HUGE 2MASX sources need to be removed, associated and added back
    ## the rest need a flag for 2MASX
    sel2mass = (lofarcat_sorted['ID_flag']==2)
    print 'adding info for {n:d} 2MASX source matches'.format(n=np.sum(sel2mass))
    # add the 2MASXJ
    names = lofarcat_sorted['2MASX_name'][sel2mass]
    names = ['2MASXJ'+n for n in names]
    
    lofarcat_sorted['ID_name'][sel2mass] = names
    lofarcat_sorted['ID_ra'][sel2mass] = lofarcat_sorted['2MASX_ra'][sel2mass]
    lofarcat_sorted['ID_dec'][sel2mass] = lofarcat_sorted['2MASX_dec'][sel2mass]
    
    
    # handle ML sources
    lLR_thresh = 0.36
    selml = (lofarcat_sorted['ID_flag']==1) & (np.log10(1+lofarcat_sorted['LR']) > lLR_thresh)
    print 'adding info for {n:d} ML source matches'.format(n=np.sum(selml))
    
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'ML_LR'))

    
    # take the PS name over the WISE name
    # why is PS name just some number ?? - pepe?
    namesP = lofarcat_sorted['LR_name_ps'][selml]
    namesW = lofarcat_sorted['LR_name_wise'][selml]
    names = [ 'PS '+str(nP) if nP != 999999  else 'AllWISE'+nW for nP,nW in zip(namesP,namesW)]
    
    
    lofarcat_sorted['ID_name'][selml] = names
    lofarcat_sorted['ID_ra'][selml] = lofarcat_sorted['LR_ra'][selml]
    lofarcat_sorted['ID_dec'][selml] = lofarcat_sorted['LR_dec'][selml]
    lofarcat_sorted['ML_LR'][selml] = lofarcat_sorted['LR'][selml]
    
    selml = (lofarcat_sorted['ID_flag']==1) & (np.log10(1+lofarcat_sorted['LR']) <= lLR_thresh)
    print 'adding info for {n:d} ML source non-matches'.format(n=np.sum(selml))
    
    lofarcat_sorted['ID_name'][selml] = 'None'

                               

    ## add LGz v0 associated sources
    # 
    lgz_select = (lgz_cat['Compoverlap']==0)&(lgz_cat['Art_prob']<0.5)&(lgz_cat['Zoom_prob']<0.5)&(lgz_cat['Blend_prob']<0.5)&(lgz_cat['Hostbroken_prob']<0.5)
    print 'Selecting {n2:d} of {n1:d} sources in the LGZv0 catalogue to add'.format(n1=len(lgz_cat),n2=np.sum(lgz_select))
    lgz_cat = lgz_cat[lgz_select]
    lgz_cat.rename_column('optRA','ID_ra')
    lgz_cat.rename_column('optDec','ID_dec')
    lgz_cat.rename_column('OptID_Name','ID_name')
    lgz_cat.rename_column('Size','LGZ_Size')
    lgz_cat.rename_column('Assoc','LGZ_Assoc')
    lgz_cat.rename_column('Assoc_Qual','LGZ_Assoc_Qual')
    lgz_cat.rename_column('ID_Qual','LGZ_ID_Qual')
    lgz_cat.add_column(Column(3*np.ones(len(lgz_cat),dtype=int),'ID_flag'))

    mergecat = vstack([lofarcat_sorted, lgz_cat])
    print 'now we have {n:d} sources'.format(n=len(mergecat))

    print 'ID_flag counts:'
    unique, counts = np.unique(mergecat['ID_flag'], return_counts=True)
    for u,c in zip(unique, counts):
        print u,c


    ## throw away extra columns
    mergecat.keep_columns(['Source_Name', 'RA', 'E_RA', 'E_RA_tot', 'DEC', 'E_DEC', 'E_DEC_tot', 'Peak_flux', 'E_Peak_flux', 'E_Peak_flux_tot', 'Total_flux', 'E_Total_flux', 'E_Total_flux_tot', 'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'Isl_rms', 'S_Code', 'Mosaic_ID', 'Isl_id', 'ID_flag', 'ID_name', 'ID_ra', 'ID_dec', 'ML_LR', 'LGZ_Size', 'LGZ_Assoc', 'LGZ_Assoc_Qual', 'LGZ_ID_Qual'])

    
    if os.path.isfile(merge_out_file):
        os.remove(merge_out_file)
    mergecat.write(merge_out_file)