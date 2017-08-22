# coding: utf-8


#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, Column
import astropy.units as u
import astropy.coordinates as ac
import utils.plot_util as pp
import os


class mask:
    
    def __init__(self, name, mask, label):
        self.name = name
        self.mask = mask
        self.label = label
        
        self.N = self.total()
        self.n = self.msum()
        self.f = self.fraction()
        self.p = self.percent()
        
        return
    
    def percent(self):
        return 100.*np.sum(self.mask)/self.N
    
    def fraction(self):
        return 1.*np.sum(self.mask)/self.N

    def msum(self):
        return np.sum(self.mask)
    
    def total(self):
        return len(self.mask)
    
    def print_frac(self):
        print '{n:6d} ({f:4.1f}%) {label:s}'.format(n=self.n, f=self.p, label=self.label)
    


path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/source_class/t1_dr1/'
lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.1.gaus.fits'
lofarcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.1.srl.fits'
psmlcat_file = path+'lofar_matched_all.fix.fits'
psmlgcat_file = path+'lofar_matched_gaus.fits'

# Gaus catalogue
lofargcat = Table.read(lofargcat_file)
# only relevant gaussians are in M or C sources
lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

# Source catalogue
lofarcat = Table.read(lofarcat_file)

# PS ML - matches for sources and gaussians
psmlcat = Table.read(psmlcat_file)
psmlgcat = Table.read(psmlgcat_file)


#sys.exit()


## match the gaussians to the sources

# quicker to generate new unique names than match on 2 columns
# get new unique source_id by combining mosaic and src id
# replace string mosaic ID with unique int (perhaps there is a more logical mapping of mosaic name to int value)
mid = lofargcat['Mosaic_ID']
mid_unique = np.unique(mid)
mid_int = np.array([np.where(mid_unique==m)[0][0] for m in mid])
# combine with Source_id for unique ID
g_src_id_new =   10000*mid_int + lofargcat['Source_id']
lofargcat.add_column(Column(g_src_id_new, 'SID'))

mid = lofarcat['Mosaic_ID']
mid_unique = np.unique(mid)
mid_int = np.array([np.where(mid_unique==m)[0][0] for m in mid])
# combine with Source_id for unique ID
src_id_new =   10000*mid_int + lofarcat['Source_id']
lofarcat.add_column(Column(src_id_new, 'SID'))



   

## get the panstarrs ML information

# join the ps ml cat  - they have identical RA/DEC (source_names were wrong)
c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
cpsml = ac.SkyCoord(psmlcat['RA'], psmlcat['DEC'], unit="deg")
f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,cpsml,nthneighbor=1)

psmlcat = psmlcat[f_nn_idx][f_nn_sep2d==0]
lofarcat = lofarcat[f_nn_sep2d==0]

lofarcat.add_column(Column(psmlcat['lr_pc_7th'], 'LR'))

# join the ps ml gaus cat  - they have identical RA/DEC (source_names were wrong)
cg = ac.SkyCoord(lofargcat['RA'], lofargcat['DEC'], unit="deg")
cpsmlg = ac.SkyCoord(psmlgcat['RA'], psmlgcat['DEC'], unit="deg")
f_nn_idx_g,f_nn_sep2d_g,f_nn_dist3d_g = ac.match_coordinates_sky(cg,cpsmlg,nthneighbor=1)

psmlgcat = psmlgcat[f_nn_idx_g][f_nn_sep2d_g==0]
lofargcat = lofargcat[f_nn_sep2d_g==0]

lofargcat.add_column(Column(psmlgcat['lr_2'], 'LR'))


#sys.exit()

add_G = False
lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=int), 'Ng'))
if add_G:
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=list), 'G_ind'))

m_S = lofarcat['S_Code'] =='S'
minds = np.where(~m_S)[0]
for i,sid in zip(minds, lofarcat['SID'][~m_S]):
    ig = np.where(lofargcat['SID']==sid)[0]
    lofarcat['Ng'][i]= len(ig)
    
    if add_G:
        lofarcat['G_ind'][i]= ig



#############################################################################

def print_frac(m,label):
    print '{n:6d} ({f:4.1f}%) {label:s}'.format(n=np.sum(m), f=100.*np.sum(m)/len(m), label=label)


    
# # pybdsm source types

# In[4]:



# ## nearest neighbour separation

# get nearest neighbour for all sources
c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,c,nthneighbor=2)

#f_nn3_idx,f_nn3_sep2d,f_nn3_dist3d = ac.match_coordinates_sky(c,c,nthneighbor=3)
#f_nn4_idx,f_nn4_sep2d,f_nn4_dist3d = ac.match_coordinates_sky(c,c,nthneighbor=4)
f_nn5_idx,f_nn5_sep2d,f_nn5_dist3d = ac.match_coordinates_sky(c,c,nthneighbor=5)
#f_nn6_idx,f_nn6_sep2d,f_nn6_dist3d = ac.match_coordinates_sky(c,c,nthneighbor=6)


lofarcat.add_column(Column(lofarcat['LR'][f_nn_idx], 'NN_LR'))
lofarcat.add_column(Column(f_nn_sep2d.to(u.arcsec).value, 'NN_sep'))
lofarcat.add_column(Column(f_nn5_sep2d.to(u.arcsec).value, 'NN5_sep'))
lofarcat.add_column(Column(lofarcat['Total_flux'][f_nn_idx], 'NN_Total_flux'))
lofarcat.add_column(Column(lofarcat['Maj'][f_nn_idx], 'NN_Maj'))


def print_classes(lofarcat):

    Ncat = len(lofarcat)

    #maskDC0 = lofarcat['DC_Maj'] == 0
    maskDC0 = lofarcat['Maj'] == 0

    m_S = lofarcat['S_Code'] == 'S'
    m_M = lofarcat['S_Code'] == 'M'
    m_C = lofarcat['S_Code'] == 'C'

    m_N1 = lofarcat['Ng'] == 1
    m_N2 = lofarcat['Ng'] == 2
    m_N3 = lofarcat['Ng'] == 3
    m_N4 = lofarcat['Ng'] == 4
    m_N5 = lofarcat['Ng'] == 5
    m_N6p = lofarcat['Ng'] > 5


    print '{n:d} sources'.format(n=Ncat)
    print_frac(m_S, 'S')
    print_frac(m_M, 'M')
    print_frac(m_C, 'C')

    print
    print '{n:d} sources'.format(n=Ncat)
    print_frac(m_N1, '1')
    print_frac(m_N2, '2')
    print_frac(m_N3, '3')
    print_frac(m_N4, '4')
    print_frac(m_N5, '5')
    print_frac(m_N6p, '6+')

    # >15 " and 10mJY -2%

    size_large = 15.           # in arcsec
    separation1 = 60.    # in arcsec
    lLR_thresh = 1.      # LR threshold
    fluxcut = 10        # in mJy

    Ncat = len(lofarcat)
    
    m_small = lofarcat['Maj'] < size_large
    m_isol = lofarcat['NN_sep'] > separation1
    m_bright = lofarcat['Total_flux'] > fluxcut

    # compact isolated S sources
    m_small_isol_S = m_small & m_isol & (m_S)
    # LR


    # compact isolated M sources
    m_small_isol_nS = m_small & m_isol & (~m_S) 
    # LR

    # compact not isolated
    m_small_nisol = m_small & ~m_isol 
    # LR/VC ?? to investigate further

    # large extended
    m_large = ~m_small

    # large bright
    m_large_bright = ~m_small & m_bright


    # VC

    # these are mutually exclusive groups containing all sources

    m_lrgood  = np.log10(1+lofarcat['LR']) >= lLR_thresh
    m_NNlrgood  = np.log10(1+lofarcat['NN_LR']) >= lLR_thresh


    # small not isolated, but has good match
    m_small_nisol_match = m_small_nisol & m_lrgood

    # small not isolated, no good match but neighbour has good match
    m_small_nisol_NNlr = m_small_nisol & (~m_lrgood) & m_NNlrgood

    # small not isolated, no good match but neighbour has good match
    m_small_nisol_nnmatch = m_small_nisol & (~m_lrgood) & (~m_NNlrgood)

    # small not isolated, but NN has good match
    m_small_nisol_NNlr = m_small_nisol & m_NNlrgood

    print 
    print '# Source classes #'
    print Ncat

    l_lrgood = 'good LR match (log LR > {s:.0f})'.format(s=lLR_thresh)
    print_frac(m_lrgood, l_lrgood)

    l_small_isol_S = 'S compact isolated (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_isol_S, l_small_isol_S)

    l_small_isol_nS ='M compact isolated (s<{s2:.0f}", NN>{nn:.0f}")'.format(s2=size_large, nn=separation1)
    print_frac(m_small_isol_nS, l_small_isol_nS)


    l_small_nisol = 'compact not isolated (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol, l_small_nisol)


    l_small_nisol_match = 'compact not isolated good LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_match, l_small_nisol_match)

    l_small_nisol_NNlr = 'compact not isolated neighbour good LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNlr, l_small_nisol_NNlr)
    
    l_small_nisol_NNlr = 'compact not isolated bad LR, but neighbour good LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNlr, l_small_nisol_NNlr)

    l_small_nisol_nnmatch = 'compact not isolated bad LR, and neighbour bad LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_nnmatch, l_small_nisol_nnmatch)

    l_large = 'large (s>{s:.0f}")'.format(s=size_large)
    print_frac(m_large, l_large)

    l_bright = 'bright (S>{s:.0f} mJy)'.format(s=fluxcut)
    print_frac(m_bright, l_bright)


    l_large_bright = 'large (s>{s:.0f}") & bright (S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large)
    print_frac(m_large_bright, l_large_bright)

    m_large_N2 = m_large & m_N2
    l_large_N2 = 'large N2 (s>{s:.0f}")'.format(s=size_large)
    print_frac(m_large_N2, l_large_N2)

    # all visual sources

    print
    label = 'compact not isolated S (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol&m_S, label)

    label = 'compact not isolated !S (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol&(~m_S), label)


print_classes(lofarcat[np.log10(1+lofarcat['LR']) >= 1.])
print_classes(lofarcat[np.log10(1+lofarcat['LR']) < 1.])


if 1:

    Ncat = len(lofarcat)

    #maskDC0 = lofarcat['DC_Maj'] == 0
    maskDC0 = lofarcat['Maj'] == 0

    m_S = lofarcat['S_Code'] == 'S'
    m_M = lofarcat['S_Code'] == 'M'
    m_C = lofarcat['S_Code'] == 'C'

    m_N1 = lofarcat['Ng'] == 1
    m_N2 = lofarcat['Ng'] == 2
    m_N3 = lofarcat['Ng'] == 3
    m_N4 = lofarcat['Ng'] == 4
    m_N5 = lofarcat['Ng'] == 5
    m_N6p = lofarcat['Ng'] > 5


    print '{n:d} sources'.format(n=Ncat)
    print_frac(m_S, 'S')
    print_frac(m_M, 'M')
    print_frac(m_C, 'C')

    print
    print '{n:d} sources'.format(n=Ncat)
    print_frac(m_N1, '1')
    print_frac(m_N2, '2')
    print_frac(m_N3, '3')
    print_frac(m_N4, '4')
    print_frac(m_N5, '5')
    print_frac(m_N6p, '6+')

    # # source classes
    # 
    # clases from draft flowchart
    # 
    # source classes - parameters & masks

    # >15 " and 10mJY -2%

    size_large = 15.           # in arcsec
    separation1 = 60.    # in arcsec
    size_huge = 25.           # in arcsec
    #separation1 = 30.    # in arcsec
    lLR_thresh = 1.      # LR threshold
    fluxcut = 10        # in mJy

    Ncat = len(lofarcat)
    
    m_all = lofarcat['RA'] > -1
    
    masterlist = []
    
    M_all = mask('m_all',
                 lofarcat['RA'] > -1,
                 'All')
    masterlist.append(M_all)

    M_huge = mask('m_huge',
                 lofarcat['Maj'] > size_huge,
                 'Huge')
    masterlist.append(M_huge)
    
    for t in masterlist:
        t.print_frac()

    m_huge = lofarcat['Maj'] > size_huge
    m_small = lofarcat['Maj'] < size_large
    m_isol = lofarcat['NN_sep'] > separation1
    m_cluster = lofarcat['NN5_sep'] < separation1
    m_bright = lofarcat['Total_flux'] > fluxcut
    m_lrgood  = np.log10(1+lofarcat['LR']) >= lLR_thresh
    m_NNlrgood  = np.log10(1+lofarcat['NN_LR']) >= lLR_thresh


    m_n_small = lofarcat['NN_Maj'] < size_large  # neighbour small
    m_n_bright = lofarcat['NN_Total_flux'] > fluxcut  # neighbour small

    # compact isolated S sources
    m_small_isol_S = m_small & m_isol & (m_S)

    # compact isolated M sources
    m_small_isol_nS = m_small & m_isol & (~m_S) 

    # compact isolated
    m_small_isol = m_small & m_isol 
    
    # compact not isolated
    m_small_nisol = m_small & ~m_isol 

    # compact not isolated clustered
    m_small_nisol_cluster = m_small & ~m_isol & m_cluster 
    
    # large extended
    m_large = ~m_small

    # large bright
    m_large_bright = ~m_small & m_bright
        
    # large bright not huge and has LR match
    m_large_bright_nhuge_lr = ~m_small & m_bright & ~m_huge & m_lrgood
    
    # large faint
    m_large_faint = ~m_small & ~m_bright

    # small not isolated, but has good match
    m_small_isol_lr = m_small_isol & m_lrgood

    # small not isolated, but has no good match
    m_small_isol_nlr = m_small_isol & ~m_lrgood

    # small not isolated, neighbour also small
    m_small_nisol_NNsmall = m_small_nisol & m_n_small
    
    # small not isolated, neighbour also small and good lr
    m_small_nisol_NNsmall_lr = m_small_nisol & m_n_small & m_lrgood
    
    # small not isolated, neighbour also small, good lr, neighbour good lr
    m_small_nisol_NNsmall_lr_NNlr = m_small_nisol & m_n_small & m_lrgood & m_NNlrgood
    
    # small not isolated, neighbour also small, good lr, neighbour bad lr
    m_small_nisol_NNsmall_lr_NNnlr = m_small_nisol & m_n_small & m_lrgood & ~m_NNlrgood
    
    # small not isolated, neighbour also small and bad lr
    m_small_nisol_NNsmall_nlr = m_small_nisol & m_n_small & ~m_lrgood
    
    # small not isolated, neighbour also small, bad lr, neighbour good lr
    m_small_nisol_NNsmall_nlr_NNlr = m_small_nisol & m_n_small & ~m_lrgood & m_NNlrgood
    
    # small not isolated, neighbour also small, bad lr, neighbour bad lr
    m_small_nisol_NNsmall_nlr_NNnlr = m_small_nisol & m_n_small & ~m_lrgood & ~m_NNlrgood
    
    # small not isolated, neighbour large
    m_small_nisol_NNlarge = m_small_nisol & ~m_n_small

    # small not isolated, neighbour large and bright
    m_small_nisol_NNlarge_NNbright = m_small_nisol & ~m_n_small & m_n_bright

    # small not isolated, neighbour large and faint
    m_small_nisol_NNlarge_NNfaint = m_small_nisol & ~m_n_small & ~m_n_bright

    # small not isolated, but has good match
    m_small_nisol_match = m_small_nisol & m_lrgood

    # small not isolated, but has no good match
    m_small_nisol_nomatch = m_small_nisol & ~m_lrgood
    
    # small not isolated, no good match but neighbour has good match
    m_small_nisol_NNlr = m_small_nisol & (~m_lrgood) & m_NNlrgood

    # small not isolated, no good match but neighbour has good match
    m_small_nisol_nnmatch = m_small_nisol & (~m_lrgood) & (~m_NNlrgood)

    # small not isolated, but NN has good match
    m_small_nisol_NNlr = m_small_nisol & m_NNlrgood


    l_lrgood = 'good LR match (log LR > {s:.0f})'.format(s=lLR_thresh)
    print_frac(m_lrgood, l_lrgood)

    l_bright = 'bright (S>{s:.0f} mJy)'.format(s=fluxcut)
    print_frac(m_bright, l_bright)


    
    m_large_N2 = m_large & m_N2
    l_large_N2 = 'large N2 (s>{s:.0f}")'.format(s=size_large)
    print_frac(m_large_N2, l_large_N2)


    print 
    print '# Source classes #'
    print Ncat


    l_large = 'large (s>{s:.0f}")'.format(s=size_large)
    print_frac(m_large, l_large)
    
    l_large_bright = ' -large (s>{s:.0f}") & bright (S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large)
    print_frac(m_large_bright, l_large_bright)

    l_large_bright_nhuge_lr = '  --large (s>{s:.0f}") & bright (S>{f:.0f} mJy) & not huge & lr'.format(f=fluxcut, s=size_large)
    print_frac(m_large_bright_nhuge_lr, l_large_bright_nhuge_lr)

    l_large_faint = ' -large (s>{s:.0f}") & faint (S<{f:.0f} mJy)'.format(f=fluxcut, s=size_large)
    print_frac(m_large_faint, l_large_faint)


    l_small = 'compact (s<{s:.0f}")'.format(s=size_large)
    print_frac(m_small, l_small)
    
    
    l_small_isol = ' -compact isolated (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_isol, l_small_isol)

    l_small_isol_lr = '  --compact isolated good LR (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_isol_lr, l_small_isol_lr)

    l_small_isol_nlr ='  --compact isolated bad LR (s<{s2:.0f}", NN>{nn:.0f}")'.format(s2=size_large, nn=separation1)
    print_frac(m_small_isol_nlr, l_small_isol_nlr)

    l_small_nisol = ' -compact not isolated (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol, l_small_nisol)

    l_small_nisol_cluster = '  --compact not isolated (s<{s:.0f}", NN<{nn:.0f}") clustered 4 within 60"'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_cluster, l_small_nisol_cluster)
    
    l_small_nisol_NNlarge = '  --compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN large (s>{s:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNlarge, l_small_nisol_NNlarge)

    l_small_nisol_NNlarge_NNbright = '  --compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN large (s>{s:.0f}") NN bright'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNlarge_NNbright, l_small_nisol_NNlarge_NNbright)

    l_small_nisol_NNlarge_NNfaint = '  --compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN large (s>{s:.0f}") NN faint'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNlarge_NNfaint, l_small_nisol_NNlarge_NNfaint)

    l_small_nisol_NNsmall = '  --compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall, l_small_nisol_NNsmall)
    
    l_small_nisol_NNsmall_lr = '   ---compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}"), good LR'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall_lr, l_small_nisol_NNsmall_lr)
    
    
    l_small_nisol_NNsmall_lr_NNlr = '    ----compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}"), good LR, NN good LR'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall_lr_NNlr, l_small_nisol_NNsmall_lr_NNlr)
    
    l_small_nisol_NNsmall_lr_NNnlr = '    ----compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}"), good LR, NN bad LR'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall_lr_NNnlr, l_small_nisol_NNsmall_lr_NNnlr)
    
    
    l_small_nisol_NNsmall_nlr = '   ---compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}"), bad LR'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall_nlr, l_small_nisol_NNsmall_nlr)
    
    
    l_small_nisol_NNsmall_nlr_NNlr = '    ----compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}"), bad LR, NN good LR'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall_nlr_NNlr, l_small_nisol_NNsmall_nlr_NNlr)
    
    l_small_nisol_NNsmall_nlr_NNnlr = '    ----compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN compact (s<{s:.0f}"), bad LR, NN bad LR'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol_NNsmall_nlr_NNnlr, l_small_nisol_NNsmall_nlr_NNnlr)
    

    #l_small_nisol_match = '  --compact not isolated good LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    #print_frac(m_small_nisol_match, l_small_nisol_match)

    #l_small_nisol_NNlr = '  --compact not isolated neighbour good LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    #print_frac(m_small_nisol_NNlr, l_small_nisol_NNlr)
    
    #l_small_nisol_NNlr = '  --compact not isolated bad LR, but neighbour good LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    #print_frac(m_small_nisol_NNlr, l_small_nisol_NNlr)

    #l_small_nisol_nnmatch = '  --compact not isolated bad LR, and neighbour bad LR (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    #print_frac(m_small_nisol_nnmatch, l_small_nisol_nnmatch)


    # all visual sources

    print
    label = 'compact not isolated S (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol&m_S, label)

    label = 'compact not isolated !S (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1)
    print_frac(m_small_nisol&(~m_S), label)

# all groups
masks = [m_large,
m_large_bright,
m_large_bright_nhuge_lr, ##
m_large_faint,
m_small,
m_small_isol,
m_small_isol_lr,
m_small_isol_nlr,
m_small_nisol,
m_small_nisol_cluster, ##
m_small_nisol_NNlarge,
m_small_nisol_NNlarge_NNbright,
m_small_nisol_NNlarge_NNfaint,
m_small_nisol_NNsmall,
m_small_nisol_NNsmall_lr,
m_small_nisol_NNsmall_lr_NNlr,
m_small_nisol_NNsmall_lr_NNnlr,
m_small_nisol_NNsmall_nlr,
m_small_nisol_NNsmall_nlr_NNlr,
m_small_nisol_NNsmall_nlr_NNnlr]

names = ['m_large',
'm_large_bright',
'm_large_bright_nhuge_lr', ##
'm_large_faint',
'm_small',
'm_small_isol',
'm_small_isol_lr',
'm_small_isol_nlr',
'm_small_nisol',
'm_small_nisol_cluster', ##
'm_small_nisol_NNlarge',
'm_small_nisol_NNlarge_NNbright',
'm_small_nisol_NNlarge_NNfaint',
'm_small_nisol_NNsmall',
'm_small_nisol_NNsmall_lr',
'm_small_nisol_NNsmall_lr_NNlr',
'm_small_nisol_NNsmall_lr_NNnlr',
'm_small_nisol_NNsmall_nlr',
'm_small_nisol_NNsmall_nlr_NNlr',
'm_small_nisol_NNsmall_nlr_NNnlr']


#for mask,name in zip(masks,names):
    #print name,'=',1.0*np.sum(mask)/len(mask)
    ## select random 100
    
def fraction(mask):
    return 1.*np.sum(mask)/len(mask)

def num(mask):
    return int(fraction(mask)*Ncat)
    
    
# make sample files
if 0:
    # mutually exclusive groups
    masks = [m_large_bright,
    m_large_bright_nhuge_lr, ## test
    m_large_faint,
    m_small_isol_lr,
    m_small_isol_nlr,
    m_small_nisol_cluster, ## test
    m_small_nisol_NNlarge,  # not exclusive anymore
    m_small_nisol_NNlarge_NNbright,
    m_small_nisol_NNlarge_NNfaint,
    m_small_nisol_NNsmall_lr_NNlr,
    m_small_nisol_NNsmall_lr_NNnlr,
    m_small_nisol_NNsmall_nlr_NNlr,
    m_small_nisol_NNsmall_nlr_NNnlr]

    names = ['m_large_bright',
    'm_large_bright_nhuge_lr', ##
    'm_large_faint',
    'm_small_isol_lr',
    'm_small_isol_nlr',
    'm_small_nisol_NNlarge',
    'm_small_nisol_NNlarge_NNbright',
    'm_small_nisol_NNlarge_NNfaint',
    'm_small_nisol_cluster', ##
    'm_small_nisol_NNsmall_lr_NNlr',
    'm_small_nisol_NNsmall_lr_NNnlr',
    'm_small_nisol_NNsmall_nlr_NNlr',
    'm_small_nisol_NNsmall_nlr_NNnlr']

    for  mask,name in zip(masks,names):
        t = lofarcat[mask]
        t = t[np.random.choice(np.arange(len(t)), 100)]
        fitsname = 'sample_'+name+'.fits'
        if os.path.exists(fitsname):
            os.system('rm '+fitsname)
        t.write(fitsname)


# make flowchart
if 1:
    import pygraphviz as pgv

    PW = 60.


    A=pgv.AGraph(directed=True, strict=True)
    A.edge_attr['arrowhead']='none'
    A.node_attr['start']='south'
    A.node_attr['end']='north'
    A.node_attr['style']='filled'
    A.node_attr['fillcolor']='white'
    A.edge_attr['color']='gray'
    A.edge_attr['tailclip']='false'
    A.edge_attr['headclip']='false'
    A.graph_attr['outputorder'] = 'edgesfirst'
    #A.graph_attr['splines'] = 'ortho'  # orthogonal
    A.graph_attr['rankdir'] = 'TB'

    A.add_node('all', label='ALL\n{n:d}'.format(n=num(m_all)), shape='parallelogram') 
    A.add_node('large', label='Large?\n>{s:.0f}"'.format(s=size_large), shape='diamond')
    A.add_node('large_bright', label='Bright?\n>{S:.0f}mJy'.format(S=fluxcut), shape='diamond')
    A.add_node('small', label='Isolated?\n>{nn:.0f}"'.format(nn=separation1), shape='diamond')
    A.add_node('small_isol', label='LR?\n ', shape='diamond')
    A.add_node('small_nisol', label='NN Large?\n>{s:.0f}"'.format(s=size_large), shape='diamond')
    A.add_node('small_isol_NNsmall', label='LR?\n ', shape='diamond')
    A.add_node('small_isol_NNsmall_NNlr', label='NN LR?\n ', shape='diamond')
    A.add_node('small_isol_NNsmall_NNnlr', label='NN LR?\n ', shape='diamond')

    A.add_node('m_large_bright', label='VC\n{n:d}'.format(n=num(m_large_bright)), shape='parallelogram', color='green') 
    A.add_node('m_large_faint', label='?\n{n:d}'.format(n=num(m_large_faint)), shape='parallelogram', color='orange') 
    A.add_node('m_small_isol_lr', label='accept LR\n{n:d}'.format(n=num(m_small_isol_lr)), shape='parallelogram', color='blue') 
    A.add_node('m_small_isol_nlr', label='accept no LR\n{n:d}'.format(n=num(m_small_isol_nlr)), shape='parallelogram', color='red') 
    A.add_node('m_small_nisol_NNlarge', label='?\n{n:d}'.format(n=num(m_small_nisol_NNlarge)), shape='parallelogram', color='orange') 
    A.add_node('m_small_nisol_NNsmall_lr_NNlr', label='accept LR\n{n:d}'.format(n=num(m_small_nisol_NNsmall_lr_NNlr)), shape='parallelogram', color='blue')
    A.add_node('m_small_nisol_NNsmall_lr_NNnlr', label='?\n{n:d}'.format(n=num(m_small_nisol_NNsmall_lr_NNnlr)), shape='parallelogram', color='orange')
    A.add_node('m_small_nisol_NNsmall_nlr_NNnlr', label='?\n{n:d}'.format(n=num(m_small_nisol_NNsmall_nlr_NNnlr)), shape='parallelogram', color='red') 
    A.add_node('m_small_nisol_NNsmall_nlr_NNlr', label='?\n{n:d}'.format(n=num(m_small_nisol_NNsmall_nlr_NNlr)), shape='parallelogram', color='orange') 

    A.add_edge('all', 'large', label='Y', penwidth=fraction(m_all)*PW)
    A.add_edge('large', 'large_bright', label='Y', penwidth=fraction(m_large)*PW)
    A.add_edge('large', 'small', label='N', penwidth=fraction(m_small)*PW)
    A.add_edge('small', 'small_isol', label='Y', penwidth=fraction(m_small_isol)*PW)
    A.add_edge('small', 'small_nisol', label='N', penwidth=fraction(m_small_nisol)*PW)
    A.add_edge('small_nisol', 'small_isol_NNsmall', label='N', penwidth=fraction(m_small_nisol_NNsmall)*PW)
    A.add_edge('small_isol_NNsmall', 'small_isol_NNsmall_NNlr', label='N', penwidth=fraction(m_small_nisol_NNsmall_nlr)*PW)
    A.add_edge('small_isol_NNsmall', 'small_isol_NNsmall_NNnlr', label='Y', penwidth=fraction(m_small_nisol_NNsmall_lr)*PW)
    A.add_edge('small_isol_NNsmall_NNnlr', 'm_small_nisol_NNsmall_lr_NNlr', label='Y', penwidth=fraction(m_small_nisol_NNsmall_lr_NNlr)*PW)
    A.add_edge('small_isol_NNsmall_NNnlr', 'm_small_nisol_NNsmall_lr_NNnlr', label='N', penwidth=fraction(m_small_nisol_NNsmall_lr_NNnlr)*PW)
    A.add_edge('small_isol_NNsmall_NNlr', 'm_small_nisol_NNsmall_nlr_NNlr', label='Y', penwidth=fraction(m_small_nisol_NNsmall_nlr_NNlr)*PW)
    A.add_edge('small_isol_NNsmall_NNlr', 'm_small_nisol_NNsmall_nlr_NNnlr', label='N', penwidth=fraction(m_small_nisol_NNsmall_nlr_NNnlr)*PW)
    A.add_edge('small_nisol', 'm_small_nisol_NNlarge', label='Y', penwidth=fraction(m_small_nisol_NNlarge)*PW)
    A.add_edge('small_isol', 'm_small_isol_lr', label='Y', penwidth=fraction(m_small_isol_lr)*PW)
    A.add_edge('small_isol', 'm_small_isol_nlr', label='N', penwidth=fraction(m_small_isol_nlr)*PW)
    A.add_edge('large_bright', 'm_large_bright', label='Y', penwidth=fraction(m_large_bright)*PW)
    A.add_edge('large_bright', 'm_large_faint', label='N', penwidth=fraction(m_large_faint)*PW)


    # adjust a graph parameter
    #A.graph_attr['epsilon']='0.001'
    print(A.string()) # print dot file to standard output
    #Optional prog=['neato'|'dot'|'twopi'|'circo'|'fdp'|'nop']
    #neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps, sccmap, tred, sfdp.
    A.layout('dot') # layout with dot
    A.draw('flow_s{s:.0f}_nn{nn:.0f}.png'.format(s=size_large,nn=separation1)) # write to file


if 0:
    # select N samples across LR space
    tt = psmlcat[m_small_nisol]
    tt.sort('lr_pc_7th')
    N = 100
    Ntt = len(tt)
    tt = tt[0:Ntt:Ntt/N]
    tt.write('sample_small_nisol.fits')
    tt = psmlcat[m_large&m_N1]
    tt.sort('lr_pc_7th')
    N = 100
    Ntt = len(tt)
    tt = tt[0:Ntt:Ntt/N]
    tt.write('sample_large_N1.fits')

def frac(n):
    global Ncat
    return 1.*n/Ncat

def fracN(f):
    global Ncat
    return f*Ncat

def convert_ax_n_to_frac(ax_n):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax_n.get_ylim()
    ax_f.set_ylim(frac(y1), frac(y2))
    ax_f.figure.canvas.draw()

def convert_ax_f_to_N(ax_f):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax_f.get_ylim()
    ax_n.set_ylim(fracN(y1), fracN(y2))
    ax_n.figure.canvas.draw()

def stack(ax, Nstep, masks, colors, alphas):
    n1 = 0
    n2 = 0
    for i in range(Nstep):
        n2 = n2+np.sum(masks[i])
        ax.fill_between([Nstep,Nstep+1], n1, n2, color=colors[i], alpha=alphas[i])
        n1 = n2
    n2 = Ncat
    ax.fill_between([Nstep,Nstep+1], n1, n2, color='C7', alpha=1.)
    return

f,ax_f = pp.paper_single_ax(TW=8, AR=0.7)
plt.subplots_adjust(right=0.85, bottom=0.3)
ax_n = ax_f.twinx()
ax_n.minorticks_on()
ax_f.callbacks.connect("ylim_changed", convert_ax_f_to_N)
ax_f.set_ylabel('fraction')
ax_n.set_ylabel('number')
ax_f.set_ylim(0,1)


# C0 = blue
# C3 = red
# C2 = green
# C7 - gray
masks = [m_large_bright, 
         m_small_isol_lr, 
         m_small_isol_nlr,
         m_small_nisol_match, 
         m_small_nisol_nomatch]
colors = ['C0',
          'C2',
          'C3',
          'C2',
          'C3']
alphas = [1,
          1,
          1,
          0.5,
          0.5]
labels = ['Large \& bright',
          'small \& isolated \& LR',
          'small \& isolated \& no LR',
          'small \& not isolated \& LR',
          'small \& not isolated \& no LR']

for i in range(len(masks)+1):
    stack(ax_n, i, masks[:i], colors[:i], alphas[:i])

#ax_f.xaxis.set_visible(False)
plt.xticks(range(1,len(masks)+2), labels, rotation=90)
ax_f.set_xlim(0,i+1)

plt.savefig('pipeline.png')

## check gaus ML


#for i in range(len(lofarcat)):
    #print lofarcat['SID'][i], np.array(lofargcat['SID'][lofarcat['G_ind'][i]])

#sys.exit()

if add_G:
    lofarmcat = lofarcat[m_M]
    sepmaxes = np.zeros(len(lofarmcat))
    classes = np.zeros(len(lofarmcat))
    for i in range(len(lofarmcat)):
        lr = lofarmcat['LR'][i] 
        alr  = np.log10(1+lr) >= lLR_thresh
        c = ac.SkyCoord(lofargcat[lofarmcat['G_ind'][i]]['RA'], lofargcat[lofarmcat['G_ind'][i]]['DEC'], unit='deg')
        sepmax = 0
        #print np.array(lofargcat[lofarmcat['G_ind'][i]]['RA'])
        for ci in c:
            #print ci.separation(c).to('arcsec')
            sepmax = np.max((sepmax, np.max(ci.separation(c).to('arcsec').value)))
        
        glr = np.array(lofargcat[lofarmcat['G_ind'][i]]['LR'] )
        aglr  = np.log10(1+glr) >= lLR_thresh
        #print lr, glr
        if alr: 
            if np.any(aglr):
                classes[i] = 1
            else:
                # accept source match
                classes[i] = 2
        elif ~alr:
            if np.any(aglr):
                classes[i] = 3
            else:
                classes[i] = 4
                # accept no source match and no gaus match
        sepmaxes[i] = sepmax
        #print alr, aglr, sepmax, lofarmcat['Maj'][i]

    f,axs = plt.subplots(2,2,sharex=True, sharey=True)
    axs = axs.flatten()
    for ic,lab in [[1,'A LR ; G LR'],[2,'A LR ; G !LR'],[3,'A !LR ; G LR'],[4,'A !LR ; G !LR']]:
        ax = axs[ic-1]
        ax.plot(lofarmcat['Maj'][classes==ic], sepmaxes[classes==ic], '.', label=lab)
        ax.legend()
        ax.set_ylabel('max G separation [arcsec]')
        ax.set_xlabel('size [arcsec]')
    plt.savefig('gaus_size_separation')


fluxcuts = np.logspace(-4, 0, 1000)
nS_fluxcuts = np.nan*np.zeros(len(fluxcuts))
for fi,fluxcut in enumerate(fluxcuts):
    m = lofarcat['Total_flux']/1e3 > fluxcut
    nS_fluxcuts[fi] = 1.*np.sum(lofarcat['S_Code'][m] == 'S') /np.sum(m)
    #nS_fluxcuts[fi] = 1.*np.sum(m)
f,ax = pp.paper_single_ax()
ax.plot(fluxcuts, nS_fluxcuts)
ax.set_ylabel('f(Single) ($S>S_{cut}$)')
ax.set_xlabel('$\log S_{cut}$ [Jy]')
plt.savefig('fraction_single_vs_S')


sizecuts = np.linspace(15, 60, 10)
fluxcuts = np.logspace(-3, 1, 1000)
f,ax = pp.paper_single_ax()
for si,sizecut in enumerate(sizecuts):
    ms = lofarcat['Maj'] > sizecut
    nS_fluxcuts = np.nan*np.zeros(len(fluxcuts))
    for fi,fluxcut in enumerate(fluxcuts):
        m = ms & (lofarcat['Total_flux']/1e3 > fluxcut)
        nS_fluxcuts[fi] = 1.*np.sum(m) /np.sum(ms)
    ax.plot(fluxcuts, nS_fluxcuts)
#ax.set_ylabel('$f(Maj>Maj_{cut})$ ($S>S_{cut}$)')
#ax.set_xlabel('$\log S_{cut}$ [Jy]')
plt.savefig('fraction_large_vs_S')



sizecuts = np.arange(10, 35, 1)
NNcuts = np.arange(20, 125, 5)
IM = np.zeros((len(sizecuts), len(NNcuts)))
fluxcuts = np.logspace(-3, 1, 1000)
f,ax = pp.paper_single_ax()
for si,sizecut in enumerate(sizecuts):
    for ni,NNcut in enumerate(NNcuts):
        m = (lofarcat['Maj'] <= sizecut) & (lofarcat['NN_sep'] >= NNcut)
        IM[si,ni] = np.sum(m)
IM = IM/Ncat
c = ax.imshow(IM.T, origin='lower', extent=(10,60, 20,120))
cbar = plt.colorbar(c)
cbar.set_label('fraction')
ax.invert_xaxis()
ax.set_xlabel(r'$<$ size [arcsec]')
ax.set_ylabel(r'$>$ NN separation [arcsec]')
plt.savefig('number_compact_isolated')

f,axs = plt.subplots(1,2,sharex=False,sharey=True,figsize=(12,6))
ax=axs[0]
ax.plot(NNcuts,IM.T)
ax.set_ylabel('fraction')
ax.set_xlabel(r'$>$ NN separation [arcsec]')
ax=axs[1]
ax.plot(sizecuts,IM)
ax.set_xlabel(r'$<$ size [arcsec]')


nb=100
# plot LR distribuion for different classes
f,ax = pp.paper_single_ax()
_ =ax.hist(np.log10(1.+lofarcat['LR']), bins=100, normed=True, log=False,histtype='step',color='k',linewidth=2,label='All')
_ =ax.hist(np.log10(1.+lofarcat['LR'][m_small_isol]), bins=100, normed=True, histtype='step', label=l_small_isol)
#_ =ax.hist(np.log10(1.+lofarcat['LR'][m_small_isol_nS]), bins=100, normed=True, histtype='step', label=l_small_isol_nS)
_ =ax.hist(np.log10(1.+lofarcat['LR'][m_small_nisol]), bins=100, normed=True, histtype='step', label=l_small_nisol)
_ =ax.hist(np.log10(1.+lofarcat['LR'][m_large]), bins=100, normed=True, histtype='step', label=l_large)
ax.legend()
ax.set_ylim(0,2)
ax.set_xlabel('$\log (1+LR)$')
ax.set_ylabel('$N$')
plt.savefig('lr_dist_classes')

# plot LR distribuion for different classes
f,ax = pp.paper_single_ax()
counts, xedges, yedges, im =ax.hist2d(np.log10(1.+lofarcat['LR']), np.log10(lofarcat['Maj']), bins=100, normed=True, vmin=0, vmax=2, label='')
cbar = plt.colorbar(im, ax=ax)
ax.legend()
#ax.set_ylim(0,2)
ax.set_xlabel('$\log (1+LR)$')
ax.set_ylabel('$\log$ Maj [arcsec]')
cbar.set_label('$N$')
plt.savefig('lr_dist_size')

f,ax = pp.paper_single_ax()
counts, xedges, yedges, im =ax.hist2d(np.log10(1.+lofarcat['LR']), np.log10(lofarcat['Total_flux']), bins=100, normed=True, vmin=0, vmax=2, label='')
cbar = plt.colorbar(im, ax=ax)
ax.legend()
#ax.set_ylim(0,2)
ax.set_xlabel('$\log (1+LR)$')
ax.set_ylabel('$\log S$ [Jy]')
cbar.set_label('$N$')
plt.savefig('lr_dist_flux')



#f,ax = pp.paper_single_ax()
#f = plt.figure()
f,axs = plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,6))
ax = axs[0]
counts, xedges, yedges, im =ax.hist2d(np.log10(lofarcat['Maj'][m_S]), np.log10(lofarcat['Total_flux'][m_S]), bins=100, label='')
cbar = plt.colorbar(im, ax=ax)
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
ax.vlines(np.log10(15.),y1,y2)
ax.hlines(np.log10(10.),x1,x2)
ax.legend()
ax.set_title('S')
#ax.set_ylim(0,2)
ax.set_xlabel('$\log $ Maj [arcsec]')
ax.set_ylabel('$\log S$ [mJy]')
cbar.set_label('$N$')
ax = axs[1]
counts, xedges, yedges, im =ax.hist2d(np.log10(lofarcat['Maj'][~m_S]), np.log10(lofarcat['Total_flux'][~m_S]), bins=100, label='')
cbar = plt.colorbar(im, ax=ax)
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
ax.vlines(np.log10(15.),y1,y2)
ax.hlines(np.log10(10.),x1,x2)
ax.legend()
ax.set_title('!S')
#ax.set_ylim(0,2)
ax.set_xlabel('$\log $ Maj [arcsec]')
#ax.set_ylabel('$\log S$ [mJy]')
cbar.set_label('$N$')
plt.savefig('lr_dist_size_flux')




#f,ax = pp.paper_single_ax()
#f = plt.figure()
f,axs = plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,6))
ax = axs[0]
counts, xedges, yedges, im =ax.hist2d((lofarcat['Maj'][m_lrgood]), (lofarcat['NN_sep'][m_lrgood]), bins=200, range=((0,50),(0,200)), label='')
cbar = plt.colorbar(im, ax=ax)
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
ax.vlines((15.),y1,y2)
ax.hlines((10.),x1,x2)
ax.legend()
ax.set_title('good LR')
#ax.set_ylim(0,2)
ax.set_xlabel('Maj [arcsec]')
ax.set_ylabel('NN separation [arcsec]')
cbar.set_label('$N$')
ax = axs[1]
counts, xedges, yedges, im =ax.hist2d((lofarcat['Maj'][~m_lrgood]), (lofarcat['NN_sep'][~m_lrgood]), bins=200, range=((0,50),(0,200)), label='')
cbar = plt.colorbar(im, ax=ax)
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
ax.vlines((15.),y1,y2)
ax.hlines((10.),x1,x2)
ax.legend()
ax.set_title('bad LR')
#ax.set_ylim(0,2)
ax.set_xlabel('Maj [arcsec]')
#ax.set_ylabel('$\log S$ [mJy]')
cbar.set_label('$N$')
plt.savefig('lr_dist_size_nnsep')





# # diagnostic plots 
# 

# In[7]:


# plot size distribution
f, ax = pp.paper_single_ax()
ax.hist(lofarcat['Maj'][~maskDC0], range=(0,80), bins=100, histtype='step', label='All')
ax.hist(lofarcat['Maj'][~maskDC0&m_S], range=(0,80), bins=100, histtype='step', label='S')
ax.hist(lofarcat['Maj'][~maskDC0&m_M], range=(0,80), bins=100, histtype='step', label='M')
ax.hist(lofarcat['Maj'][~maskDC0&m_C], range=(0,80), bins=100, histtype='step', label='C')
ax.set_xlabel('Major Axis [arcsec]')
ax.set_ylabel('N')
ax.legend()
plt.savefig('size_dist_classes')





# In[8]:


# plot nearest neighbour distribution
f,ax = pp.paper_single_ax()
ax.hist(f_nn_sep2d.to('arcsec').value, bins=100, histtype='step', label='All')
ax.hist(f_nn_sep2d.to('arcsec').value[m_S], bins=100, histtype='step', label='S')
ax.set_xlabel('Nearest source [arcsec]')
ax.set_ylabel('N')
ax.legend()
plt.savefig('NNdist_dist')


# In[9]:



# 2D histogram : size-nearest neighbour distance
# for 'S' sources
f,ax = pp.paper_single_ax()
X =  f_nn_sep2d.to('arcsec').value[~maskDC0&m_S]
Y = lofarcat['Maj'][~maskDC0&m_S]
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
ax.hlines(size_large,x1,x2,colors='k',linestyle='dashed')
ax.vlines(separation1,y1,y2,colors='k',linestyle='dashed')
ax.set_xlabel('NN separation [arcsec]')
ax.set_ylabel('DCmaj [arcsec]')
ax.contour(xc, yc, H2)
plt.savefig('size_NNdist_dist_s')



# In[10]:


# and 'M' sources
f,ax = pp.paper_single_ax()
X =  f_nn_sep2d.to('arcsec').value[~maskDC0&m_M]
Y = lofarcat['Maj'][~maskDC0&m_M]
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
ax.hlines(size_large,x1,x2,colors='k',linestyle='dashed')
ax.vlines(separation1,y1,y2,colors='k',linestyle='dashed')
ax.set_xlabel('NN separation [arcsec]')
ax.set_ylabel('DCmaj [arcsec]')
ax.contour(xc, yc, H2)
plt.savefig('size_NNdist_dist_m')



# In[ ]:



