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


class Mask:
    '''Mask store a boolean mask and associated information necessary for the flowchart
    '''
    
    def __init__(self, mask, label, trait, level=0, verbose=True, masterlist=None, qlabel=None, color=None):
        self.mask = mask
        self.label = label
        if qlabel is not None :
            self.qlabel = qlabel
        else:
            self.qlabel = label
        if isinstance(trait,str):
            self.traits = list([trait])
        else:
            self.traits = list(trait)
            
        self.name = '_'.join(self.traits)
            
        self.color = color
        self.level = level
        
        self.N = self.total()
        self.n = self.msum()
        self.f = self.fraction()
        self.p = self.percent()
        
        self.has_children = False
        self.has_parent = False
        
        self.Nchildren = 0
        self.children = None
        self.parent = None
        
        if masterlist is not None:
            masterlist.append(self)
        
        if verbose:
            self.print_frac()
        
        return
    
    def percent(self):
        return 100.*np.sum(self.mask)/self.N
    
    def fraction(self):
        return 1.*np.sum(self.mask)/self.N

    def msum(self):
        return np.sum(self.mask)
    
    def total(self):
        return len(self.mask)
    
    def print_frac(self, vformat=True):
        '''vformat = True will print with formatted spaces indicative of the hierarchical structure
        '''
        if vformat and self.level > 0:
            vv = ' '*self.level + '-'*self.level
        else:
            vv = ' '
        print '{n:6d} ({f:5.1f}%){vv:s}{label:s}'.format(vv=vv, n=self.n, f=self.p, label=self.label)
        
    def __str__(self):
        return self.name
        
    def submask(self, joinmask, label, newtrait, edgelabel='Y', verbose=True, qlabel=None, masterlist=None, color=None):
        '''create a new submask based on this instance -- join masks with AND
        # qlabel  is the question that will be asked
        # edgelabel is the answer to the question asked to get here
        '''
        newmask = self.mask & joinmask
        newtraits = list(self.traits)  # copy list of traits - lists are mutable!!
        newtraits.append(newtrait)     # append new trait onto copy
        newlevel = self.level + 1
        
        childmask = Mask(newmask, label, newtraits, level=newlevel, masterlist=masterlist, verbose=verbose, qlabel=qlabel, color=color)  
        
        childmask.has_parent = True
        childmask.parent = self
        
        childmask.edgelabel = edgelabel  
        
        if not self.has_children:
            self.has_children = True
            self.children = [childmask]
            self.Nchildren = 1
        else:
            newchildren = list(self.children)  # copy list of traits - lists are mutable!!
            newchildren.append(childmask)
            self.children = newchildren
            self.Nchildren = len(newchildren)
            
        return childmask
    
    # make sample files
    def make_sample(self, cat, Nsample=250):
        '''create a random subsample of the masked catalogue 'cat'
        '''
        
        t = cat[self.mask]
        if Nsample is None:
            Nsample = len(t)
        Nsample = np.min((Nsample, len(t)))
        if Nsample ==0 : return
        t = t[np.random.choice(np.arange(len(t)), Nsample)]
        fitsname = 'sample_'+self.name+'.fits'
        if os.path.exists(fitsname):
            os.remove(fitsname)
        t.write(fitsname)
        
        return
    
    def is_disjoint(self, othermask):
        
        assert isinstance(othermask, Mask), 'need to compare to another Mask instance'
        
        if np.sum((self.mask) & (othermask.mask)) == 0:
            return True
        else:
            return False
        
        return
        
def Masks_disjoint_complete(masklist):
    '''test whether a list of masks is disjoint and complete
    '''
    Z = np.zeros(len(masklist[0].mask), dtype=bool)
    O = np.ones(len(masklist[0].mask), dtype=bool)
    for t in masklist:
        Z = Z & t.mask
        O = O | t.mask
    
    return np.all(O) and np.all(~Z)



### Required INPUTS
# lofar source catalogue, gaussian catalogue and ML catalogues for each

#path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/source_class/t1_dr1/'
#lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.1.gaus.fits'
#lofarcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.1.srl.fits'
#psmlcat_file = path+'lofar_matched_all.fix.fits'
#psmlgcat_file = path+'lofar_matched_gaus.fits'


path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.gaus.fits'
lofarcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fits'
psmlcat_file = path+'lofar_pw.fits'
psmlgcat_file = path+'lofar_gaus_pw.fits'




# Gaus catalogue
lofargcat = Table.read(lofargcat_file)
# only relevant gaussians are in M or C sources
lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

# Source catalogue
lofarcat = Table.read(lofarcat_file)

# PS ML - matches for sources and gaussians
psmlcat = Table.read(psmlcat_file)
psmlgcat = Table.read(psmlgcat_file)



## match the gaussians to the sources

## quicker to generate new unique names than match on 2 columns
## get new unique source_id by combining mosaic and src id
## replace string mosaic ID with unique int (perhaps there is a more logical mapping of mosaic name to int value)
#mid = lofargcat['Mosaic_ID']
#mid_unique = np.unique(mid)
#mid_int = np.array([np.where(mid_unique==m)[0][0] for m in mid])
## combine with Source_id for unique ID
#g_src_id_new =   10000*mid_int + lofargcat['Source_Name']
#lofargcat.add_column(Column(g_src_id_new, 'SID'))

#mid = lofarcat['Mosaic_ID']
#mid_unique = np.unique(mid)
#mid_int = np.array([np.where(mid_unique==m)[0][0] for m in mid])
## combine with Source_id for unique ID
#src_id_new =   10000*mid_int + lofarcat['Source_Name']
#lofarcat.add_column(Column(src_id_new, 'SID'))



   

## get the panstarrs ML information

# join the ps ml cat  - they have identical RA/DEC (source_names were wrong)
c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
cpsml = ac.SkyCoord(psmlcat['RA'], psmlcat['DEC'], unit="deg")
f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,cpsml,nthneighbor=1)

#psmlcat = psmlcat[f_nn_idx][f_nn_sep2d==0]
#lofarcat = lofarcat[f_nn_sep2d==0]

# note the large sources are missing from the ML catalogue
lrcol = np.zeros(len(lofarcat),dtype=float)
lrcol[f_nn_sep2d==0] = psmlcat['lr'][f_nn_idx][f_nn_sep2d==0]

#lofarcat.add_column(Column(psmlcat['lr_pc_7th'], 'LR'))
lofarcat.add_column(Column(lrcol, 'cLR'))
lrcol[np.isnan(lrcol)] = 0
lofarcat.add_column(Column(lrcol, 'LR'))


# join the ps ml gaus cat  - they have identical RA/DEC (source_names were wrong)
cg = ac.SkyCoord(lofargcat['RA'], lofargcat['DEC'], unit="deg")
cpsmlg = ac.SkyCoord(psmlgcat['RA'], psmlgcat['DEC'], unit="deg")
f_nn_idx_g,f_nn_sep2d_g,f_nn_dist3d_g = ac.match_coordinates_sky(cg,cpsmlg,nthneighbor=1)

# note the large sources are missing from the ML catalogue
lrgcol = np.zeros(len(lofargcat),dtype=float)
lrgcol[f_nn_sep2d_g==0] = psmlgcat['lr'][f_nn_idx_g][f_nn_sep2d_g==0]

#lofarcat.add_column(Column(psmlcat['lr_pc_7th'], 'LR'))
lofargcat.add_column(Column(lrgcol, 'LR'))

add_G = False   # add the gaussian information
lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=int), 'Ng'))
if add_G:
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=list), 'G_ind'))

m_S = lofarcat['S_Code'] =='S'
minds = np.where(~m_S)[0]
for i,sid in zip(minds, lofarcat['Source_Name'][~m_S]):
    ig = np.where(lofargcat['Source_Name']==sid)[0]
    lofarcat['Ng'][i]= len(ig)
    
    if add_G:
        lofarcat['G_ind'][i]= ig



## get 2MASX information
xsc_file = path+'2MASX_hetdex.fits'
xsc = Table.read(xsc_file)
xsc['r_ext'][xsc['designation']=='13174820+4702571    '] = 20. # fix known issue with one source - has bad size, use sdss size!
xsc['r_ext'][xsc['designation']=='11244075+5218078    '] = 20. # fix known issue with one source - has bad size, use sdss size!
xsc['r_ext'][xsc['designation']=='11435865+5600295    '] = 12. # fix known issue with one source - has bad size, use sdss size!
xsc['r_ext'][xsc['designation']=='11095404+4859120    '] = 13. # fix known issue with one source - has bad size, use sdss size!
xsc['r_ext'][xsc['designation']=='11393585+5555286    '] = 13. # fix known issue with one source - has bad size, use sdss size!


c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
cxsc = ac.SkyCoord(xsc['ra'], xsc['dec'], unit="deg")
f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,cxsc,nthneighbor=1)

xsc_nn = xsc[f_nn_idx]




def accept_match_2mass(mask, lcat, xcat, plot=False, selname=None):
    from matplotlib import patches

    
    if plot:
        f,ax = pp.paper_single_ax()
    if selname is not None:
        nn = lcat['Source_Name']==selname
        mask = mask & nn
        
    idx = np.arange(len(lcat))
    iinds = idx[mask]
    inellipse = np.zeros(len(lcat ), dtype=bool)
    for i in iinds:
        
        tl = lcat[i]
        tx = xcat[i]
            
        # assumning flat sky here...
        g_ell_center = (tx['ra'], tx['dec'])
        r_a = (tx['r_ext']+ tl['Maj'] )  / 3600. #to deg
        r_b = r_a*tx['k_ba']
        angle = 90.-tx['k_phi']  #(anticlockwise from x-axis)
        rangle = angle *np.pi/180.

        cos_angle = np.cos(np.radians(180.-angle))
        sin_angle = np.sin(np.radians(180.-angle))

        xc = tl['RA'] - g_ell_center[0]
        yc = tl['DEC'] - g_ell_center[1]

        xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle 

        rad_cc = (xct/r_a)**2. + (yct/r_b)**2.
        
        if rad_cc <= 1:
            inellipse[i] = 1
            
        if plot:
            g_ellipse = patches.Ellipse(g_ell_center, 2*r_a, 2*r_b, angle=angle, fill=False, ec='green', linewidth=2)

            ax.plot(g_ell_center[0], g_ell_center[1], 'g.')
            ax.plot(tl['RA'], tl['DEC'], 'k.')
            ax.add_patch(g_ellipse)
            
            l_ellipse = patches.Ellipse((tl['RA'], tl['DEC']), 2*tl['Maj']/ 3600., 2*tl['Min']/ 3600., angle=90.-tl['PA'], fill=False, ec='blue', linewidth=2)
            ax.add_patch(l_ellipse)
            
            if rad_cc <= 1:
                ax.plot(tl['RA'], tl['DEC'], 'r+')
            
            mell = g_ellipse.contains_point((tl['RA'], tl['DEC']), radius=tl['Maj']/3600)
            if mell:
                ax.plot(tl['RA'], tl['DEC'], 'gx')
            
            ax.plot(g_ell_center[0], g_ell_center[1], 'g+')
            ax.plot(tl['RA'], tl['DEC'], 'b.')
    if plot:
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.invert_xaxis()
        ax.axis('equal')
        
    return inellipse
        
xmatch0 = f_nn_sep2d.value*u.deg < np.array(xsc_nn['r_ext'])*u.arcsec
xmatch1 = f_nn_sep2d.value*u.deg < np.array(xsc_nn['r_ext'] + lofarcat['Maj'])*u.arcsec

inellipse = accept_match_2mass(xmatch1, lofarcat, xsc_nn)
xmatch = xmatch1 & inellipse

Xhuge =  xmatch & (xsc_nn['r_ext'] >= 240.)
XLarge =  xmatch & (xsc_nn['r_ext'] >= 60.) & (xsc_nn['r_ext'] < 240.)
Xlarge =  xmatch & (xsc_nn['r_ext'] >= 20.) & (xsc_nn['r_ext'] < 60.)
Xsmall =  xmatch0 & (xsc_nn['r_ext'] >= 0.) & (xsc_nn['r_ext'] < 20.)

if '2MASX' not in lofarcat.colnames:
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=bool),'2MASX'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype='S20'),'2MASX_name'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=float),'2MASX_ra'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=float),'2MASX_dec'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=float),'2MASX_size'))
for m in [Xhuge, XLarge, Xlarge, Xsmall]:
    lofarcat['2MASX'][m]  = m[m]
    lofarcat['2MASX_name'][m]  = xsc_nn['designation'][m]
    lofarcat['2MASX_ra'][m]  = xsc_nn['ra'][m]
    lofarcat['2MASX_dec'][m]  = xsc_nn['dec'][m]
    lofarcat['2MASX_size'][m]  = xsc_nn['r_ext'][m]


#t=M_all.submask(huge , '2MASX_huge', '2MASX_huge')
#t.make_sample(lofarcat)
#t=M_all.submask(Large , '2MASX_Large', '2MASX_Large')
#t.make_sample(lofarcat)
#t=M_all.submask(large , '2MASX_large', '2MASX_large')
#t.make_sample(lofarcat)
#t=M_all.submask(small , '2MASX_small', '2MASX_small')
#t.make_sample(lofarcat)

#accept_match_2mass(xmatch1, lofarcat, xsc_nn, selname='ILTJ113935.922+555529.21', plot=True)


## get artefact information

# for now, no artefacts
artefact = np.zeros(len(lofarcat),dtype=bool)
if 'aretefact' not in lofarcat.colnames:
    lofarcat.add_column(Column(artefact,'artefact'))
    

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
   
artefact = np.zeros(len(lofarcat),dtype=bool)
selind1 = np.where((lofarcat['Maj'] < 8.) & (lofarcat['Total_flux'] > 100.))[0]
c = ac.SkyCoord(lofarcat['RA'][selind1], lofarcat['DEC'][selind1], unit="deg")
    

#############################################################################


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


########################################################


# make samples

# # source classes
# 
# clases from draft flowchart
# 
# source classes - parameters & masks

# >15 " and 10mJY -2%

size_large = 15.           # in arcsec
separation1 = 60.          # in arcsec
size_huge = 25.            # in arcsec
#separation1 = 30.          # in arcsec
lLR_thresh = 1.            # LR threshold
fluxcut = 10               # in mJy

Ncat = len(lofarcat)

#m_all = lofarcat['RA'] > -1

masterlist = []

M_all = Mask(lofarcat['RA'] > -1,
                'All',
                'all',
                qlabel='Large?\n(s>{s:.0f}")'.format(s=size_large),
                masterlist=masterlist)


# large 
M_large = M_all.submask(lofarcat['Maj'] > size_large,
                    'large (s>{s:.0f}")'.format(s=size_large),
                    'large',
                    qlabel='Bright?\n(S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                    masterlist=masterlist)

# large bright
M_large_bright = M_large.submask(lofarcat['Total_flux'] > fluxcut,
                    'large (s>{s:.0f}") & bright (S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                    'bright',
                    qlabel='Not huge?\nLR?',
                    color='green',
                    masterlist=masterlist)

# large bright not huge with lr
M_large_bright_nhuge_lr = M_large_bright.submask((lofarcat['Maj'] <= size_huge) & (np.log10(1+lofarcat['LR']) <= lLR_thresh),
                    'large (s>{s:.0f}") & bright (S>{f:.0f} mJy) & not huge (s<={s2:.0f}") & lr'.format(f=fluxcut, s=size_large, s2=size_huge),
                    'nhuge_lr',
                    qlabel='accept LR??',
                    color='blue',
                    masterlist=masterlist)

# large bright not huge with lr
M_large_bright_huge = M_large_bright.submask(~((lofarcat['Maj'] <= size_huge) & (np.log10(1+lofarcat['LR']) <= lLR_thresh)),
                    'large (s>{s:.0f}") & bright (S>{f:.0f} mJy) & not(not huge (s<={s2:.0f}") & lr)'.format(f=fluxcut, s=size_large, s2=size_huge),
                    'huge',
                    qlabel='VC',
                    color='green',
                    edgelabel='N',
                    masterlist=masterlist)

# large faint
M_large_faint = M_large.submask(lofarcat['Total_flux'] <= fluxcut,
                    'large (s>{s:.0f}") & faint (S<={f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                    'faint',
                    edgelabel='N',
                    qlabel='TBC?',
                    color='orange',
                    masterlist=masterlist)

# compact 
M_small = M_all.submask(lofarcat['Maj'] <= size_large,
                    'compact (s<{s:.0f}")'.format(s=size_large),
                    'small',
                    edgelabel='N',
                    qlabel='Isolated?\n(NN>{nn:.0f}")'.format(nn=separation1),
                    masterlist=masterlist)
# compact isolated
M_small_isol = M_small.submask(lofarcat['NN_sep'] > separation1,
                    'compact isolated (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1),
                    'isol',
                    qlabel='S?',
                    masterlist=masterlist)


# compact isolated
M_small_isol_S = M_small_isol.submask(lofarcat['S_Code'] == 'S',
                    'compact isolated (s<{s:.0f}", NN>{nn:.0f}") S'.format(s=size_large, nn=separation1),
                    'S',
                    qlabel='LR?',
                    masterlist=masterlist)


# compact isolated good lr
M_small_isol_S_lr = M_small_isol_S.submask(np.log10(1+lofarcat['LR']) > lLR_thresh,
                    'compact isolated good LR (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1),
                    'lr',
                    color='blue',
                    qlabel='Accept LR',
                    masterlist=masterlist)

# compact isolated badd lr
M_small_isol_S_nlr = M_small_isol_S.submask(np.log10(1+lofarcat['LR']) <= lLR_thresh,
                    'compact isolated bad LR (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1),
                    'nlr',
                    edgelabel='N',
                    color='red',
                    qlabel='Accept no LR',
                    masterlist=masterlist)

# compact isolated
M_small_isol_nS = M_small_isol.submask(lofarcat['S_Code'] != 'S',
                    'compact isolated (s<{s:.0f}", NN>{nn:.0f}") !S'.format(s=size_large, nn=separation1),
                    'nS',
                    edgelabel='N',
                    color='orange',
                    qlabel='TBC?',
                    masterlist=masterlist)


# compact not isolated
M_small_nisol = M_small.submask(lofarcat['NN_sep'] <= separation1,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1),
                    'nisol',
                    edgelabel='N',
                    qlabel='NN Large?',
                    masterlist=masterlist)

# compact not isolated, nnlarge
M_small_nisol_NNlarge = M_small_nisol.submask(lofarcat['NN_Maj'] > size_large,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN large (s>{s:.0f}")'.format(s=size_large, nn=separation1),
                    'NNlarge',
                    edgelabel='Y',
                    qlabel='NN Bright?',
                    masterlist=masterlist)


# compact not isolated, nnlarge
M_small_nisol_NNlarge_NNbright = M_small_nisol_NNlarge.submask(lofarcat['NN_Total_flux'] > fluxcut,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN large (s>{s:.0f}") NN bright (S>{f:.0f} mJy)'.format(s=size_large, nn=separation1, f=fluxcut),
                    'NNbright',
                    edgelabel='Y',
                    color='green',
                    qlabel='NN in VC list',
                    masterlist=masterlist)

# compact not isolated, nnlarge
M_small_nisol_NNlarge_NNfaint = M_small_nisol_NNlarge.submask(lofarcat['NN_Total_flux'] <= fluxcut,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN large (s>{s:.0f}") NN faint (S<={f:.0f} mJy)'.format(s=size_large, nn=separation1, f=fluxcut),
                    'NNfaint',
                    edgelabel='N',
                    color='orange',
                    qlabel='TBC?',
                    masterlist=masterlist)


# compact not isolated, nnsmall
M_small_nisol_NNsmall = M_small_nisol.submask(lofarcat['NN_Maj'] <= size_large,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}")'.format(s=size_large, nn=separation1),
                    'NNsmall',
                    edgelabel='N',
                    qlabel='LR?',
                    masterlist=masterlist)

# compact not isolated, nnsmall, lr
M_small_nisol_NNsmall_lr = M_small_nisol_NNsmall.submask(np.log10(1+lofarcat['LR']) > lLR_thresh,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR'.format(s=size_large, nn=separation1),
                    'lr',
                    edgelabel='Y',
                    qlabel='NN LR?',
                    masterlist=masterlist)

# compact not isolated, nnsmall, lr, NNlr
M_small_nisol_NNsmall_lr_NNlr = M_small_nisol_NNsmall_lr.submask(np.log10(1+lofarcat['NN_LR']) > lLR_thresh,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR, NN good lr'.format(s=size_large, nn=separation1),
                    'NNlr',
                    edgelabel='Y',
                    color='blue',
                    qlabel='accept LR?',
                    masterlist=masterlist)

# compact not isolated, nnsmall, lr, NNnlr
M_small_nisol_NNsmall_lr_NNnlr = M_small_nisol_NNsmall_lr.submask(np.log10(1+lofarcat['NN_LR']) <= lLR_thresh,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR, NN bad lr'.format(s=size_large, nn=separation1),
                    'NNnlr',
                    edgelabel='N',
                    color='orange',
                    qlabel='TBC?',
                    masterlist=masterlist)

# compact not isolated, nnsmall, nlr
M_small_nisol_NNsmall_nlr = M_small_nisol_NNsmall.submask(np.log10(1+lofarcat['LR']) <= lLR_thresh,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR'.format(s=size_large, nn=separation1),
                    'nlr',
                    edgelabel='N',
                    qlabel='NN LR?',
                    masterlist=masterlist)

# compact not isolated, nnsmall, nlr, NNlr
M_small_nisol_NNsmall_nlr_NNlr = M_small_nisol_NNsmall_nlr.submask(np.log10(1+lofarcat['NN_LR']) > lLR_thresh,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN good lr'.format(s=size_large, nn=separation1),
                    'NNlr',
                    edgelabel='Y',
                    color='orange',
                    qlabel='TBC?',
                    masterlist=masterlist)

# compact not isolated, nnsmall, nlr, NNnlr
M_small_nisol_NNsmall_nlr_NNnlr = M_small_nisol_NNsmall_nlr.submask(np.log10(1+lofarcat['NN_LR']) <= lLR_thresh,
                    'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr'.format(s=size_large, nn=separation1),
                    'NNnlr',
                    edgelabel='N',
                    color='red',
                    qlabel='TBC?',
                    masterlist=masterlist)


# other masks

#maskDC0 = lofarcat['DC_Maj'] == 0
maskDC0 = lofarcat['Maj'] == 0

M_S = Mask(lofarcat['S_Code'] == 'S',
            'S',
            'single',
            verbose=True)
M_M = Mask(lofarcat['S_Code'] == 'M',
            'M',
            'multiple',
            verbose=True)
M_C = Mask(lofarcat['S_Code'] == 'C',
            'C',
            'complex',
            verbose=True)

M_Ngaus = []
for i in range(1,6):
    M_Ngaus.append(Mask(lofarcat['Ng'] == i,
                        'Ng='+str(i),
                        'Ng='+str(i),
                        verbose=True           ))

M_huge = Mask(lofarcat['Maj'] > size_huge,
                'huge',
                'huge')

M_small = Mask(lofarcat['Maj'] <= size_large,
                'small',
                'small')

M_isol = Mask(lofarcat['NN_sep'] > separation1,
                'Isolated',
                'isol')

M_cluster = Mask(lofarcat['NN5_sep'] < separation1,
                'Clustered (5 sources within sep1)',
                'clustered')

M_bright = Mask(lofarcat['Total_flux'] > fluxcut,
                'bright',
                'bright')
M_nlr = Mask(np.log10(1+lofarcat['LR']) > lLR_thresh,
                'LR good',
                'lr')
M_lr = Mask(np.log10(1+lofarcat['LR']) <= lLR_thresh,
                'LR bad',
                'nlr')

    

# make a test sample for each final mask
makesample = 0
if makesample:
    for t in masterlist:
        if not t.has_children :
            print t.name
            t.make_sample(lofarcat)

# test that the final masks are indeed mutually disjoint and cover all sources
endlist = []
for t in masterlist:
    if not t.has_children:
        endlist.append(t)
if not Masks_disjoint_complete(endlist):
    print 'WARNING: children aren\'t disjoint and complete'


# make flowchart from list of masks
plot_flowchart = True
plot_verbose = False
try:
    import pygraphviz as pgv
except ImportError:
    print 'no pygraphviz; cannot make visual flowchart'
    plot_flowchart = False
if plot_flowchart:

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

    A.add_node('start', label='ALL\n{n:d}'.format(n=M_all.N), shape='parallelogram') 
    A.add_node('m_all', label='Large?'.format(n=M_all.N), shape='diamond') 

    A.add_edge('start', 'm_all', label='Y', penwidth=M_all.f*PW)
    for t in masterlist:
        
        if t.has_children:
            shape='diamond'         # intermediate point is a question
        else:
            shape='parallelogram'   # end point is a final mask
        label=t.qlabel + '\n' + str(t.n)
        if t.color:
            c = t.color
        else:
            c = 'black'
        # add node
        A.add_node(t.name, label=label, shape=shape, color=c)
        
        # add edge to parent
        if t.has_parent:
            A.add_edge(t.parent.name, t.name, label=t.edgelabel, penwidth=t.f*PW)

    if plot_verbose:
        print(A.string()) # print dot file to standard output

    # make the flowchart
    #Optional prog=['neato'|'dot'|'twopi'|'circo'|'fdp'|'nop']
    #neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps, sccmap, tred, sfdp.
    A.layout('dot') # layout with dot
    A.draw('flow_s{s:.0f}_nn{nn:.0f}.png'.format(s=size_large,nn=separation1)) # write to file




## TESTING ##
## check gaus ML


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
_ =ax.hist(np.log10(1.+lofarcat['LR'][M_small_isol_S.mask]), bins=100, normed=True, histtype='step', label=M_small_isol_S.label)
#_ =ax.hist(np.log10(1.+lofarcat['LR'][m_small_isol_nS]), bins=100, normed=True, histtype='step', label=l_small_isol_nS)
_ =ax.hist(np.log10(1.+lofarcat['LR'][M_small_nisol.mask]), bins=100, normed=True, histtype='step', label=M_small_nisol.label)
_ =ax.hist(np.log10(1.+lofarcat['LR'][M_large.mask]), bins=100, normed=True, histtype='step', label=M_large.label)
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
counts, xedges, yedges, im =ax.hist2d(np.log10(lofarcat['Maj'][M_S.mask]), np.log10(lofarcat['Total_flux'][m_S]), bins=100, label='')
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
counts, xedges, yedges, im =ax.hist2d(np.log10(lofarcat['Maj'][~M_S.mask]), np.log10(lofarcat['Total_flux'][~m_S]), bins=100, label='')
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
counts, xedges, yedges, im =ax.hist2d((lofarcat['Maj'][M_lr.mask]), (lofarcat['NN_sep'][M_lr.mask]), bins=200, range=((0,50),(0,200)), label='')
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
counts, xedges, yedges, im =ax.hist2d((lofarcat['Maj'][~M_lr.mask]), (lofarcat['NN_sep'][~M_lr.mask]), bins=200, range=((0,50),(0,200)), label='')
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


# plot size distribution
f, ax = pp.paper_single_ax()
ax.hist(lofarcat['Maj'][~maskDC0], range=(0,80), bins=100, histtype='step', label='All')
ax.hist(lofarcat['Maj'][~maskDC0&M_S.mask], range=(0,80), bins=100, histtype='step', label='S')
ax.hist(lofarcat['Maj'][~maskDC0&M_M.mask], range=(0,80), bins=100, histtype='step', label='M')
ax.hist(lofarcat['Maj'][~maskDC0&M_C.mask], range=(0,80), bins=100, histtype='step', label='C')
ax.set_xlabel('Major Axis [arcsec]')
ax.set_ylabel('N')
ax.legend()
plt.savefig('size_dist_classes')


# plot nearest neighbour distribution
f,ax = pp.paper_single_ax()
ax.hist(lofarcat['NN_sep'], bins=100, histtype='step', label='All')
ax.hist(lofarcat['NN_sep'][M_S.mask], bins=100, histtype='step', label='S')
ax.set_xlabel('Nearest source [arcsec]')
ax.set_ylabel('N')
ax.legend()
plt.savefig('NNdist_dist')


# 2D histogram : size-nearest neighbour distance
# for 'S' sources
f,ax = pp.paper_single_ax()
X =  lofarcat['NN_sep'][~maskDC0&M_S.mask]
Y = lofarcat['Maj'][~maskDC0&M_S.mask]
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


# and 'M' sources
f,ax = pp.paper_single_ax()
X =  lofarcat['NN_sep'][~maskDC0&M_M.mask]
Y = lofarcat['Maj'][~maskDC0&M_M.mask]
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



