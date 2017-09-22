import glob
import os

ddir = 'test_gaus/'

fitslist = glob.glob(ddir+'*.fits')

for f in fitslist:
    fdir = f.replace('.fits','')
    filelist = glob.glob(fdir+'/*')
    print filelist
    cmd = 'montage '+' '.join(filelist)+'  -tile 5x -geometry 256x256+1+1 '+fdir+'.pdf'
    os.system(cmd)

