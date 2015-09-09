#sync all coadd and SE images / psfs / scamp headers needed for a given
#tile and band
#requires desdb, and make sure desdb/desdb/bin is in your $PATH
import multi_epoch_tools.file_tools as file_tools
import sys
import subprocess
import numpy as np

usage="""
python sync_tile.py <tilename> <band> <mode>
mode should be one of 'all','coadd' or 'red'"""

try:
    tilename,band,mode = sys.argv[1],sys.argv[2],sys.argv[3]
except IndexError:
    print usage
    exit()

if mode not in ['all','coadd','red']:
    raise KeyError("third argument should be mode, one of 'all','coadd','red'")

if mode=='all' or mode=='coadd':
    coadd=file_tools.setup_tile(tilename, band=band, sync=True)
else:
    coadd=file_tools.setup_tile(tilename, band=band, sync=False)

if mode=='all' or mode=='red':
    #Now SE images
    #Print out the runs needed to get all the single exposures for this tile
    runs=np.array([s['run'] for s in coadd.srclist])
    expnames=np.array([s['expname'] for s in coadd.srclist])
    exp_unique,inds_unique=np.unique(expnames,return_index=True)
    print exp_unique,inds_unique
    for exp,run in zip(exp_unique,runs[inds_unique]):
        subprocess.call(["des-sync-red",run,exp])
