#sync all coadd and SE images / psfs / scamp headers needed for a given
#tile and band
#requires desdb, and make sure desdb/desdb/bin is in your $PATH
import multi_epoch_tools.file_tools as file_tools
import sys
import subprocess
import numpy as np
import argparse
import pickle
usage="""
python sync_tile.py <tilename> <band> <mode>
mode should be one of 'all','coadd' or 'red'"""
"""
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
"""
def main(args):

    if args.info_dict is not None:
        try:
            coadd_dict=pickle.load(open(args.info_dict,'rb'))
        except IOError:
            coadd_dict={}
        
    for band in args.bands:
        print coadd_dict
        coadd=file_tools.setup_tile(args.tilename,band=args.bands[0],sync=False)

        if args.info_dict is not None:
            try:
                coadd_dict[args.tilename][band]={}
            except KeyError:
                coadd_dict[args.tilename]={}
                coadd_dict[args.tilename][band]={}
            for key,val in coadd.iteritems():
                coadd_dict[args.tilename][band][key]=val
            coadd_dict[args.tilename][band]['srclist']=coadd.srclist

        if args.mode in ['all','coadd'] and band==args.bands[0]:
            file_tools.sync_coadd(coadd['coadd_run'])

        if args.mode in ['all','red']:
            #Now SE images
            #Print out the runs needed to get all the single exposures for this tile
            runs=np.array([s['run'] for s in coadd.srclist])
            expnames=np.array([s['expname'] for s in coadd.srclist])
            exp_unique,inds_unique=np.unique(expnames,return_index=True)
            print exp_unique,inds_unique
            for exp,run in zip(exp_unique,runs[inds_unique]):
                subprocess.call(["des-sync-red",run,exp])
    
    if args.info_dict is not None:
        with open(args.info_dict,'wb') as f:
            pickle.dump(coadd_dict, f)


parser=argparse.ArgumentParser()
parser.add_argument('tilename',type=str,help='DES tile name')
parser.add_argument('mode',type=str,help='coadd, red, all or just info')
parser.add_argument('bands',type=str,nargs='*')
parser.add_argument('--info_dict',type=str,default=None)
args=parser.parse_args()

if __name__=="__main__":
    assert args.mode in ['coadd','red','all','info']
    main(args)
