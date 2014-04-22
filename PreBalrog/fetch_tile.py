#!/usr/bin/env python

import utils
import argparse
import os
import sys
import subprocess


def checkdir(dir):
    if not os.path.exists(dir):
        mkdir(dir)


def mkdir(dir):
    subprocess.call(['mkdir', dir])


def funpack(ifile):
    cimage = file
    ofile = ifile.replace('.fits.fz', '.fits')
    if os.path.exists(ofile):
        subprocess.call( ['rm', ofile] )
    #print ifile
    #print ofile
    subprocess.call( ['funpack', '-O', ofile, ifile] )


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument( "-b", "--band", help="Filter band", type=str, default='i')
    parser.add_argument( "-t", "--tile", help="Coadd tile name", type=str, required=True)
    parser.add_argument( "-d", "--directory", help="Directory where your output will go. Tiles will have subdirectories under this.", default=os.environ['SVA1_TILE'])
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    
    args = GetArgs()
    cur = utils.get_cursor()
    cur.execute("select path, tilename, imagename from sva1_coadd_file where tilename='%s' and band='%s'" %(args.tile, args.band))
    path, tile, name = cur.fetchone()

    filepath = os.path.join( utils.db_specs.file_archive, path)
    catpath = filepath.replace('.fits.fz', '_cat.fits')
    psfcatpath = filepath.replace('.fits.fz', '_psfcat.fits')
    psfpath = psfcatpath.replace('.fits', '.psf')
    remotepaths = [filepath,catpath,psfcatpath,psfpath]

    topdir = os.path.join(args.directory, tile)
    checkdir(topdir)

    subdirs = [] 
    outpaths = []
    sdirs = ['image','cat','psfcat','psf']
    for sub in sdirs:
        s = os.path.join( topdir,sub )
        subdirs.append(s)
    for subdir,remotepath in zip(subdirs,remotepaths):
        checkdir(subdir)
        name = remotepath.split('/')[-1]
        outpaths.append( os.path.join(subdir,name) )
    
    wgetcommon = 'wget --no-check-certificate -q'
    for remotepath, outpath in zip(remotepaths,outpaths):
        os.system('%s %s -O %s' %(wgetcommon,remotepath,outpath) )

    cimage = outpaths[0]
    funpack(cimage)
    
