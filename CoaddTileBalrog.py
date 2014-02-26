#!/usr/bin/env python

import time
import sys
import os
import subprocess
from multiprocessing import Pool
import argparse
import numpy as np
import pyfits
import balrog


def Run(args):
    index = args[0]
    opts = args[1]

    tiledir = os.path.join(opts.datadir, opts.tile)
    filelabel = '%s_%s' %(opts.tile, opts.band)
    outdir = os.path.join(tiledir, 'Balrog', opts.band, str(index))
    imagein = os.path.join(tiledir, 'image', '%s.fits'%(filelabel))
    psfin = os.path.join(tiledir, 'psf', '%s_psfcat.psf'%(filelabel))

    bargs = ['./balrog.py',
             '--outdir', outdir,
             '--imagein', imagein,
             '--psfin', psfin,
             '--ngal', opts.ngal]

    #print ' '.join(bargs); sys.stdout.flush()
    subprocess.call( bargs )


def RunBalrog(opts):
    nsim = opts.ntot / opts.ngal
    if opts.ntot % opts.ngal != 0:
        nsim += 1

    index = opts.index
    args = []
    for i in range(nsim):
        arg = (index, opts)
        args.append(arg)
        index += 1
    '''
    pool = Pool(opts.nthreads)
    result = pool.map(Run, args)
    '''
    Run(args[0])


def GetOpts():
    parser = argparse.ArgumentParser()
    parser.add_argument( "-t", "--tile", help="Tilename", type=str, required=True)
    parser.add_argument( "-b", "--band", help="Filter band", type=str, required=True)
    parser.add_argument( "-ng", "--ngal", help="Number of simulated galaxies per simulation", type=int, default=1e3)
    parser.add_argument( "-nt", "--ntot", help="Total number of simulated galaxies. This will determine how many simulations to run", type=int, default=5e5)
    parser.add_argument( "-i", "--index", help="Where to start the index for labelling the simulations (esssentially exists so you can 'continue where you left off')", type=int, default=0)
    parser.add_argument( "-nc", "--nthreads", help="Number of threads to run with", type=int, default=1)
    parser.add_argument( "-d", "--datadir", help="Directory where the tiles live", type=str, default='/n/des/suchyta.1/balrog_sva1/')
    args = parser.parse_args()
    
    return args


if __name__ == "__main__":

    opts = GetOpts()
    RunBalrog(opts)

