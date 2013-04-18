#!/usr/bin/env python

import subprocess
import os
import numpy

from run_sextractor import *

def run_sex(args):
    """
    arg_str = ''
    for i in range( len(args) ):
        for j in range( len(args[i]) ):
            arg_str = arg_str + '%s ' %args[i][j]
    arg_str = arg_str.strip()

    #os.system( './run_sextractor.py %s' %arg_str )
    """

    this_dir = '/'.join( os.path.realpath(__file__).split('/')[:-1] )
    exc = os.path.join( this_dir, 'run_sextractor.py' )
    args.insert( 0, exc )
    subprocess.call(args)
    return 'Done'


if __name__=='__main__':

    cmd_args= ['-i', '/n/cima/cima02/des/SV/suchyta/cluster_coadds/sptw2/images/sptw2_r.13.fits',
               '-w', '/n/cima/cima02/des/SV/suchyta/cluster_coadds/sptw2/weights/sptw2_r.13.fits']
    status = run_sex( cmd_args )
    print status
