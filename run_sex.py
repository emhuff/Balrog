#!/usr/bin/env python

import run_sextractor
import subprocess
import os
import numpy

def run_sex(args):
    """
    arg_str = ''
    for i in range( len(args) ):
        for j in range( len(args[i]) ):
            arg_str = arg_str + '%s ' %args[i][j]
    arg_str = arg_str.strip()

    #os.system( './run_sextractor.py %s' %arg_str )
    """
    args.insert( 0, './run_sextractor.py' )
    subprocess.call(args)
    return 'Done'


if __name__=='__main__':

    cmd_args= ['-i', '/n/cima/cima02/des/SV/suchyta/cluster_coadds/sptw2/images/sptw2_r.13.fits',
               '-w', '/n/cima/cima02/des/SV/suchyta/cluster_coadds/sptw2/weights/sptw2_r.13.fits']
    status = run_sex( cmd_args )
    print status
