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

    outdir = os.path.join(opts.outdir, str(index))
    if not os.path.exists(outdir):
        subprocess.call(['mkdir', outdir])

    imageout = os.path.join(outdir, 'tempimage.fits')
    weightout = os.path.join(outdir, 'tempweight.fits')
    psfout = os.path.join(outdir, 'temppsf.psf')
    catalogtruth = os.path.join(outdir, 'truth.fits')
    catalogmeasured = os.path.join(outdir, 'measured.fits')
    assocfile = os.path.join(outdir, 'assoc.txt')

    thisdir = os.path.dirname( os.path.realpath(__file__) )
    balrog = os.path.join(thisdir,'balrog.py')

    bargs = [balrog,
            '--imagein', opts.tempimage,
            '--imageout', imageout,
            '--weightin', opts.tempweight,
            '--weightout', weightout,
            '--psfin', opts.temppsf,
            '--psfout', psfout,
            '--catalogtruth', catalogtruth,
            '--catalogmeasured', catalogmeasured,
            '--ngal', str(opts.ngal),
            '--band', opts.band,
            '--zeropoint', str(opts.zeropoint),
            '--sexpath', opts.sexpath]

    if opts.noassoc:
        bargs.append('--noassoc')
    else:
        bargs.append('--assoc')
        bargs.append(assocfile)

    if opts.clean:
        bargs.append('--clean')

    if opts.noempty:
        bargs.append('--noempty')

    subprocess.call( bargs )


def RunBalrog(opts):
    index = opts.index
    args = []
    for i in range(opts.nsim):
        arg = (index, opts)
        index += 1
        args.append(arg)
    pool = Pool(opts.nthreads)
    result = pool.map(Run, args)


def SubsampleImages(opts):
    balrog.WritePsf(opts, opts.inpsf, opts.temppsf)
    image, head = balrog.ReadImage(opts,opts.inimage)
    balrog.WriteImage(image, opts.tempimage, header=head)
    image, head = balrog.ReadImage(opts,opts.inweight)
    balrog.WriteImage(image, opts.tempweight, header=head)


def ParseArgs(opts):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    if opts.dir==None:
        updir = '/'.join(thisdir.split('/')[:-1])
        opts.dir = os.path.join(updir, 'output')
        if not os.path.exists(opts.dir):
            subprocess.call(['mkdir', opts.dir])
        opts.dir = os.path.join(opts.dir, 'magnification_sextractor')
        if not os.path.exists(opts.dir):
            subprocess.call(['mkdir', opts.dir])
        opts.dir = os.path.join(opts.dir, '%s'%(opts.clus))
        if not os.path.exists(opts.dir):
            subprocess.call(['mkdir', opts.dir])

    opts.outdir = os.path.join(opts.dir, opts.band)
    if not os.path.exists(opts.outdir):
        subprocess.call(['mkdir', opts.outdir])

    imstr = '%s_%s.%s' %(opts.clus,opts.band,opts.code)
    opts.inimage = os.path.join(opts.datadir, opts.clus, 'images', '%s.fits'%(imstr))
    opts.inweight = os.path.join(opts.datadir, opts.clus, 'weights', '%s.fits'%(imstr))
    opts.inpsf = os.path.join(opts.datadir, opts.clus, 'psfmodels', '%s_det_%s.%s.faint.deg%i'%(imstr,opts.detband,opts.detcode,opts.psfdeg), '%s.psfcat.psf'%(imstr))

    num = opts.nsim - 1
    index = '%i_to_%i' %(opts.index,num)
    opts.temppsf = os.path.join(opts.outdir, 'psf_%s.psf'%(index))
    opts.tempimage = os.path.join(opts.outdir, 'image_%s.fits'%(index))
    opts.tempweight = os.path.join(opts.outdir, 'weight_%s.fits'%(index))
    opts.zeropoint = pyfits.open(opts.inimage)[0].header['AVG_ZP']

    opts.nosim_imageout = opts.tempimage
    opts.nosim_weightout = opts.tempweight
    opts.nosim_catalogmeasured = opts.nosim_imageout.replace('image_','catalog_')
    opts.psfout = opts.temppsf
    opts.assoc = None

    configdir = os.path.join(thisdir, 'astro_config')
    if opts.sexconfig==None:
        opts.sexconfig = os.path.join(configdir, 'sex.config')
    if opts.sexemptyparam==None:
        opts.sexemptyparam = os.path.join(configdir, 'sex.param')
    if opts.sexnnw==None:
        opts.sexnnw = os.path.join(configdir, 'sex.nnw')
    if opts.sexconv==None:
        opts.sexconv = os.path.join(configdir, 'sex.conv')

    return opts


def GetOpts():
    parser = argparse.ArgumentParser()
    parser.add_argument( "-dd", "--datadir", help="input image area", type=str, default='/n/des/suchyta.1/des/SV/SV_clusters_project_home/coadd_products/')
    parser.add_argument( "-c", "--clus", help="which cluster image to use", type=str, default='rxj')
    parser.add_argument( "-deg", "--psfdeg", help="psf degree to use", type=int, default=8)
    parser.add_argument( "-mb", "--band", help="measurement filter band to use", type=str, default='i', choices=['g','r','i','z'])
    parser.add_argument( "-mc", "--code", help="measurement filter band code to use", type=str, default='44')
    parser.add_argument( "-db", "--detband", help="detection filter band to use", type=str, default='rriizz')
    parser.add_argument( "-dc", "--detcode", help="detection filter band code to use", type=str, default='43.42.44.40.44.40')
    
    parser.add_argument( "-ns", "--nsim", help="number of simulation realizations", type=int, default=100)
    parser.add_argument( "-ng", "--ngal", help="number of galxies per simulation realization", type=int, default=40)
    parser.add_argument( "-xmin", "--xmin", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=8001)
    parser.add_argument( "-xmax", "--xmax", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=20000)
    parser.add_argument( "-ymin", "--ymin", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=8001)
    parser.add_argument( "-ymax", "--ymax", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=20000)
    parser.add_argument( "-ps", "--pixscale", help="pixel scale in arcsec/pixel", type=float, default=0.263)

    parser.add_argument( "-spp", "--sexpath", help='Path for sextractor binary', type=str, default='/n/des/suchyta.1/des/SV/SV_clusters_project_home/pipeline/astro_software/sextractor-trunk-multithreaded_install/bin/sex')
    parser.add_argument( "-a", "--noassoc", help="don't run sextractor in association mode", action='store_true')
    parser.add_argument( "-ne", "--noempty", help="Don't do sextractor run that doesn't have simulated galaxies", action="store_true")

    parser.add_argument( "-sc", "--sexconfig", help='Sextractor config file', type=str, default=None)
    parser.add_argument( "-sn", "--sexnnw", help='Sextractor neural network file', type=str, default=None)
    parser.add_argument( "-sv", "--sexconv", help='Sextractor filter convolution file', type=str, default=None)
    parser.add_argument( "-sep", "--sexemptyparam", help='Sextractor param file used when image without simulated galaxies is extracted. This is basically irrelevant, since all the run is inteded to do check for things in image before simulating. Not doing model fitting is faster, so the default is to use one that does not do model fitting.', type=str, default=None)


    parser.add_argument( "-d", "--dir", help="label for what directory to write output to", type=str, default=None)
    parser.add_argument( "-i", "--index", help="index to start at for labeling output", type=int, default=0)
    parser.add_argument( "-cl", "--clean", help="Clean images", action='store_true')

    parser.add_argument( "-nt", "--nthreads", help="number of threads to use", type=int, default=1)

    args = parser.parse_args()
    args = ParseArgs(args)
    return args


if __name__ == "__main__":

    opts = GetOpts()

    SubsampleImages(opts)
    RunBalrog(opts)
    balrog.RunSextractor(opts, [], nosim=True)

