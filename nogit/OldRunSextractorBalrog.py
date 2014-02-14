#!/usr/bin/env python

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
    image = args[2]
    weight = args[3]
    psf = args[4]

    if opts.where=='constant':
        xmin = opts.wherex
        xmax = xmin + opts.size - 1
        ymin = opts.wherey
        ymax = ymin + opts.size - 1
    elif opts.where=='random':
        np.random.seed()
        subxmax = opts.xmax - opts.xmin + 1
        xmin = np.random.random_integers(1, subxmax-opts.size+ 1)
        xmax = xmin + opts.size - 1
        subymax = opts.ymax - opts.ymin + 1
        ymin = np.random.random_integers(1, subymax-opts.size+ 1)
        ymax = ymin + opts.size - 1
    
    odir = os.path.join(opts.dir, opts.band)
    if not os.path.exists(odir):
        subprocess.call(['mkdir', odir])
    outdir = os.path.join(odir, str(index))
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
            '--imagein', image,
            '--imageout', imageout,
            '--weightin', weight,
            '--weightout', weightout,
            '--psfin', psf,
            '--psfout', psfout,
            '--catalogtruth', catalogtruth,
            '--catalogmeasured', catalogmeasured,
            '--xmin', str(xmin),
            '--xmax', str(xmax),
            '--ymin', str(ymin),
            '--ymax', str(ymax),
            '--ngal', str(opts.ngal),
            '--band', opts.band,
            '--sexpath', opts.sexpath]

    if opts.noassoc:
        bargs.append('--noassoc')
    else:
        bargs.append('--assoc')
        bargs.append(assocfile)

    if opts.nomodel:
        bargs.append('--nomodel')

    if opts.clean:
        bargs.append('--clean')

    subprocess.call( bargs )


def Subsample(inimage, opts, label):
    if label=='psf':
        outimage = 'temppsf.psf'
        balrog.WritePsf(opts, inimage, outimage)
    else:
        image, ihead = balrog.ReadImage(opts,inimage)
        outimage = 'temp%s.fits' %label
        balrog.WriteImage(image, outimage, header=ihead)

    return outimage


def SubsampleImages(opts):
    imstr = '%s_%s.%s' %(opts.clus,opts.band,opts.code)
    inimage = os.path.join(opts.datadir, opts.clus, 'images', '%s.fits'%(imstr))
    inweight = os.path.join(opts.datadir, opts.clus, 'weights', '%s.fits'%(imstr))
    inpsf = os.path.join(opts.datadir, opts.clus, 'psfmodels', '%s_det_%s.%s.faint.deg%i'%(imstr,opts.detband,opts.detcode,opts.psfdeg), '%s.psfcat.psf'%(imstr))

    subimage = Subsample(inimage, opts, 'image')
    subweight = Subsample(inweight, opts, 'weight')
    subpsf = Subsample(inpsf, opts, 'psf')

    return subimage,subweight,subpsf


def RunBalrog(opts, subimage, subweight, subpsf):
    index = opts.index
    args = []
    for i in range(opts.nsim):
            arg = (index, opts, subimage, subweight, subpsf)
            index += 1
            args.append(arg)
    pool = Pool(opts.nthreads)
    result = pool.map(Run, args)
    #Run(args[0])


def GetOpts():
    parser = argparse.ArgumentParser()
    parser.add_argument( "-c", "--clus", help="which cluster image to use", type=str, default='rxj')
    parser.add_argument( "-deg", "--psfdeg", help="psf degree to use", type=int, default=8)
    parser.add_argument( "-ps", "--pixscale", help="pixel scale in arcsec/pixel", type=float, default=0.263)
    
    parser.add_argument( "-ns", "--nsim", help="number of simulation realizations", type=int, default=100)
    parser.add_argument( "-s", "--size", help="simulation image is size x size pixels", type=int, default=2000)
    parser.add_argument( "-ng", "--ngal", help="number of galxies per simulation realization", type=int, default=40)
    parser.add_argument( "-w", "--where", help="how to choose image locations", type=str, default='random', choices=['random','constant'])
    parser.add_argument( "-wx", "--wherex", help="constant image location x (bottom left corner), in the smaller tempimage", type=int, default=1)
    parser.add_argument( "-wy", "--wherey", help="constant image location y (bottom left corner), in the smaller tempimage", type=int, default=1)

    parser.add_argument( "-a", "--noassoc", help="don't run sextractor in association mode", action='store_true')
    parser.add_argument( "-nm", "--nomodel", help="don't fit a model", action='store_true')
    parser.add_argument( "-cl", "--clean", help="Clean images", action='store_true')

    parser.add_argument( "-mb", "--band", help="measurement filter band to use", type=str, default='i', choices=['g','r','i','z'])
    parser.add_argument( "-mc", "--code", help="measurement filter band code to use", type=str, default='44')
    parser.add_argument( "-db", "--detband", help="detection filter band to use", type=str, default='rriizz')
    parser.add_argument( "-dc", "--detcode", help="detection filter band code to use", type=str, default='43.42.44.40.44.40')

    parser.add_argument( "-xmin", "--xmin", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=8000)
    parser.add_argument( "-xmax", "--xmax", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=20000)
    parser.add_argument( "-ymin", "--ymin", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=8000)
    parser.add_argument( "-ymax", "--ymax", help="only use a certain piece of the image-->[ymin:ymax, xmin:xmax]. Use to need less memory for opening image", type=int, default=20000)

    parser.add_argument( "-d", "--dir", help="label for what directory to write to", type=str, default=None)
    parser.add_argument( "-dd", "--datadir", help="input image area", type=str, default='/n/des/suchyta.1/des/SV/SV_clusters_project_home/coadd_products/')
    parser.add_argument( "-i", "--index", help="index to start at for labeling output", type=int, default=0)
    parser.add_argument( "-nt", "--nthreads", help="number of threads to use", type=int, default=1)
    parser.add_argument( "-sp", "--sexpath", help='Path for sextractor binary', type=str, default='sex')
    args = parser.parse_args()

    thisdir = os.path.dirname( os.path.realpath(__file__) )
    if args.dir==None:
        updir = '/'.join(thisdir.split('/')[:-1])
        args.dir = os.path.join(updir, 'output')
        if not os.path.exists(args.dir):
            subprocess.call(['mkdir', args.dir])
        args.dir = os.path.join(args.dir, 'magnification_sextractor')
        if not os.path.exists(args.dir):
            subprocess.call(['mkdir', args.dir])
        args.dir = os.path.join(args.dir, '%s'%(args.clus))
        if not os.path.exists(args.dir):
            subprocess.call(['mkdir', args.dir])

    return args



if __name__ == "__main__":
    opts = GetOpts()
    print 'Subsampling original image' 
    subimage,subweight,subpsf = SubsampleImages(opts)
    print 'Done subsampling'
    print 'Running Balrog'
    RunBalrog(opts, subimage, subweight, subpsf)
    print 'Done running Balrog' 
