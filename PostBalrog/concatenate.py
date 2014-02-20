#!/usr/bin/env python

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from esutil import htm



if __name__ == "__main__":

    infile = None
    supfile = None
    start = 0
    end = 299
    outfile = 'catalogs/rxj_i.44_index=0-299.fits'
    outsupfile = outfile.replace('.fits', '_sup.fits')
    topdir = '../../output/magnification_sextractor/rxj/i/'


    h = htm.HTM()
    matchrad = 1.0 / 3600.0


    if infile==None:
        started = False
        nums = np.empty(0)
        areas = np.empty(0)
        cindex = np.empty(0)
    else:
        started = True
        hdu = pyfits.open(infile)[1]
        data = hdu.data
        header = hdu.header
        j = 0
        newkeys = {}
        while True:
            key = 'V%i' %j
            if key not in header.keys():
                break
            newkeys[key] = header[key]
            j += 1

        sup = pyfits.open(supfile)[1].data
        nums = sup['NSIM']
        areas = sup['AREA']
        cindex = sup['CINDEX']


    #ind = subprocess.check_output( ['ls',topdir] ).strip().split()
    ind = np.arange(start,end)
    for index in ind:

        dir = os.path.join( topdir, str(index) )
        mcat = os.path.join( dir, 'measured.fits' )
        ncat = os.path.join( dir, 'measured.nosim.fits')

        if not os.path.exists(mcat):
            continue
        if not os.path.exists(ncat):
            continue

        try:
            dhdu = pyfits.open(mcat)[2]
            d = dhdu.data
            dh = dhdu.header
        except:
            print 'error openning measured %i' %index
            continue
        try:
            nhdus = pyfits.open(ncat)
            n = nhdus[2].data
        except:
            print 'error opening nosim %i' %index

        if len(n)!=0:
            mn, md, radius = h.match( n['ALPHAWIN_J2000'], n['DELTAWIN_J2000'], d['ALPHAWIN_J2000'], d['DELTAWIN_J2000'], matchrad, maxmatch=1 )
            d = np.delete(d, md)

        msize = len( d )
        if not started:
            data = d
            j = 0
            newkeys = {}
            while True:
                key = 'V%i' %j
                if key not in dh.keys():
                    break
                newkeys[key] = dh[key]
                j += 1
            started = True
        else:
            data = np.concatenate( (data,d) )

        tcat = os.path.join( dir, 'truth.fits' )
        if not os.path.exists(tcat):
            continue
        hdus = pyfits.open(tcat)
        num = hdus[1].header['NSIM']
        dx = hdus[1].header['XEND'] - hdus[1].header['XSTART'] + 1
        dy = hdus[1].header['YEND'] - hdus[1].header['YSTART'] + 1
        area = dx*dy

        i = int(index)
        n = np.array( [num]*msize )
        a = np.array( [area]*msize )
        ci = np.array( [i]*msize )
        nums = np.append(nums, n)
        areas = np.append(areas, a)
        cindex = np.append(cindex, ci)
        

    tbhdu = pyfits.new_table(data)
    for newkey in newkeys:
        tbhdu.header[newkey] = newkeys[newkey]
    phdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([phdu,tbhdu])
    if os.path.exists(outfile):
        subprocess.call(['rm',outfile])
    hdulist.writeto(outfile)

    col1 = pyfits.Column(name='NSIM', format='J', array=nums)
    col2 = pyfits.Column(name='AREA', format='E', array=areas)
    col3 = pyfits.Column(name='CINDEX', format='J', array=cindex)
    cols = pyfits.ColDefs([col1,col2,col3])
    tbhdu = pyfits.new_table(cols)
    phdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([phdu,tbhdu])
    if os.path.exists(outsupfile):
        subprocess.call(['rm',outsupfile])
    hdulist.writeto(outsupfile)

