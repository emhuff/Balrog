#!/usr/bin/env python

import subprocess
import os
import pyfits

# region = [ (BL row, BL col),  (TR row, TR col) ]. Use index, i.e. starting at 0
def get_subimage(image, region):
    return image[ region[0][0]:(region[1][0]+1), region[0][1]:(region[1][1]+1) ]


if __name__=='__main__':

    file = '/n/cima/cima02/des/SV/suchyta/cluster_coadds/sptw2/images/sptw2_r.13.fits'
    hdus = pyfits.open(file)
    image = hdus[0].data
    hdus.close()

    subimage = get_subimage(image, [ (10000,10000), (12000,12000) ] )
    phdu = pyfits.PrimaryHDU(subimage)
    phdu.writeto('test_subimage.fits')
