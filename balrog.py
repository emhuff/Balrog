#!/usr/bin/env python

import galsim
import numpy as np
import pyfits
import math
import sys
from optparse import OptionParser
import sextractor_engine 

def defineParameters():
    ngal = 50
    parameters_list = []
    for i in np.arange(ngal):
        parameters = {}
        parameters['Sersic index'] = 4.0
        parameters['half light radius'] = 2.5
        parameters['flux'] = np.random.lognormal(mean=5.,sigma=1.)
        parameters['g1'] = 0.2*np.random.randn()
        parameters['g2'] = 0.2*np.random.randn()
        parameters_list.append(parameters)
    return parameters_list

    
def defineCalibration():
    '''
    This sets the number of counts per unit flux. Used to set the
    noise of the added galaxies.
    '''
    calibParams = {}
    calibParams['gain'] = 10000.0
    return calibParams

def getBigImage(file='example.fits',subRegion= (None,None,None,None)):
    '''
    Takes a filename for the large image to be simulated. Reads in a
    smaller piece of that image defined by subRegion (which is a
    galsim._galsim.bounds object), and returns the subRegion as a GalSim
    image object.
    '''
    bigImage = galsim.fits.read(file)
    pixelScale = 0.1
    bigImage.setScale(pixelScale)
    bigCenter = bigImage.bounds.center()
    subBounds = bigImage.bounds
    if subRegion[0]  > 0:
        subBounds.xmin = subRegion[0]
    if subRegion[1]  > 0:
        subBounds.xmax = subRegion[1]
    if subRegion[2]  > 0:
        subBounds.ymin = subRegion[2]
    if subRegion[3]  > 0:
        subBounds.ymax = subRegion[3]
    bigImage = bigImage[subBounds]
    subCenter = bigImage.bounds.center()
    centOffset = subCenter - bigCenter
    offset = galsim.PositionD(centOffset.x,centOffset.y)
    return bigImage, offset

if __name__ == "__main__":
    '''
    Takes an input image .fits file name, and output file name, reads in the zeroth
    extension, adds Sersic galaxy images, and then writes the output
    to another .fits file.
    USAGE:
       SimImage inputFileName  outputFileName

    TODO:
       SExtractor interface.
    '''
    parser = OptionParser()
    parser.add_option("-i", "--image", action="store",dest="ImageFile",
                      type="string",help="Image file to be read",default="example.fits")
    parser.add_option("-o", "--output",action="store",type="string",
                      dest="OutputFile",default="example_output.fits",
                      help="file to put the output simulated image in.")
    parser.add_option("-w", "--weightin",action="store",type="string",
                      dest="WeightMapIn",default="example_weight_in.fits",
                      help="file to read the weight map from.")
    parser.add_option("-v", "--weightout",action="store",type="string",
                      dest="WeightMapOut",default="example_weight_out.fits",
                      help="file to write the extracted weight map to.")
    parser.add_option("-t", "--catalogout",action="store",type="string",
                      dest="CatalogOutFile",default="example_catalog.fits",
                      help="file to write the SExtractor catalog to.")
    parser.add_option("-p", "--psfmodel",action="store",type="string",
                      dest="PSFExFile",default="example.psfcat.psf",
                      help="File containing PSFEx psf model to use.")
    parser.add_option("--xmin",action="store",type="int",default="-1",
                      help="minimum column of extracted subImage, unit-indexed",dest="xmin")
    parser.add_option("--xmax",action="store",type="int",default="-1",
                      help="maximum column of extracted subImage, unit-indexed",dest="xmax")
    parser.add_option("--ymin",action="store",type="int",default="-1",
                      help="minimum row of extracted subImage, unit-indexed",dest="ymin")
    parser.add_option("--ymax",action="store",type="int",default="-1",
                      help="minimum row of extracted subImage, unit-indexed",dest="ymax")

    (opts,args ) = parser.parse_args()

    '''
    Work out whether the user wants to extract a piece of the original
    image to do the simulations on.
    '''
    subRegion = (opts.xmin,opts.xmax,opts.ymin,opts.ymax)

    rng = galsim.UniformDeviate()
    parameters_list = defineParameters()
    calib = defineCalibration()
    bigImage, offset = getBigImage(opts.ImageFile,subRegion=subRegion)
    psfmodel = galsim.des.DES_PSFEx(opts.PSFExFile)
    center = bigImage.bounds.center()
    

    for parameters in parameters_list:
        sersicObj = galsim.Sersic(n=parameters['Sersic index'],half_light_radius=
                                  parameters['half light radius'],flux = parameters['flux'])
        sersicObj.applyShear(g1=parameters['g1'],g2=parameters['g2'])
        x = np.random.random_sample()*bigImage.array.shape[0]
        y = np.random.random_sample()*bigImage.array.shape[1]

        ix = int(np.floor(x+.05))
        iy = int(np.floor(y+.05))
        dx = x-ix
        dy = y-iy
        pos = galsim.PositionD(x,y) * bigImage.getScale()
        sersicObj.applyShift(dx*bigImage.getScale(),dy*bigImage.getScale())
        # Convolve with the pixel.
        pix = galsim.Pixel(bigImage.getScale())
        psf = psfmodel.getPSF(pos-offset,bigImage.getScale())
        finalPSF = galsim.Convolve([pix,psf])
        sersicObjConv = galsim.Convolve([finalPSF,sersicObj])
        smallImage = sersicObjConv.draw(dx=bigImage.getScale())
        smallImage.addNoise(galsim.CCDNoise(gain=calib['gain'],read_noise=0))
        smallImage.setCenter(ix,iy)
        bounds = smallImage.bounds & bigImage.bounds
        
        bigImage[bounds] += smallImage[bounds]

    
    bigImage.write(opts.OutputFile)
    # Next, write the subImage weightmap.
    subWeight, Wcent = getBigImage(opts.WeightMapIn,subRegion=subRegion)
    subWeight.write(opts.WeightMapOut)

    eng = sextractor_engine.SextractorEngine(IMAGE=opts.OutputFile,
                                             WEIGHT_IMAGE=opts.WeightMapOut,
                                             CHECKIMAGE_TYPE='BACKGROUND',
                                             CATALOG_NAME=opts.CatalogOutFile)
    eng.run()
