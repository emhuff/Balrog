#!/usr/bin/env python

import galsim
import numpy as np
import pyfits
import math
import sys

def defineParameters():
    ngal = 20
    parameters_list = []
    for i in np.arange(ngal):
        parameters = {}
        parameters['Sersic index'] = 4.0
        parameters['half light radius'] = 2.5
        parameters['flux'] = np.random.lognormal(mean=5.,sigma=1.)
        parameters['xpos'] = np.random.random_sample()*1000
        parameters['ypos'] = np.random.random_sample()*1000
        parameters['g1'] = 0.1*np.random.randn()
        parameters['g2'] = 0.1*np.random.randn()
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

def getBigImage(file='example.fits',subRegion= None):
    '''
    Takes a filename for the large image to be simulated. Reads in a
    smaller piece of that image defined (somehow) by subRegion, and
    returns the subRegion as a GalSim imageF object.
    '''
    bigImage = galsim.fits.read(file)
    pixelScale = 0.1
    bigImage.setScale(pixelScale)
    return bigImage

if __name__ == "__main__":
    '''
    Takes an input image .fits file name, and output file name, reads in the zeroth
    extension, adds Sersic galaxy images, and then writes the output
    to another .fits file.
    TODO:
       Implement subRegion definition
       Implement psf model and convolution.
    '''
    
    if len(sys.argv) > 1:
        ImageFile = sys.argv[1]
    else:
        ImageFile = 'example.fits'
    if len(sys.argv) > 2:
        OutputFile = sys.argv[2]
    else:
        OutputFile = 'example_with_simimages.fits'
        
    rng = galsim.UniformDeviate()
    parameters_list = defineParameters()
    calib = defineCalibration()
    bigImage = getBigImage(ImageFile)
    center = bigImage.bounds.center()
    for parameters in parameters_list:
        sersicObj = galsim.Sersic(n=parameters['Sersic index'],half_light_radius=
                                  parameters['half light radius'],flux = parameters['flux'])
        sersicObj.applyShear(g1=parameters['g1'],g2=parameters['g2'])
        x = parameters['xpos']
        y = parameters['ypos']
        ix = int(np.floor(x+.05))
        iy = int(np.floor(y+.05))
        dx = x-ix
        dy = y-iy
        pos = galsim.PositionD(x,y) * bigImage.getScale()
        smallImage = sersicObj.draw(dx=bigImage.getScale())
        smallImage.addNoise(galsim.CCDNoise(gain=calib['gain'],read_noise=0))
        smallImage.setCenter(ix,iy)
        bounds = smallImage.bounds & bigImage.bounds

        bigImage[bounds] += smallImage[bounds]
    
    bigImage.write(OutputFile)

    
