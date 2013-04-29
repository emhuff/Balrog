#!/usr/bin/env python

import galsim
import numpy as np
import pyfits
import math
import sys

def defineParameters():
    ngal = 100
    parameters_list = []
    for i in np.arange(ngal):
        parameters = {}
        parameters['Sersic index'] = 1.0
        parameters['half light radius'] = 2.5
        parameters['flux'] = 1000.
        parameters['xpos'] = np.random.random_sample()*2000
        parameters['ypos'] = np.random.random_sample()*2000
        parameters_list.append(parameters)
    return parameters_list

def defineCalibration():
    pass

def getBigImage(file,subRegion= None):
    '''
    Takes a filename for the large image to be simulated. Reads in a
    smaller piece of that image defined (somehow) by subRegion, and
    returns the subRegion as a GalSim imageF object.
    '''
    
    pixelScale = 0.1
    imageSize = 2000
    bigImage = galsim.ImageF(imageSize,imageSize)
    bigImage.setScale(pixelScale)
    return bigImage

if __name__ == "__main__":
    '''
    
    '''
    ImageFile = sys.argv[1]    
    print ImageFile ## TEMP; Remove during final implementation
    rng = galsim.UniformDeviate()
    parameters_list = defineParameters()
    bigImage = getBigImage(ImageFile)
    center = bigImage.bounds.center()
    for parameters in parameters_list:
        sersicObj = galsim.Sersic(n=parameters['Sersic index'],half_light_radius=
                                  parameters['half light radius'],flux = parameters['flux'])
        sersicObj.applyShear(g1=0.2,g2=0)
        x = parameters['xpos']
        y = parameters['ypos']
        ix = int(np.floor(x+.05))
        iy = int(np.floor(y+.05))
        dx = x-ix
        dy = y-iy
        pos = galsim.PositionD(x,y) * bigImage.getScale()
        smallImage = sersicObj.draw(dx=bigImage.getScale())
        smallImage.addNoise(galsim.CCDNoise(gain=1000,read_noise=0))
        smallImage.setCenter(ix,iy)
        bounds = smallImage.bounds & bigImage.bounds
        bigImage[bounds] += smallImage[bounds]
    
    bigImage.write('test.fits')

    
