import galsim
import numpy as np
import pyfits
import math

if __name__ == "__main__":
    '''
    SimImage takes as argument a list of dictionaries containing
    Sersic parameters and shears, and some other parameters that
    define output and working directory locations.
    '''
    rng = galsim.UniformDeviate()
    pixel_scale = 0.1
    parameters = {}
    parameters['Sersic index'] = 1.0
    parameters['half light radius'] = 2.5
    parameters['flux'] = 10000.
    sersicObj = galsim.Sersic(n=parameters['Sersic index'],half_light_radius=parameters['half light radius'], flux = parameters['flux'])
    sersicObj.applyShear(g1=0.2,g2=0)
    image_size = 2000
    bigImage = galsim.ImageF(image_size,image_size)
    bigImage.setScale(pixel_scale)
    center = bigImage.bounds.center()
    x = rng() * (image_size-1) + 1
    y = rng() * (image_size-1) + 1
    ix = int(np.floor(x+.05))
    iy = int(np.floor(y+.05))
    dx = x-ix
    dy = y-iy
    pos = galsim.PositionD(x,y) * pixel_scale
    smallImage = sersicObj.draw(dx=pixel_scale)
    smallImage.setCenter(ix,iy)
    bounds = smallImage.bounds & bigImage.bounds
    bigImage[bounds] += smallImage[bounds]
    
    bigImage.write('test.fits')

