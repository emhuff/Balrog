#!/usr/bin/env python

import galsim
import numpy as np
import pyfits
import pywcs
import math
import sys
from optparse import OptionParser
import sextractor_engine


def defineParameters(x=None,y=None):
    '''
This assignes parameters (some randomly generated) of an object centered at (x,y)
'''
    parameters = {}
    parameters['Sersic index'] = 1.0 # Why 1?
    parameters['half light radius'] = 0.6 # arcsec 
    parameters['flux'] = 10000*np.random.rand()+100. # Random generated flux
    parameters['g1'] = 0.2*np.random.randn() # Randomly generated gravitational shear (weak lensing?)
    parameters['g2'] = 0.2*np.random.randn() # Randomly generated gravitational shear (weak lensing?)
    parameters['x'] = x # x coordinate of center
    parameters['y'] = y # y coordinate of center
    return parameters

    
def defineCalibration():
# Counts of what?
    '''
This sets the number of counts per unit flux. Used to set the
noise of the added galaxies.
'''
    calibParams = {}
    calibParams['gain'] = 10000.0
    calibParams['pixel_scale'] = 0.263
    return calibParams

def writeFitsCatalog(catalog,fileName):
    '''
This assumes that catalog is a list of dictionaries.
It writes out a Fits catalog version of the catalog given to location
fileName.
'''
    columns = []
    for key in catalog[0]: # Cycles through all the keys in the first dictionary (all the same?)
        # Make an array from all the entries with this key name
        arr = []
        for entry in catalog: # For each dictionary in the catalog
            arr.append(entry[key]) # Make an array of all the elements of a certain key
        arr = np.array(arr)
        # Make a FITS column with the name of the key and the array of elements
        columns.append(pyfits.Column(name=key,format='E',array=arr))
    tbhdu = pyfits.new_table(pyfits.ColDefs(columns)) # Make a table of a Header Data Unit
    hdu = pyfits.PrimaryHDU() # Primary Header Data Unit
    thdulist = pyfits.HDUList([hdu,tbhdu])
    thdulist.writeto(fileName,clobber=True) #Write out the HDU list (2 long?)

def writeFitsImage(image,outFile,wcs=None):
    '''
Write a GalSim Image to extension 0 of a .fits header.
Optionally, write the wcs.
'''
    imArr = image.array
    header = wcs.to_header() # World Coordinate System. But what is it?
    hdu = pyfits.PrimaryHDU(imArr,header=header)
    hdu.writeto(outFile,clobber=True)
    
# If this is returning a subregion, then why is it named getBigImage?
def getBigImage(file='example.fits',subRegion= (None,None,None,None),calibration=None):
    '''
Takes a filename for the large image to be simulated. Reads in a
smaller piece of that image defined by subRegion (which is a
galsim._galsim.bounds object), and returns the subRegion as a GalSim
image object.
We want to preserve the wcs.
'''

    bigImage = galsim.fits.read(file)
    if calibration is not None:
        pixelScale = calibration['pixel_scale']
    else:
        pixelScale = 0.27
    bigImage.setScale(pixelScale)

    # Define the bounds of the sub image
    subBounds = bigImage.bounds
    if subRegion[0] > 0:
        subBounds.xmin = subRegion[0]
    if subRegion[1] > 0:
        subBounds.xmax = subRegion[1]
    if subRegion[2] > 0:
        subBounds.ymin = subRegion[2]
    if subRegion[3] > 0:
        subBounds.ymax = subRegion[3]

    bigImage = bigImage[subBounds]# Crop bigImage to the new bounds

    # Preserve the offset pixel coordinates, because the PSF model is with respect to the original 
    offset = galsim.PositionD(subRegion[0],subRegion[2])
    #bigImage.setOrigin(1,1)

    # Extracts image WCS
    hdulist = pyfits.open(file)
    #hdulist[0].header['CRPIX1'] -= subRegion[0]
    #hdulist[0].header['CRPIX2'] -= subRegion[2]
    wcs = pywcs.WCS(hdulist[0].header)
    hdulist.close()
    
    return bigImage, offset, wcs

if __name__ == "__main__":
    '''
Takes an input image .fits file name, and output file name, reads in the zeroth
extension, adds Sersic galaxy images, and then writes the output
to another .fits file.
USAGE:
SimImage inputFileName outputFileName

TODO:
SExtractor interface.
'''
    # Terminal command flags
    parser = OptionParser()
    parser.add_option("-i", "--image", action="store",dest="ImageFile",
                      type="string",help="Image file to be read",default="example.fits")
    parser.add_option("-o", "--output",action="store",type="string",
                      dest="OutputFile",default="example_output.fits",
                      help="File to write the output simulated image to.")
    parser.add_option("-w", "--weightin",action="store",type="string",
                      dest="WeightMapIn",default="example_weight_in.fits",
                      help="File to read the weight map from.")
    parser.add_option("-v", "--weightout",action="store",type="string",
                      dest="WeightMapOut",default="example_weight_out.fits",
                      help="File to write the extracted weight map to.")
    parser.add_option("--catalogout",action="store",type="string",
                      dest="CatalogOutFile",default="example_catalog.fits",
                      help="File to write the SExtractor catalog to.")
    parser.add_option("--catalogIn",action="store",type="string",
                      dest="CatalogInFile",default="input_catalog.fits",
                      help="File to write the input (simulated) catalog to. Must not exist.")
    parser.add_option("-p", "--psfmodel",action="store",type="string",
                      dest="PSFExFile",default="example.psfcat.psf",
                      help="File containing PSFEx psf model to use.")
    parser.add_option("--xmin",action="store",type="int",default="0",
                      help="Minimum column of extracted subImage, unit-indexed",dest="xmin")
    parser.add_option("--xmax",action="store",type="int",default="-1",
                      help="Maximum column of extracted subImage, unit-indexed",dest="xmax")
    parser.add_option("--ymin",action="store",type="int",default="0",
                      help="Minimum row of extracted subImage, unit-indexed",dest="ymin")
    parser.add_option("--ymax",action="store",type="int",default="-1",
                      help="Minimum row of extracted subImage, unit-indexed",dest="ymax")
    parser.add_option("--ngal",action="store",type="int",default="50",
                      help="Number of simulated galaxies to add.",dest="ngal")


    (opts,args ) = parser.parse_args()

    '''
Work out whether the user wants to extract a piece of the original
image to do the simulations on.
'''
    subRegion = (opts.xmin,opts.xmax,opts.ymin,opts.ymax)# Limits the area being observed

    rng = galsim.UniformDeviate() # Random number generator, never gets used?
    calib = defineCalibration()# Image calibration
    bigImage, offset, wcs = getBigImage(opts.ImageFile,subRegion=subRegion,calibration=calib)
    psfmodel = galsim.des.DES_PSFEx(opts.PSFExFile)# Uses the DES PSF model
    center = bigImage.bounds.center()
    inputCatalog =[] # Catalog of galaxies added to the image

    for i in range(opts.ngal): # For every simulated galaxy to add to the image
        x = np.random.random_sample()*bigImage.array.shape[1]+subRegion[0] # Random x position in subRegion
        y = np.random.random_sample()*bigImage.array.shape[0]+subRegion[2] # Random y position in subRegion

        # Dictionary of parameter values (this is where the galaxy is defined)
        parameters = defineParameters(x=x-subRegion[0],y=y-subRegion[2])
        inputCatalog.append(parameters) # Add this galaxy's information to the catalog of galaxies

        # Make the galaxy
        sersicObj = galsim.Sersic(n=parameters['Sersic index'],
                                  half_light_radius=parameters['half light radius'],
                                  flux = parameters['flux'],
                                  trunc=5*parameters['half light radius'])
        sersicObj.applyShear(g1=parameters['g1'],g2=parameters['g2'])
        
        ix = int(np.floor(x)) # Integer floor of x
        iy = int(np.floor(y)) # Integer floor of y
        dx = x-ix # Decimal of x
        dy = y-iy # Decimal of y
        pos = galsim.PositionD(x,y) # Double position of x and y
        # Shifts the galaxy in the image, as otherwise the galaxy would only have a random integer position
        sersicObj.applyShift(dx,dy)

        # Make the pixel model
        pix = galsim.Pixel(bigImage.getScale())

        # Build psf model. This is where things usually go wrong.
        psf = psfmodel.getPSF(pos,bigImage.getScale())
        psf.setFlux(1.)

        # Convolve the galaxy with the psf
        sersicObj = galsim.Convolve([psf,sersicObj])

        # Convolve result with the pixel model
        sersicObjConv = galsim.Convolve([pix,sersicObj])

        # Define the size of the image containing a simulated galaxy
        postageStampSize = int(max(np.ceil(6*parameters['half light radius']),25))
        smallImage = galsim.ImageD(postageStampSize,postageStampSize)# Image to contain galaxy
        smallImage = sersicObjConv.draw(dx=bigImage.getScale(),image=smallImage)# Put the galaxy into the image
        smallImage.addNoise(galsim.CCDNoise(gain=calib['gain'],read_noise=0))# add noise to the image
        smallImage.setCenter(ix,iy)
        bounds = smallImage.bounds & bigImage.bounds
        bigImage[bounds] += smallImage[bounds]

    # Record the catalog of generated objects.
    if opts.ngal > 0:
        writeFitsCatalog(inputCatalog,opts.CatalogInFile) 
    # Write the subImage file.
    #bigImage.write(opts.OutputFile)
    writeFitsImage(bigImage,opts.OutputFile,wcs)
    
    # Resize the weight map
    imageRegion = (subRegion[0], bigImage.array.shape[1]+subRegion[0], subRegion[1], bigImage.array.shape[0]+subRegion[2])
    subWeight, Wcent, Ewcs = getBigImage(opts.WeightMapIn,subRegion=imageRegion)


    subWeight.write(opts.WeightMapOut)
    
    # Runs SourceExtractor, which generates catalogs from the image.
    eng = sextractor_engine.SextractorEngine(IMAGE=opts.OutputFile,
                                             WEIGHT_IMAGE=opts.WeightMapOut,
                                             CHECKIMAGE_TYPE='SEGMENTATION,BACKGROUND',
                                             CATALOG_NAME=opts.CatalogOutFile)
    eng.run()
