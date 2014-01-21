#!/usr/bin/env python

import sys
import os
import subprocess
import argparse

import galsim
import numpy as np
import pyfits
import pywcs

import sextractor_engine
from model_class import *


def defineCalibration():
    '''
    This sets the number of counts per unit flux. Used to set the
    noise of the added galaxies.
    '''
    calibParams = {}
    calibParams['gain'] = 10000.0
    calibParams['pixel_scale'] = 0.263
    return calibParams


def WriteComponents(sample, outfile, subRegion):
    columns = []
    for key in sample.galaxy.keys():
        name = '%s' %(key)
        arr = sample.galaxy[key]
        if key=='x':
            arr = arr - subRegion[0] + 1
        if key=='y':
            arr = arr - subRegion[2] + 1
        col = pyfits.Column(name=name, array=arr,format='E')
        columns.append(col)
    for i in range(len(sample.component)):
        for key in sample.component[i].keys():
            name = '%s_%i' %(key,i)
            col = pyfits.Column(name=name, array=sample.component[i][key],format='E')
            columns.append(col)
    tbhdu = pyfits.new_table(pyfits.ColDefs(columns))
    phdu = pyfits.PrimaryHDU()
    hdus = pyfits.HDUList([phdu,tbhdu])
    if os.path.exists(outfile):
        subprocess.call(['rm',outfile])
    hdus.writeto(outfile)


def writeFitsImage(image,outFile,wcs=None,header=None):
    '''
    Write a GalSim Image to extension 0 of a .fits header.
    Optionally, write the wcs.
    '''
    imArr = image.array
    #header = wcs.to_header() #L- World Coordinate System. But what is it?
    hdu = pyfits.PrimaryHDU(imArr,header=header)
    if os.path.exists(outFile):
        subprocess.call(['rm', outFile])
    hdu.writeto(outFile)


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

    subBounds = bigImage.bounds
    if subRegion[0] > 0:
        subBounds.xmin = subRegion[0]
    if subRegion[1] > 0:
        subBounds.xmax = subRegion[1]
    if subRegion[2] > 0:
        subBounds.ymin = subRegion[2]
    if subRegion[3] > 0:
        subBounds.ymax = subRegion[3]
	
    bigImage = bigImage[subBounds]
    offset = galsim.PositionD(subRegion[0],subRegion[2])

    hdulist = pyfits.open(file)
    hdulist[0].header['CRPIX1'] -= subRegion[0] - 1
    hdulist[0].header['CRPIX2'] -= subRegion[2] - 1
    wcs = pywcs.WCS(hdulist[0].header)

    header = hdulist[0].header
    hdulist.close()
    
    return bigImage, offset, wcs, header


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument( "-ii", "--imagein", help="image to put simulated galaxies into", type=str, required=True)
    parser.add_argument( "-io", "--imageout", help="output image with simulated galaxies in it", type=str, required=True)
    parser.add_argument( "-wi", "--weightin", help="weight map of image simulated galaxies are put into", type=str, required=True)
    parser.add_argument( "-wo", "--weightout", help="weight map of output image with simulated galaxies in it", type=str, required=True)
    parser.add_argument( "-pi", "--psfin", help="psf of original image", type=str, required=True)
    parser.add_argument( "-po", "--psfout", help="modified psf with different coordinate keywords appropriopriate for smaller image", type=str, required=True)
    parser.add_argument( "-ct", "--catalogtruth", help="catalog with parameters of simulated galaxies", type=str, required=True)
    parser.add_argument( "-cm", "--catalogmeasured", help="catalog measured from image with simulated galaxies in it", type=str, required=True)
    parser.add_argument( "-cs", "--catalogsample", help="catalog used to sample parameter distriubtions from", type=str, default=None)
    parser.add_argument( "-xmin", "--xmin", help="minimum column of extracted subimage, starts at 1", type=int, default=1)
    parser.add_argument( "-xmax", "--xmax", help="maximum column of extracted subimage", type=int, default=-1)
    parser.add_argument( "-ymin", "--ymin", help="minimum row of extracted subimage, starts at 1", type=int, default=1)
    parser.add_argument( "-ymax", "--ymax", help="maximum row of extracted subimage", type=int, default=-1)
    parser.add_argument( "-ngal", "--ngal", help="number of simulated galaxies", type=int, default=50)

    parser.add_argument( "-nt", "--noisetype", help="type of noise to add to image", type=str, default='epoch', choices=['epoch','calib'])
    parser.add_argument( "-ne", "--nepoch", help="~ number of epochs the image is made from", type=float, default=8)

    args = parser.parse_args()
    if args.xmax==-1:
        args.xmax = pyfits.open(args.imagein)[0].header['NAXIS1']
    if args.ymax==-1:
        args.ymax = pyfits.open(args.imagein)[0].header['NAXIS2']

    return args


def GetSimulatedGalaxies(opts):
    #simulatedgals = Exponential(ngal=opts.ngal, catalog=opts.catalogsample)
    simulatedgals = deVaucouleur(ngal=opts.ngal, catalog=opts.catalogsample)
    
    '''
    off = 50
    xrange = np.linspace( subRegion[0]+50, subRegion[1]-50, np.sqrt(opts.ngal))
    yrange = np.linspace( subRegion[2]+50, subRegion[3]-50, np.sqrt(opts.ngal))
    xarr, yarr = np.meshgrid(xrange,yrange)
    xarr = xarr.flatten()
    yarr = yarr.flatten()
    simulatedgals.GalaxyRule(key='x', rule=Rule(type='array', array=xarr) )
    simulatedgals.GalaxyRule(key='y', rule=Rule(type='array', array=yarr) )
    '''

    #simulatedgals.ComponentRule(key='flux', rule=Rule(type='catalog', column='FLUX_AUTO', joint=True) )
    #simulatedgals.ComponentRule(key='flux', rule=Rule(type='value', value=2.0e3) )
    #simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='value', value=2.5 ) )
    #simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='value', value=0.6/pixscale) )
    #simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='catalog', column='FLUX_RADIUS', joint=True) )
    #simulatedgals.ComponentRule(key='axisratio',rule=Rule(type='value', value=1.0 ) )
    #simulatedgals.ComponentRule(key='beta',rule=Rule(type='value', value=0.0 ) )

    simulatedgals.GalaxyRule(key='x', rule=Rule(type='uniform', minimum=subRegion[0], maximum=subRegion[1]) )
    simulatedgals.GalaxyRule(key='y', rule=Rule(type='uniform', minimum=subRegion[2], maximum=subRegion[3]) )
    simulatedgals.ComponentRule(key='flux', rule=Rule(type='uniform', minimum=1.0e2, maximum=1.0e4) )
    #simulatedgals.ComponentRule(key='flux', rule=Rule(type='value', value=5.0e3) )
    simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='uniform', minimum=0.3/pixscale, maximum=1.3/pixscale) )
    #simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='value', value=0.6/pixscale) )
    simulatedgals.ComponentRule(key='axisratio',rule=Rule(type='uniform', minimum=0.33, maximum=1) )
    #simulatedgals.ComponentRule(key='axisratio',rule=Rule(type='value', value=0.9 ) )
    #simulatedgals.ComponentRule(key='beta',rule=Rule(type='uniform', minimum=-90, maximum=90) )
    simulatedgals.ComponentRule(key='beta',rule=Rule(type='value', value=0.0 ) )

    simulatedgals.Sample()
    return simulatedgals


if __name__ == "__main__":

    opts = GetArgs()
    subRegion = (opts.xmin,opts.xmax,opts.ymin,opts.ymax)
    calib = defineCalibration()

    bigImage, offset, wcs, header = getBigImage(opts.imagein,subRegion=subRegion,calibration=calib)
    pixscale = bigImage.getScale()
    psfmodel = galsim.des.DES_PSFEx(opts.psfin)
    bigImage.array[:,:] = 0.0

    psfhdus = pyfits.open(opts.psfin)
    psfhdus[1].header['POLZERO1'] = psfhdus[1].header['POLZERO1'] - (opts.xmin - 1)
    psfhdus[1].header['POLZERO2'] = psfhdus[1].header['POLZERO2'] - (opts.ymin - 1)
    if os.path.exists(opts.psfout):
        subprocess.call(['rm', opts.psfout])
    psfhdus.writeto(opts.psfout)
    
    simulatedgals = GetSimulatedGalaxies(opts)
    
    for i in range(opts.ngal):
        x = float(simulatedgals.galaxy['x'][i])
        y = float(simulatedgals.galaxy['y'][i])
        scales = np.zeros(len(simulatedgals.component))
        for j in range(len(simulatedgals.component)):
            reff = simulatedgals.component[j]['halflightradius'][i] 
            #reff = reff / pow(simulatedgals.component[j]['axisratio'][i],simulatedgals.component[j]['sersicindex'][i])
            scales[j] = reff
            reff = reff * pixscale

            sersicObj = galsim.Sersic(n=float(simulatedgals.component[j]['sersicindex'][i]), half_light_radius=float(reff), flux=float(simulatedgals.component[j]['flux'][i]))
            intrinsic_shear = galsim.Shear(q=simulatedgals.component[j]['axisratio'][i], beta=simulatedgals.component[j]['beta'][i]*galsim.degrees )
            sersicObj.applyShear(intrinsic_shear)

            if j==0:
                combinedObj = sersicObj
            else:
                combinedObj = combinedObj + sersicObj

        lensing_shear = galsim.Shear(g1=simulatedgals.galaxy['g1'][i], g2=simulatedgals.galaxy['g2'][i])
        combinedObj.applyShear(lensing_shear)
        combinedObj.applyMagnification(simulatedgals.galaxy['magnification'][i])

        minsize = 53
        psize = 6.0 * np.amax(scales) * 1.0/simulatedgals.component[j]['axisratio'][i]
        psize = max( psize, minsize )
        psize = int(psize)
        if psize%2==0:
            psize += 1
        postageStampSize = psize

        ix = int(np.floor(x)) 
        iy = int(np.floor(y))
        dx = x-ix
        dy = y-iy
        pos = galsim.PositionD(x,y)
        combinedObj.applyShift(dx*pixscale,dy*pixscale)
        
        psf = psfmodel.getPSF(pos,bigImage.getScale())
        psf.setFlux(1.)
        combinedObj = galsim.Convolve([psf,combinedObj])
        pix = galsim.Pixel(bigImage.getScale())
        combinedObjConv = galsim.Convolve([pix,combinedObj])

        try:
            smallImage = galsim.ImageD(postageStampSize,postageStampSize)
            smallImage = combinedObjConv.draw(dx=bigImage.getScale(),image=smallImage)
            
            if opts.noisetype=='epoch':
                cut = (smallImage.array > 0.001)
                smallImage.array[cut] = np.random.normal( smallImage.array[cut], np.sqrt(smallImage.array[cut]/opts.nepoch) )
                #smallImage.array[:,:] = np.random.normal( smallImage.array[:,:], np.sqrt(smallImage.array[:,:]/opts.nepoch) )
                #smallImage.array[:,:] = np.random.poisson( smallImage.array[:,:] )
            elif opts.noisetype=='calib':
                #smallImage.addNoise(galsim.CCDNoise(gain=calib['gain'],read_noise=0))
                pass

            smallImage.setCenter(ix,iy)
            bounds = smallImage.bounds & bigImage.bounds
            bigImage[bounds] += smallImage[bounds]
        except:
            print 'bad', simulatedgals.component[0]['sersicindex'][i],simulatedgals.component[0]['halflightradius'][i],simulatedgals.component[0]['axisratio'][i], psize
    
    #bigImage.addNoise(galsim.CCDNoise(gain=calib['gain'],read_noise=0))
    shapey = bigImage.array.shape[0]
    shapex = bigImage.array.shape[1]
    #rands = np.random.normal( 0,8.0e-2, shapex*shapey )
    rands = np.random.normal( 0,5.0e-1, shapex*shapey )
    rands = np.reshape( rands, (shapey,shapex) )
    bigImage.array[:,:] = bigImage.array[:,:] + rands

    if opts.ngal > 0:
        WriteComponents(simulatedgals, opts.catalogtruth, subRegion)

    # Write the subImage file.
    writeFitsImage(bigImage,opts.imageout,wcs,header)

    # Next, write the subImage weightmap.
    subWeight, Wcent, Ewcs,Wheader = getBigImage(opts.weightin,subRegion=subRegion)
    writeFitsImage(subWeight,opts.weightout,Ewcs,Wheader)
  
    this_dir = os.path.dirname( os.path.realpath(__file__) )
    config_dir = os.path.join(this_dir, 'astro_config')
    config_file = os.path.join(config_dir, 'catalog_sex.config')
    param_file = os.path.join(config_dir, 'catalog_model_sex.param')
    #param_file = os.path.join(config_dir, 'sex.param')

    eng = sextractor_engine.SextractorEngine(IMAGE='%s,%s'%(opts.imageout,opts.imageout),
                                             WEIGHT_IMAGE='%s,%s'%(opts.weightout,opts.weightout),
                                             CATALOG_NAME=opts.catalogmeasured, 
                                             c=config_file,
                                             PARAMETERS_NAME=param_file)
    eng.config['MAG_ZEROPOINT'] = 30.2087061678777
    eng.config['PSF_NAME'] = '%s,%s' %(opts.psfout,opts.psfout)
    #eng.config['PSF_NAME'] = '%s,%s' %(opts.psfin,opts.psfin)
    #eng.config['BACK_TYPE'] = 'MANUAL'
    #eng.config['BACK_VALUE'] = 0
    #eng.config['DETECT_THRESH'] = 0.5
    #eng.config['ANALYSIS_THRESH'] = 0.5
    #eng.config['CHECKIMAGE_TYPE'] = 'SEGMENTATION,BACKGROUND_RMS,BACKGROUND,APERTURES'
    #eng.config['CHECKIMAGE_NAME'] = 'simulation_files/seg.fits,simulation_files/backrms.fits,simulation_files/back.fits,simulation_files/apps.fits'
    eng.run()
