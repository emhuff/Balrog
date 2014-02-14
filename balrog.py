#!/usr/bin/env python

import sys
import time
import os
import subprocess
import argparse
import numpy as np
import pyfits
import galsim
import sextractor_engine
from model_class import *
from config import *



def WriteCatalog(sample, outfile, opts):
    columns = []
    for key in sample.galaxy.keys():
        name = '%s' %(key)
        arr = sample.galaxy[key]
        if key=='x':
            arr = arr - opts.xmin + 1
        if key=='y':
            arr = arr - opts.ymin + 1
        col = pyfits.Column(name=name, array=arr,format='E')
        columns.append(col)
    for i in range(len(sample.component)):
        for key in sample.component[i].keys():
            name = '%s_%i' %(key,i)
            if key.find('halflightradius')!=-1:
                col = pyfits.Column(name=name, array=sample.component[i][key]/np.sqrt(sample.component[i]['axisratio']), format='E')
            else:
                col = pyfits.Column(name=name, array=sample.component[i][key],format='E')
            columns.append(col)
    tbhdu = pyfits.new_table(pyfits.ColDefs(columns))
    tbhdu.header['XSTART'] = opts.xmin
    tbhdu.header['XEND'] = opts.xmax
    tbhdu.header['YSTART'] = opts.ymin
    tbhdu.header['YEND'] = opts.ymax
    tbhdu.header['NSIM'] = opts.ngal

    phdu = pyfits.PrimaryHDU()
    hdus = pyfits.HDUList([phdu,tbhdu])
    if os.path.exists(outfile):
        subprocess.call(['rm',outfile])
    hdus.writeto(outfile)

    if opts.assoc!=None:
        data = tbhdu.data
        d = []
        for name in data.columns.names:
            d.append( data[name] )
        d = tuple(d)
        np.savetxt(opts.assoc, np.dstack(d)[0], fmt='%.5f')
        return data.columns.names

    else:
        return None


def CopyAssoc(opts, assocnames):
    mhdus = pyfits.open(opts.catalogmeasured, mode='update')
    mhead = mhdus[2].header
    for i in range(len(assocnames)):
        mhead[ 'V%i'%i ] = assocnames[i]
    mhdus.close() 


def WriteImage(image,outFile,header=None):
    imArr = image.array
    hdu = pyfits.PrimaryHDU(imArr,header=header)
    if os.path.exists(outFile):
        subprocess.call(['rm', outFile])
    hdu.writeto(outFile)


def ReadImage(opts,file): 
    bigImage = galsim.fits.read(file)
    bigImage.setScale(opts.pixscale)
    subBounds = bigImage.bounds
    subBounds.xmin = opts.xmin
    subBounds.xmax = opts.xmax
    subBounds.ymin = opts.ymin
    subBounds.ymax = opts.ymax
    bigImage = bigImage[subBounds]

    hdulist = pyfits.open(file)
    try:
        hdulist[0].header['CRPIX1'] -= opts.xmin - 1
        hdulist[0].header['CRPIX2'] -= opts.ymin - 1
    except:
        pass
    header = hdulist[0].header
    hdulist.close()
    
    return bigImage, header


def WritePsf(opts, psfin, psfout):
    psfhdus = pyfits.open(psfin)
    psfhdus[1].header['POLZERO1'] = psfhdus[1].header['POLZERO1'] - (opts.xmin - 1)
    psfhdus[1].header['POLZERO2'] = psfhdus[1].header['POLZERO2'] - (opts.ymin - 1)
    if os.path.exists(psfout):
        subprocess.call(['rm', psfout])
    psfhdus.writeto(psfout)


def ParseArgs(args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    configdir = os.path.join(thisdir, 'astro_config')

    if args.sexconfig==None:
        args.sexconfig = os.path.join(configdir, 'sex.config')
    if args.sexparam==None:
        args.sexparam = os.path.join(configdir, 'bulge.param')
        #args.sexparam = os.path.join(configdir, 'disk.param')
        #args.sexparam = os.path.join(configdir, 'bulge+disk.param')
    if args.sexemptyparam==None:
        args.sexemptyparam = os.path.join(configdir, 'sex.param')
    if args.sexnnw==None:
        args.sexnnw = os.path.join(configdir, 'sex.nnw')
    if args.sexconv==None:
        args.sexconv = os.path.join(configdir, 'sex.conv')

    if args.imagein==None:
        args.imagein = os.path.join(defdir, 'example.fits')
    if args.imageout==None:
        args.imageout = args.imagein.replace('.fits', '.out.fits')

    if args.weightin==None:
        args.weightin = os.path.join(defdir, 'example_weight.fits')
    if args.weightout==None:
        args.weightout = args.weightin.replace('.fits', '.out.fits')

    if args.psfin==None:
        args.psfin = os.path.join(defdir, 'example.psf')
    if args.psfout==None:
        args.psfout = args.psfin.replace('.psf', '.out.psf')

    if not args.noassoc:
        if args.assoc==None:
            args.assoc = os.path.join(defdir, 'example_assoc.txt')

    if args.catalogtruth==None:
        args.catalogtruth = os.path.join(defdir, 'truthcat.fits')
    if args.catalogmeasured==None:
        args.catalogmeasured = os.path.join(defdir, 'measuredcat.fits')

    if (args.xmax==-1) and (args.ymax==-1) and (args.xmin==1) and (args.ymin==1):
        args.subsample = False
    else:
        args.subsample = True

    if args.xmax==-1:
        args.xmax = pyfits.open(args.imagein)[0].header['NAXIS1']
    if args.ymax==-1:
        args.ymax = pyfits.open(args.imagein)[0].header['NAXIS2']

    try:
        args.gain = float(args.gain)
    except:
        try:
            args.gain = pyfits.open(args.imagein)[0].header[args.gain]
        except:
            args.gain = 1.0

    args.nosim_imageout = args.imageout.replace('.fits','.nosim.fits')
    args.nosim_weightout = args.weightout.replace('.fits','.nosim.fits')
    args.nosim_catalogmeasured = args.catalogmeasured.replace('.fits','.nosim.fits')

    args.frame = 0

    return args


def InsertSimulatedGalaxies(bigImage, simulatedgals, psizes, psfmodel, opts):
    for i in range(opts.ngal):
        postageStampSize = int(psizes[i])
        combinedObj = simulatedgals.GetPSFConvolved(opts, psfmodel, i)
        pix = galsim.Pixel(opts.pixscale)
        combinedObjConv = galsim.Convolve([pix,combinedObj])

        smallImage = galsim.ImageD(postageStampSize,postageStampSize)
        smallImage = combinedObjConv.draw(dx=opts.pixscale,image=smallImage)
        smallImage.addNoise(galsim.CCDNoise(gain=opts.gain,read_noise=0))

        ix = int(simulatedgals.galaxy['x'][i])
        iy = int(simulatedgals.galaxy['y'][i])
        smallImage.setCenter(ix,iy)

        bounds = smallImage.bounds & bigImage.bounds
        bigImage[bounds] += smallImage[bounds]

    return bigImage



def RunSextractor(opts, assocnames, nosim=False):
    if nosim:
        catalogmeasured = opts.nosim_catalogmeasured
        imageout = opts.nosim_imageout
        weightout = opts.nosim_weightout
    else:
        catalogmeasured = opts.catalogmeasured
        imageout = opts.imageout
        weightout = opts.weightout

    config_file = opts.sexconfig
    if not nosim:
        param_file = opts.sexparam
    else:
        param_file = opts.sexemptyparam

    if opts.assoc!=None:
        pfile = catalogmeasured.replace('.fits', '_sex.params')
        txt = open(param_file).read().strip()

        lines = txt.split('\n')
        todelete = []
        for i in range(len(lines)):
            line = lines[i]
            if line=='':
                continue
            line = line.strip()
            if line=='':
                continue
            if line[0] =='#':
                continue
            if line.startswith('VECTOR_ASSOC('):
                todelete.append(i)
        lines = np.array(lines)
        lines = np.delete(lines, todelete)
        txt = '\n'.join(lines)
        start = 'VECTOR_ASSOC(%i)' %(len(assocnames))
        txt = '%s\n%s' %(start,txt)
        #txt = txt.replace('#VECTOR_ASSOC(2)', 'VECTOR_ASSOC(%i)' %(len(assocnames)) )

        stream = open(pfile, 'w')
        stream.write(txt)
        stream.close()
        param_file = pfile

    eng = sextractor_engine.SextractorEngine(IMAGE='%s,%s'%(imageout,imageout),
                                             WEIGHT_IMAGE='%s,%s'%(weightout,weightout),
                                             CATALOG_NAME=catalogmeasured, 
                                             c=config_file,
                                             PARAMETERS_NAME=param_file,
                                             STARNNW_NAME=opts.sexnnw,
                                             FILTER_NAME=opts.sexconv)
    eng.Path(opts.sexpath)
    eng.config['MAG_ZEROPOINT'] = opts.zeropoint
    eng.config['PSF_NAME'] = '%s,%s' %(opts.psfout,opts.psfout)

    if opts.assoc!=None:
        ind = range(1, len(assocnames)+1)
        inds = []
        for i in ind:
            inds.append(str(i))
            if assocnames[i-1] == 'x':
                x = i
            if assocnames[i-1] == 'y':
                y = i
        eng.config['ASSOC_NAME'] = opts.assoc
        eng.config['ASSOC_PARAMS'] = '%i,%i' %(x,y)
        eng.config['ASSOC_DATA'] = ','.join(inds)
        eng.config['ASSOC_RADIUS'] = '2.0'
        eng.config['ASSOC_TYPE'] = 'NEAREST'
        eng.config['ASSOCSELEC_TYPE'] = 'MATCHED'
    
    extraconfig = {}
    SextractorConfigs(opts, extraconfig)
    for key in extraconfig.keys():
        eng.config[key] = extraconfig[key]

    eng.run()


def Cleanup(opts):
    files = [opts.imageout, opts.psfout, opts.weightout, opts.nosim_imageout, opts.nosim_weightout]
    for file in files:
        if os.path.exists(file):
            subprocess.call(['rm',file])


def GetSimulatedGalaxies(opts):
    rules = SimRules()
    SimulationRules(opts,rules)
    simulatedgals = DefineRules(opts, x=rules.x, y=rules.y, g1=rules.g1, g2=rules.g2, magnification=rules.magnification, nProfiles=rules.nProfiles, axisratio=rules.axisratio, beta=rules.beta, halflightradius=rules.halflightradius, magnitude=rules.magnitude, sersicindex=rules.sersicindex)
    return DoSampling(opts,simulatedgals)
 

def DoSampling(opts, simulatedgals):
    psfmodel = galsim.des.DES_PSFEx(opts.psfin)
    WritePsf(opts, opts.psfin, opts.psfout)
    psizes = simulatedgals.Sample(opts, psfmodel)
    #print np.average(psizes)
    return (simulatedgals,psizes,psfmodel)


class SimRules():
    def __init__(self):
        self.x = None
        self.y = None
        self.g1 = None
        self.g2 = None
        self.magnification = None

        self.nProfiles = 1
        self.axisratio = [None]
        self.beta = [None]
        self.halflightradius = [None]
        self.magnitude = [None]
        self.sersicindex = [None]


def GetArgs():
    parser = argparse.ArgumentParser()

    # Input and (temporary) output Images
    parser.add_argument( "-ii", "--imagein", help="Image to put simulated galaxies into", type=str, default=None)
    parser.add_argument( "-io", "--imageout", help="Output image with simulated galaxies in it", type=str, default=None)
    parser.add_argument( "-wi", "--weightin", help="Weight map of image simulated galaxies are put into", type=str, default=None)
    parser.add_argument( "-wo", "--weightout", help="Weight map of output image with simulated galaxies in it", type=str, default=None)
    parser.add_argument( "-pi", "--psfin", help="PSF of original image without simulated galaxies", type=str, default=None)
    parser.add_argument( "-po", "--psfout", help="Modified psf with different coordinate keywords appropriopriate for (smaller) image with simulated galaxies", type=str, default=None)
    parser.add_argument( "-c", "--clean", help="Delete output images", action="store_true")

    # Output catalogs
    parser.add_argument( "-ct", "--catalogtruth", help="Catalog with parameters of simulated galaxies", type=str, default=None)
    parser.add_argument( "-cm", "--catalogmeasured", help="Catalog measured by sextractor from image with simulated galaxies in it", type=str, default=None)

    # Properties you want your simulated image to have
    parser.add_argument( "-xmin", "--xmin", help="Minimum column of extracted subimage, starts at 1. Default uses entire image.", type=int, default=1)
    parser.add_argument( "-xmax", "--xmax", help="Maximum column of extracted subimage. Default uses entire image.", type=int, default=-1)
    parser.add_argument( "-ymin", "--ymin", help="Minimum row of extracted subimage, starts at 1. Default uses entire image.", type=int, default=1)
    parser.add_argument( "-ymax", "--ymax", help="Maximum row of extracted subimage. Default uses entire image.", type=int, default=-1)
    parser.add_argument( "-ngal", "--ngal", help="Number of simulated galaxies", type=int, default=50)
    parser.add_argument( "-ps", "--pixscale", help="Pixel scale in arcsec/pixel", type=float, default=0.263)
    parser.add_argument( "-gain", "--gain", help="Gain, needed for adding noise. Can be a float, or a keyword from the image header. (Default is to read from image header keyword 'GAIN'. If that fails, default gain is set to 1)", default='GAIN')
    parser.add_argument( "-zp", "--zeropoint", help="Zeropoint used to convert simulated magnitude to a simulated image fluxes. Sextractor will run also with this zeropoint. Default = 30", type=float, default=30)

    ##### Parameters when placing simulated galaxies into the images.
    parser.add_argument( "-ft", "--fluxthresh", help="Flux value where to cutoff the postage stamp", type=float, default=0.1)
    parser.add_argument( "-inc", "--inc", help="Increment if postage stamp size needs to be increased", type=int, default=50)
    parser.add_argument( "-ms", "--minsize", help="Minimum postage stamp size", type=int, default=100)
    #parser.add_argument( "-inc", "--inc", help="Increment if postage stamp size needs to be increased", type=int, default=50)
    #parser.add_argument( "-ms", "--minsize", help="Minimum postage stamp size", type=int, default=80)

    # How to run sextractor
    parser.add_argument( "-spp", "--sexpath", help='Path for sextractor binary', type=str, default='sex')
    parser.add_argument( "-sc", "--sexconfig", help='Sextractor config file', type=str, default=None)
    parser.add_argument( "-sp", "--sexparam", help='Sextractor param file', type=str, default=None)
    parser.add_argument( "-sn", "--sexnnw", help='Sextractor neural network S/G file', type=str, default=None)
    parser.add_argument( "-sv", "--sexconv", help='Sextractor filter convolution file', type=str, default=None)

    parser.add_argument( "-ne", "--noempty", help="Don't do sextractor run that doesn't have simulated galaxies.", action="store_true")
    parser.add_argument( "-sep", "--sexemptyparam", help='Sextractor param file used when image without simulated galaxies is extracted. Which parameters to extract is basically irrelevant, since all the run is inteded to do check for things in image before simulating. Not doing model fitting is faster, so the default is to use one that does not do model fitting.', type=str, default=None)

    parser.add_argument( "-na", "--noassoc", help="Don't do association mode matching in sextractor", action="store_true")
    parser.add_argument( "-as", "--assoc", help='Unless you flag --noassoc sextractor runs in associaion mode, extracting only the simulated positions. The association list file is saved as the given string.', type=str, default=None)

    CustomArguments(parser) 
    args = parser.parse_args()

    args = ParseArgs(args)
    CustomParseArgs(args)
    return args



if __name__ == "__main__":

    # Parse command line options
    opts = GetArgs()

    # Get the subsampled flux and weightmap images
    bigImage, header = ReadImage(opts, opts.imagein)
    subWeight, Wheader = ReadImage(opts, opts.weightin)

    # Get simulated galaxy sample. Write it to a truth catalog. If associating, get a list keeping track of what quantity is which index in the output vector.
    catalog,psizes,psfmodel = GetSimulatedGalaxies(opts)
    assocnames = WriteCatalog(catalog, opts.catalogtruth, opts)

    # Run sextractor over the image without any simulated galaxies. This is to make sure no simulated galaxies "found" are actually brigter, larger galaxies in the data image prior
    # to simulation. The match will be done between the nosim catalog and the simulation catalog later (in concatenate.py)
    if not opts.noempty:
        if opts.subsample:
            WriteImage(bigImage,opts.nosim_imageout,header)
            WriteImage(subWeight,opts.nosim_weightout,Wheader)
        else:
            if os.path.lexists(opts.nosim_imageout):
                subprocess.call( ['rm', opts.nosim_imageout] )
            if os.path.lexists(opts.nosim_weightout):
                subprocess.call( ['rm', opts.nosim_weightout] )
            subprocess.call( ['ln', '-s', opts.imagein, opts.nosim_imageout] )
            subprocess.call( ['ln', '-s', opts.weightin, opts.nosim_weightout] )
        RunSextractor(opts, assocnames, nosim=True)

    # Insert simulated galaxies. Write out the flux and weight images with simulated galaxies in them.
    bigImage = InsertSimulatedGalaxies(bigImage, catalog, psizes, psfmodel, opts)
    WriteImage(bigImage,opts.imageout,header)
    WriteImage(subWeight,opts.weightout,Wheader)

    #Run sextractor over the simulated image. Write association parameter labels to measured catalog if necessary
    RunSextractor(opts, assocnames)
    if opts.assoc!=None:
        CopyAssoc(opts, assocnames)

    #If chosen, Clean up image files you don't need anymore
    if opts.clean:
        Cleanup(opts)

