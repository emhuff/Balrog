#!/usr/bin/env python

import copy
import sys
import time
import os
import subprocess
import argparse
import numpy as np
import astropy.io.fits as pyfits
import galsim
import galsim.des
import sextractor_engine
from model_class import *
from config import *



def WriteCatalog(sample, outfile, cmdline_opts):
    columns = []
    for key in sample.galaxy.keys():
        name = '%s' %(key)
        arr = sample.galaxy[key]
        if key=='x':
            arr = arr - cmdline_opts.xmin + 1
        if key=='y':
            arr = arr - cmdline_opts.ymin + 1
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
    tbhdu.header['XSTART'] = cmdline_opts.xmin
    tbhdu.header['XEND'] = cmdline_opts.xmax
    tbhdu.header['YSTART'] = cmdline_opts.ymin
    tbhdu.header['YEND'] = cmdline_opts.ymax
    tbhdu.header['NSIM'] = cmdline_opts.ngal

    phdu = pyfits.PrimaryHDU()
    hdus = pyfits.HDUList([phdu,tbhdu])
    if os.path.exists(outfile):
        subprocess.call(['rm',outfile])
    hdus.writeto(outfile)

    if derived_opts.assoc!=None:
        data = tbhdu.data
        d = []
        for name in data.columns.names:
            d.append( data[name] )
        d = tuple(d)
        np.savetxt(derived_opts.assoc, np.dstack(d)[0], fmt='%.5f')
        return data.columns.names

    else:
        return None


def CopyAssoc(cmdline_opts, assocnames):
    mhdus = pyfits.open(derived_opts.catalogmeasured, mode='update')
    mhead = mhdus[2].header
    for i in range(len(assocnames)):
        mhead[ 'V%i'%i ] = assocnames[i]
    mhdus.close() 


'''
def WriteImage(image,outFile,header=None):
    imArr = image.array
    hdu = pyfits.PrimaryHDU(imArr,header=header)
    if os.path.exists(outFile):
        subprocess.call(['rm', outFile])
    hdu.writeto(outFile)


def ReadImage(cmdline_opts,file,ext): 
    bigImage = galsim.fits.read(file, hdu=ext)
    if bigImage.wcs==galsim.PixelScale(1):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'fiducialwcs.fits')
        bigImage.wcs = galsim.GSFitsWCS(file)
        cmdline_opts.wcshead = file
    cmdline_opts.galsimwcs = bigImage.wcs

    subBounds = galsim.BoundsI(cmdline_opts.xmin,cmdline_opts.xmax,cmdline_opts.ymin,cmdline_opts.ymax)
    bigImage = bigImage[subBounds]

    hdulist = pyfits.open(file)
    try:
        hdulist[0].header['CRPIX1'] -= cmdline_opts.xmin - 1
        hdulist[0].header['CRPIX2'] -= cmdline_opts.ymin - 1
    except:
        pass
    header = hdulist[0].header
    hdulist.close()
    
    return bigImage, header
'''


def ReadImages(cmdline_opts, derived_opts):
    image = galsim.fits.read(cmdline_opts.imagein, hdu=cmdline_opts.imageext)
    weight = galsim.fits.read(cmdline_opts.weightin, hdu=cmdline_opts.weightext)
    if image.wcs==galsim.PixelScale(1):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'fiducialwcs.fits')
        image.wcs = galsim.GSFitsWCS(file)
        weight.wcs = image.wcs
        derived_opts.wcshead = file
    derived_opts.galsimwcs = image.wcs

    subBounds = galsim.BoundsI(cmdline_opts.xmin,cmdline_opts.xmax,cmdline_opts.ymin,cmdline_opts.ymax)
    image = image[subBounds]
    weight = weight[subBounds]
    psfmodel = galsim.des.DES_PSFEx(cmdline_opts.psfin, derived_opts.wcshead)

    return image, weight, psfmodel


def WriteImages(cmdline_opts, derived_opts, image, weight, nosim=False):
    if nosim:
        imageout = derived_opts.nosim_imageout
        weightout = derived_opts.nosim_weightout
    else:
        imageout = derived_opts.imageout
        weightout = derived_opts.weightout

    if weightout==imageout:
        galsim.fits.writeMulti(image_list=[image,weight], file_name=imageout)
    else:
        galsim.fits.write(image=image, file_name=imageout)
        galsim.fits.write(image=weight, file_name=weightout)

    if not derived_opts.psf_written:
        WritePsf(cmdline_opts, cmdline_opts.psfin, derived_opts.psfout)
        opts.psf_written = True


def WritePsf(cmdline_opts, psfin, psfout):
    psfhdus = pyfits.open(psfin)
    psfhdus[1].header['POLZERO1'] = psfhdus[1].header['POLZERO1'] - (cmdline_opts.xmin - 1)
    psfhdus[1].header['POLZERO2'] = psfhdus[1].header['POLZERO2'] - (cmdline_opts.ymin - 1)
    if os.path.exists(psfout):
        subprocess.call(['rm', psfout])
    psfhdus.writeto(psfout)



def InsertSimulatedGalaxies(bigImage, simulatedgals, psizes, psfmodel, cmdline_opts, derived_opts):
    for i in range(cmdline_opts.ngal):
        postageStampSize = int(psizes[i])
        combinedObjConv = simulatedgals.GetPSFConvolved(derived_opts, psfmodel, i)

        ix = int(simulatedgals.galaxy['x'][i])
        iy = int(simulatedgals.galaxy['y'][i])
        smallImage = galsim.Image(postageStampSize,postageStampSize)
        smallImage.setCenter(ix,iy)
        smallImage.wcs = bigImage.wcs
        smallImage = combinedObjConv.draw(smallImage)
        smallImage.addNoise(galsim.CCDNoise(gain=cmdline_opts.gain,read_noise=0))
        bounds = smallImage.bounds & bigImage.bounds
        bigImage[bounds] += smallImage[bounds]

    return bigImage


def RunSextractor(cmdline_opts, derived_opts, assocnames, nosim=False):
    if nosim:
        catalogmeasured = derived_opts.nosim_catalogmeasured
        imageout = derived_opts.nosim_imageout
        weightout = derived_opts.nosim_weightout
    else:
        catalogmeasured = derived_opts.catalogmeasured
        imageout = derived_opts.imageout
        weightout = derived_opts.weightout

    config_file = cmdline_opts.sexconfig
    if not nosim:
        param_file = cmdline_opts.sexparam
    else:
        param_file = cmdline_opts.sexemptyparam

    if derived_opts.assoc!=None:
        pfile = DefaultName(catalogmeasured, '.fits', '.sex.params', derived_opts.sexdir)
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

    eng = sextractor_engine.SextractorEngine(IMAGE='%s[%i],%s[%s]'%(imageout,derived_opts.outimageext,imageout,derived_opts.outimageext),
                                             WEIGHT_IMAGE='%s[%i],%s[%i]'%(weightout,derived_opts.outweightext,weightout,derived_opts.outweightext),
                                             CATALOG_NAME=catalogmeasured, 
                                             c=config_file,
                                             PARAMETERS_NAME=param_file,
                                             STARNNW_NAME=cmdline_opts.sexnnw,
                                             FILTER_NAME=cmdline_opts.sexconv)
    eng.Path(cmdline_opts.sexpath)
    eng.config['MAG_ZEROPOINT'] = cmdline_opts.zeropoint
    eng.config['PSF_NAME'] = '%s,%s' %(derived_opts.psfout,derived_opts.psfout)

    if derived_opts.assoc!=None:
        ind = range(1, len(assocnames)+1)
        inds = []
        for i in ind:
            inds.append(str(i))
            if assocnames[i-1] == 'x':
                x = i
            if assocnames[i-1] == 'y':
                y = i
        eng.config['ASSOC_NAME'] = derived_opts.assoc
        eng.config['ASSOC_PARAMS'] = '%i,%i' %(x,y)
        eng.config['ASSOC_DATA'] = ','.join(inds)
        eng.config['ASSOC_RADIUS'] = '2.0'
        eng.config['ASSOC_TYPE'] = 'NEAREST'
        eng.config['ASSOCSELEC_TYPE'] = 'MATCHED'
    
    extraconfig = {}
    SextractorConfigs(derived_opts.user_version, extraconfig)
    for key in extraconfig.keys():
        eng.config[key] = extraconfig[key]

    logfile = DefaultName(catalogmeasured, '.fits', '.log.txt', derived_opts.logdir)
    eng.run(logfile=logfile)


def NosimRunSextractor(cmdline_opts, derived_opts, bigImage, subweight, assocnames):
    if derived_opts.subsample:
        WriteImages(cmdline_opts, derived_opts, bigImage, subWeight, nosim=True)
    else:
        if os.path.lexists(derived_opts.nosim_imageout):
            subprocess.call( ['rm', derived_opts.nosim_imageout] )
        if os.path.lexists(derived_opts.nosim_weightout):
            subprocess.call( ['rm', derived_opts.nosim_weightout] )
        if os.path.lexists(derived_opts.psfout):
            subprocess.call( ['rm', derived_opts.psfout] )

        subprocess.call( ['ln', '-s', cmdline_opts.imagein, derived_opts.nosim_imageout] )
        subprocess.call( ['ln', '-s', cmdline_opts.psfin, derived_opts.psfout] )
        derived_opts.psf_written = True
        if derived_opts.nosim_weightout!=derived_opts.nosim_imageout:
            subprocess.call( ['ln', '-s', cmdline_opts.weightin, derived_opts.nosim_weightout] )

    RunSextractor(cmdline_opts, derived_opts, assocnames, nosim=True)


def Cleanup(cmdline_opts,derived_opts):
    files = [derived_opts.imageout, derived_opts.psfout, derived_opts.weightout, derived_opts.nosim_imageout, derived_opts.nosim_weightout]
    for file in files:
        if os.path.exists(file):
            subprocess.call(['rm',file])


def GetSimulatedGalaxies(cmdline_opts, derived_opts, psfmodel):
    rules = SimRules()
    cmdline_opts_copy = copy.copy(cmdline_opts)
    CustomParseArgs(cmdline_opts_copy)
    SimulationRules(cmdline_opts_copy,rules)
    derived_opts.user_version = cmdline_opts_copy
    simulatedgals = DefineRules(cmdline_opts, x=rules.x, y=rules.y, g1=rules.g1, g2=rules.g2, magnification=rules.magnification, nProfiles=rules.nProfiles, axisratio=rules.axisratio, beta=rules.beta, halflightradius=rules.halflightradius, magnitude=rules.magnitude, sersicindex=rules.sersicindex)
    return DoSampling(cmdline_opts,derived_opts,simulatedgals,psfmodel)


def DoSampling(cmdline_opts, derived_opts, simulatedgals, psfmodel):
    psizes = simulatedgals.Sample(cmdline_opts, derived_opts, psfmodel)
    return (simulatedgals,psizes)



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


class DerivedArgs():
    def __init__(self,args):
        self.imgdir = os.path.join(args.outdir, 'balrog_image')
        self.catdir = os.path.join(args.outdir, 'balrog_cat')
        self.logdir = os.path.join(args.outdir, 'balrog_log')
        self.sexdir = os.path.join(args.outdir, 'balrog_sexconfig')

        self.imageout = DefaultName(args.imagein, '.fits', '.sim.fits', self.imgdir)
        self.weightout = self.imageout
        if args.weightin!=args.imagein:
            self.weightout = DefaultName(args.weightin, '.fits', '.weight.sim.fits', self.imgdir)
        self.psfout = DefaultName(args.psfin, '.psf', '.psf', self.imgdir)
        self.catalogtruth = DefaultName(args.imagein, '.fits', '.truthcat.sim.fits', self.catdir)
        self.catalogmeasured = DefaultName(args.imagein, '.fits', '.measuredcat.sim.fits', self.catdir)
        self.assoc = None
        if not args.noassoc:
            self.assoc = DefaultName(args.imagein, '.fits', '.assoc.txt', self.sexdir)

        CreateDir(args.outdir)
        CreateSubDir(self.imgdir)
        CreateSubDir(self.catdir)
        CreateSubDir(self.logdir)
        CreateSubDir(self.sexdir)

        #self.frame = 0
        self.psf_written = False
        self.wcshead = args.imagein
        length = len('.sim.fits')
        ext = '.nosim.fits'
        self.nosim_imageout = '%s%s' %(self.imageout[:-length],ext)
        self.nosim_weightout = '%s%s' %(self.weightout[:-length],ext)
        self.nosim_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)

        self.outimageext = 0
        self.outweightext = 0
        if self.weightout==self.imageout:
            self.outweightext = self.outimageext + 1

        self.subsample = True
        if args.xmin==1 and args.ymin==1 and args.xmax==pyfits.open(args.imagein)[0].header['NAXIS1'] and args.ymax == pyfits.open(args.imagein)[0].header['NAXIS2']:
            self.subsample = False



def GetOpts():
    parser = argparse.ArgumentParser()
    DefaultArgs(parser)
    CustomArgs(parser) 

    cmdline_args = parser.parse_args()
    ParseDefaultArgs(cmdline_args)
    derived_args = DerivedArgs(cmdline_args)
    #CustomParseArgs(args)

    return [cmdline_args, derived_args]


def DefaultName(startfile, lookfor, replacewith, outdir):
    file = os.path.basename(startfile)
    file = os.path.join(outdir, file)
    length = len(lookfor)
    if file.endswith(lookfor):
        fstr = file[:-length]
    else:
        fstr = file
    return '%s%s' %(fstr, replacewith)


def CreateSubDir(subdir):
    if not os.path.exists(subdir):
        subprocess.call(['mkdir', subdir])


def CreateDir(dir):
    full = False
    while dir[0]=='/':
        dir = dir[1:]
        full = True
    while dir[-1]=='/':
        dir = dir[-1]
    dirs = dir.strip().split('/')
    if full:
        subdir = '/'
    else:
        subdir = './'

    for dir in dirs:
        subdir = os.path.join(subdir,dir)
        if not os.path.exists(subdir):
            subprocess.call( ['mkdir', subdir] )


def ParseDefaultArgs(args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    indir = os.path.join(defdir, 'input')
    outdir = os.path.join(defdir, 'output')
    configdir = os.path.join(thisdir, 'astro_config')

    if args.sexconfig==None:
        args.sexconfig = os.path.join(configdir, 'sex.config')
    if args.sexparam==None:
        args.sexparam = os.path.join(configdir, 'bulge.param')
    if args.sexemptyparam==None:
        args.sexemptyparam = os.path.join(configdir, 'sex.param')
    if args.sexnnw==None:
        args.sexnnw = os.path.join(configdir, 'sex.nnw')
    if args.sexconv==None:
        args.sexconv = os.path.join(configdir, 'sex.conv')

    if args.outdir==None:
        args.outdir = outdir
    if args.imagein==None:
        args.imagein = os.path.join(indir, 'example.fits')
    if args.weightin==None:
        args.weightin = args.imagein
        if args.weightext == None:
            args.weightext = args.imageext + 1
    if args.weightext==None:
        args.weightext = 0
    if args.psfin==None:
        args.psfin = os.path.join(indir, 'example.psf')


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

    return args



def DefaultArgs(parser):
    # Input and (temporary) output Images
    parser.add_argument( "-od", "--outdir", help="Directory where to put output. By default, output files will be named automatically based on the input file name. Each default filename can overridden with the appropriate command line option.", default=None)

    parser.add_argument( "-ii", "--imagein", help="Input image to draw simulated galaxies into", type=str, default=None)
    parser.add_argument( "-ie", "--imageext", help="FITS extension where the image lives in the input file", type=int, default=0)

    parser.add_argument( "-wi", "--weightin", help="Weight map of input image", type=str, default=None)
    parser.add_argument( "-we", "--weightext", help="FITS extension where the weight map lives in the input weight file", type=int, default=None)
    
    parser.add_argument( "-pi", "--psfin", help="PSF of thin input image, to be convolved with simulated galaxies", type=str, default=None)
    parser.add_argument( "-c", "--clean", help="Delete output image, weight, and PSF files", action="store_true")

    # Properties you want your simulated image to have
    parser.add_argument( "-xmin", "--xmin", help="Minimum column of extracted subimage, indexing ranges from (1,NumPixelsX)", type=int, default=1)
    parser.add_argument( "-xmax", "--xmax", help="Maximum column of extracted subimage, indexing ranges from (1,NumPixelsX)", type=int, default=-1)
    parser.add_argument( "-ymin", "--ymin", help="Minimum row of extracted subimage, indexing ranges from (1,NumPixelsY)", type=int, default=1)
    parser.add_argument( "-ymax", "--ymax", help="Maximum row of extracted subimage, indexing ranges from (1,NumPixelsY)", type=int, default=-1)
    parser.add_argument( "-ngal", "--ngal", help="Number of simulated galaxies", type=int, default=50)
    parser.add_argument( "-gain", "--gain", help="Gain, needed for adding noise. Can be a float, or a keyword from the image header. (Default is to read from image header keyword 'GAIN'. If that fails, default gain is set to 1)", default='GAIN')
    parser.add_argument( "-zp", "--zeropoint", help="Zeropoint used to convert simulated magnitude to a simulated image fluxes. Sextractor will run also with this zeropoint. Default = 30", type=float, default=30)

    ##### Parameters when placing simulated galaxies into the images.
    parser.add_argument( "-ft", "--fluxthresh", help="Flux value where to cutoff the postage stamp", type=float, default=0.01)
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




if __name__ == "__main__":

    # Parse command line options
    cmdline_opts, derived_opts = GetOpts()

    # Get the subsampled flux and weightmap images, along with the PSF model
    bigImage, subWeight, psfmodel = ReadImages(cmdline_opts, derived_opts)

    # Get simulated galaxy sample. Write it to a truth catalog. If associating, get a list keeping track of what quantity is which index in the output vector.
    catalog, psizes = GetSimulatedGalaxies(cmdline_opts, derived_opts, psfmodel)
    assocnames = WriteCatalog(catalog, derived_opts.catalogtruth, cmdline_opts)

    # Run sextractor over the image without any simulated galaxies. This is to make sure no simulated galaxies "found" are actually brigter, larger galaxies in the data image prior to simulation.
    if not cmdline_opts.noempty:
        NosimRunSextractor(cmdline_opts, derived_opts, bigImage, subWeight, assocnames)

    # Insert simulated galaxies. Write out the flux and weight images with simulated galaxies in them.
    bigImage = InsertSimulatedGalaxies(bigImage, catalog, psizes, psfmodel, cmdline_opts, derived_opts)
    WriteImages(cmdline_opts, derived_opts, bigImage, subWeight)

    #Run sextractor over the simulated image. Write association parameter labels to measured catalog if necessary
    RunSextractor(cmdline_opts, derived_opts, assocnames)
    if derived_opts.assoc!=None:
        CopyAssoc(cmdline_opts, assocnames)

    #If chosen, clean up image files you don't need anymore
    if cmdline_opts.clean:
        Cleanup(cmdline_opts, derived_opts)

