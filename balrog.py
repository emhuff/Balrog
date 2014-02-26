#!/usr/bin/env python

import copy
import datetime
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



def WriteCatalog(sample, BalrogSetup, txt=None, fits=False):
    columns = []
    for key in sample.galaxy.keys():
        name = '%s' %(key)
        arr = sample.galaxy[key]
        if key=='x':
            arr = arr - BalrogSetup.xmin + 1
        if key=='y':
            arr = arr - BalrogSetup.ymin + 1
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
    tbhdu.header['XSTART'] = BalrogSetup.xmin
    tbhdu.header['XEND'] = BalrogSetup.xmax
    tbhdu.header['YSTART'] = BalrogSetup.ymin
    tbhdu.header['YEND'] = BalrogSetup.ymax
    tbhdu.header['NSIM'] = BalrogSetup.ngal

    if fits:
        phdu = pyfits.PrimaryHDU()
        hdus = pyfits.HDUList([phdu,tbhdu])
        if os.path.exists(BalrogSetup.catalogtruth):
            subprocess.call(['rm',BalrogSetup.catalogtruth])
        hdus.writeto(BalrogSetup.catalogtruth)

    if txt!=None:
        data = tbhdu.data
        d = []
        for name in data.columns.names:
            d.append( data[name] )
        d = tuple(d)
        np.savetxt(txt, np.dstack(d)[0], fmt='%.5f')
        BalrogSetup.assocnames = data.columns.names



def CopyAssoc(BalrogSetup, outfile):
    mhdus = pyfits.open(outfile, mode='update')
    mhead = mhdus[2].header
    for i in range(len(BalrogSetup.assocnames)):
        mhead[ 'V%i'%i ] = BalrogSetup.assocnames[i]
    mhdus.close() 



def ReadImages(BalrogSetup):
    image = galsim.fits.read(BalrogSetup.imagein, hdu=BalrogSetup.imageext)
    weight = galsim.fits.read(BalrogSetup.weightin, hdu=BalrogSetup.weightext)
    if image.wcs==galsim.PixelScale(1):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'fiducialwcs.fits')
        image.wcs = galsim.GSFitsWCS(file)
        weight.wcs = image.wcs
        BalrogSetup.wcshead = file
    wcs = image.wcs

    subBounds = galsim.BoundsI(BalrogSetup.xmin,BalrogSetup.xmax,BalrogSetup.ymin,BalrogSetup.ymax)
    image = image[subBounds]
    weight = weight[subBounds]
    psfmodel = galsim.des.DES_PSFEx(BalrogSetup.psfin, BalrogSetup.wcshead)

    return image, weight, psfmodel, wcs


def WriteImages(BalrogSetup, image, weight, nosim=False):
    if nosim:
        imageout = BalrogSetup.nosim_imageout
        weightout = BalrogSetup.nosim_weightout
    else:
        imageout = BalrogSetup.imageout
        weightout = BalrogSetup.weightout

    if weightout==imageout:
        galsim.fits.writeMulti(image_list=[image,weight], file_name=imageout)
    else:
        galsim.fits.write(image=image, file_name=imageout)
        galsim.fits.write(image=weight, file_name=weightout)

    if not BalrogSetup.psf_written:
        WritePsf(BalrogSetup, BalrogSetup.psfin, BalrogSetup.psfout)
        opts.psf_written = True


def WritePsf(BalrogSetup, psfin, psfout):
    psfhdus = pyfits.open(psfin)
    psfhdus[1].header['POLZERO1'] = psfhdus[1].header['POLZERO1'] - (BalrogSetup.xmin - 1)
    psfhdus[1].header['POLZERO2'] = psfhdus[1].header['POLZERO2'] - (BalrogSetup.ymin - 1)
    if os.path.exists(psfout):
        subprocess.call(['rm', psfout])
    psfhdus.writeto(psfout)



def InsertSimulatedGalaxies(bigImage, simulatedgals, psfmodel, BalrogSetup, wcs):
    psizes, athresh = simulatedgals.GetPSizes(BalrogSetup, wcs)

    t0 = datetime.datetime.now()
    rt = long( t0.microsecond )

    simulatedgals.galaxy['sum'] = np.copy(simulatedgals.galaxy['x'])
    for i in range(BalrogSetup.ngal):

        postageStampSize = int(psizes[i])
        combinedObjConv = simulatedgals.GetPSFConvolved(psfmodel, i, wcs, athresh)

        ix = int(simulatedgals.galaxy['x'][i])
        iy = int(simulatedgals.galaxy['y'][i])
        smallImage = galsim.Image(postageStampSize,postageStampSize)
        smallImage.setCenter(ix,iy)
        smallImage.wcs = bigImage.wcs
        smallImage = combinedObjConv.draw(smallImage)

        t1 = datetime.datetime.now()
        dt = t1 - t0
        micro = long( (dt.days*24*60*60 + dt.seconds)*1.0e6 + dt.microseconds ) + rt

        smallImage.addNoise(galsim.CCDNoise(gain=BalrogSetup.gain,read_noise=0,rng=galsim.BaseDeviate(rt)))
        bounds = smallImage.bounds & bigImage.bounds

        sum = np.sum(smallImage[bounds].array.flatten())
        simulatedgals.galaxy['sum'][i] = sum

        bigImage[bounds] += smallImage[bounds]

    return bigImage


def IsValidLine(line):
    if line=='':
        return False
    line = line.strip()
    if line=='':
        return False
    if line[0] =='#':
        return False
    return True


def ParamTxtWithoutAssoc(param_file):
    txt = open(param_file).read().strip()

    lines = txt.split('\n')
    todelete = []
    for i in range(len(lines)):
        line = lines[i]
        if not IsValidLine(line):
            continue
        if line.startswith('VECTOR_ASSOC('):
            todelete.append(i)
    lines = np.array(lines)
    lines = np.delete(lines, todelete)
    txt = '\n'.join(lines)
    return txt


def WriteParamFile(BalrogSetup, catalogmeasured, nosim):
    if not nosim:
        param_file = BalrogSetup.sexparam
    else:
        param_file = BalrogSetup.sexemptyparam

    pfile = DefaultName(catalogmeasured, '.fits', '.sex.params', BalrogSetup.sexdir)
    txt = ParamTxtWithoutAssoc(param_file)
    if not BalrogSetup.noassoc:
        start = 'VECTOR_ASSOC(%i)' %(len(BalrogSetup.assocnames))
        txt = '%s\n%s' %(start,txt)
    stream = open(pfile, 'w')
    stream.write(txt)
    stream.close()
    return pfile


def WriteConfigFile(BalrogSetup, config_file, catalogmeasured):
    cfile = DefaultName(catalogmeasured, '.fits', '.sex.config', BalrogSetup.sexdir)
    txt = open(config_file).read().strip()
    lines = txt.split('\n')
    todelete = []
    for i in range(len(lines)):
        line = lines[i]
        if not IsValidLine(line):
            continue
        if line.find('ASSOC')!=-1:
            todelete.append(i)
    if len(todelete)==0:
        return config_file
    lines = np.array(lines)
    lines = np.delete(lines, todelete)
    txt = '\n'.join(lines)
    stream = open(cfile, 'w')
    stream.write(txt)
    stream.close()
    return cfile


def AutoConfig(autologfile, BalrogSetup, imageout, weightout, catalogmeasured, config_file, param_file, afile, eng):
    out = open(autologfile, 'w')
    eng.Path(BalrogSetup.sexpath)

    eng.config['IMAGE'] = '%s[%i],%s[%s]' %(imageout,BalrogSetup.outimageext,imageout,BalrogSetup.outimageext)
    out.write('IMAGE %s[%i],%s[%s]\n' %(imageout,BalrogSetup.outimageext,imageout,BalrogSetup.outimageext) )
    eng.config['WEIGHT_IMAGE'] = '%s[%i],%s[%i]' %(weightout,BalrogSetup.outweightext,weightout,BalrogSetup.outweightext)
    out.write('WEIGHT_IMAGE %s[%i],%s[%i]\n' %(weightout,BalrogSetup.outweightext,weightout,BalrogSetup.outweightext) )
    eng.config['CATALOG_NAME'] = catalogmeasured
    out.write('CATALOG_NAME %s\n' %(catalogmeasured) )
    eng.config['c'] = config_file
    out.write('c %s\n' %(config_file) )
    eng.config['PARAMETERS_NAME'] = param_file
    out.write('PARAMETERS_NAME %s\n' %(param_file) )
    eng.config['STARNNW_NAME'] = BalrogSetup.sexnnw
    out.write('STARNNW_NAME %s\n' %(BalrogSetup.sexnnw) )
    eng.config['FILTER_NAME'] = BalrogSetup.sexconv
    out.write('FILTER_NAME %s\n'  %(BalrogSetup.sexconv) )
    eng.config['MAG_ZEROPOINT'] = BalrogSetup.zeropoint
    out.write('MAG_ZEROPOINT %s\n' %(BalrogSetup.zeropoint) )
    eng.config['PSF_NAME'] = '%s,%s' %(BalrogSetup.psfout,BalrogSetup.psfout)
    out.write('PSF_NAME %s,%s\n' %(BalrogSetup.psfout,BalrogSetup.psfout) )

    if not BalrogSetup.noassoc:
        ind = range(1, len(BalrogSetup.assocnames)+1)
        inds = []
        for i in ind:
            inds.append(str(i))
            if BalrogSetup.assocnames[i-1] == 'x':
                x = i
            if BalrogSetup.assocnames[i-1] == 'y':
                y = i
        eng.config['ASSOC_NAME'] = afile
        out.write('ASSOC_NAME %s\n' %(afile) )
        eng.config['ASSOC_PARAMS'] = '%i,%i' %(x,y)
        out.write('ASSOC_PARAMS %i,%i\n' %(x,y) )
        eng.config['ASSOC_DATA'] = ','.join(inds)
        out.write('ASSOC_DATA %s\n' %(','.join(inds)) )
        eng.config['ASSOC_RADIUS'] = '2.0'
        out.write('ASSOC_RADIUS 2.0\n')
        eng.config['ASSOC_TYPE'] = 'NEAREST'
        out.write('ASSOC_TYPE NEAREST\n')
        eng.config['ASSOCSELEC_TYPE'] = 'MATCHED'
        out.write('ASSOCSELEC_TYPE MATCHED\n')
   
    out.close()


def RunSextractor(BalrogSetup, ExtraSexConfig, catalog, nosim=False):

    if nosim:
        catalogmeasured = BalrogSetup.nosim_catalogmeasured
        imageout = BalrogSetup.nosim_imageout
        weightout = BalrogSetup.nosim_weightout
        autologfile = BalrogSetup.nosim_sexautolog
        logfile = BalrogSetup.nosim_sexlog
        afile = BalrogSetup.assoc_nosimfile
    else:
        catalogmeasured = BalrogSetup.catalogmeasured
        imageout = BalrogSetup.imageout
        weightout = BalrogSetup.weightout
        autologfile = BalrogSetup.sexautolog
        logfile = BalrogSetup.sexlog
        afile = BalrogSetup.assoc_simfile

    if not BalrogSetup.noassoc:
        WriteCatalog(catalog, BalrogSetup, txt=afile, fits=False)

    param_file = WriteParamFile(BalrogSetup, catalogmeasured, nosim)
    config_file = BalrogSetup.sexconfig
    if not BalrogSetup.noassoc:
        config_file = WriteConfigFile(BalrogSetup, config_file, catalogmeasured)

    eng = sextractor_engine.SextractorEngine()
    for key in ExtraSexConfig.keys():
        eng.config[key] = ExtraSexConfig[key]

    AutoConfig(autologfile, BalrogSetup, imageout, weightout, catalogmeasured, config_file, param_file, afile, eng)
    eng.run(logfile=logfile)

    if not BalrogSetup.noassoc:
        CopyAssoc(BalrogSetup, catalogmeasured)


def NosimRunSextractor(BalrogSetup, bigImage, subweight, ExtraSexConfig, catalog):
    if BalrogSetup.subsample:
        WriteImages(BalrogSetup, bigImage, subWeight, nosim=True)
    else:
        if os.path.lexists(BalrogSetup.nosim_imageout):
            subprocess.call( ['rm', BalrogSetup.nosim_imageout] )
        if os.path.lexists(BalrogSetup.nosim_weightout):
            subprocess.call( ['rm', BalrogSetup.nosim_weightout] )
        if os.path.lexists(BalrogSetup.psfout):
            subprocess.call( ['rm', BalrogSetup.psfout] )

        subprocess.call( ['ln', '-s', BalrogSetup.imagein, BalrogSetup.nosim_imageout] )
        subprocess.call( ['ln', '-s', BalrogSetup.psfin, BalrogSetup.psfout] )
        BalrogSetup.psf_written = True
        if BalrogSetup.nosim_weightout!=BalrogSetup.nosim_imageout:
            subprocess.call( ['ln', '-s', BalrogSetup.weightin, BalrogSetup.nosim_weightout] )

    RunSextractor(BalrogSetup, ExtraSexConfig, catalog, nosim=True)


def Cleanup(BalrogSetup):
    files = [BalrogSetup.imageout, BalrogSetup.psfout, BalrogSetup.weightout, BalrogSetup.nosim_imageout, BalrogSetup.nosim_weightout]
    for file in files:
        if os.path.exists(file):
            subprocess.call(['rm',file])


def UserDefinitions(cmdline_args, BalrogSetup):
    rules = SimRules()
    ExtraSexConfig = {}
    cmdline_args_copy = copy.copy(cmdline_args)

    CustomParseArgs(cmdline_args_copy)
    SimulationRules(cmdline_args_copy,rules)
    SextractorConfigs(cmdline_args_copy, ExtraSexConfig)

    LogCmdlineOpts(cmdline_args, cmdline_args_copy, BalrogSetup)
    LogExtraSexConfig(ExtraSexConfig, BalrogSetup)
    rules = DefineRules(BalrogSetup, x=rules.x, y=rules.y, g1=rules.g1, g2=rules.g2, magnification=rules.magnification, nProfiles=rules.nProfiles, axisratio=rules.axisratio, beta=rules.beta, halflightradius=rules.halflightradius, magnitude=rules.magnitude, sersicindex=rules.sersicindex)
    return rules, ExtraSexConfig


def GetSimulatedGalaxies(BalrogSetup, simgals):
    simgals.Sample(BalrogSetup)
    LogSimRules(simgals, BalrogSetup)
    return simgals



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
        if not args.noassoc:
            self.assoc_simfile = DefaultName(args.imagein, '.fits', '.assoc.sim.txt', self.sexdir)
            self.assoc_nosimfile = DefaultName(args.imagein, '.fits', '.assoc.nosim.txt', self.sexdir)

        self.psf_written = False
        self.wcshead = args.imagein
        length = len('.sim.fits')
        ext = '.nosim.fits'
        self.nosim_imageout = '%s%s' %(self.imageout[:-length],ext)
        self.nosim_weightout = '%s%s' %(self.weightout[:-length],ext)
        self.nosim_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)

        self.cmdlinelog = DefaultName(args.imagein, '.fits', '.cmdline_arguments.log.txt', self.logdir)
        self.derivedlog = DefaultName(args.imagein, '.fits', '.derived_arguments.log.txt', self.logdir)
        self.extrasexlog = DefaultName(args.imagein, '.fits', '.sextractor_config_override.log.txt', self.logdir)
        self.sexautolog = DefaultName(args.imagein, '.fits', '.sextractor_config_auto.sim.log.txt', self.logdir)
        self.sexlog = DefaultName(self.catalogmeasured, '.fits', '.log.txt', self.logdir)
        self.nosim_sexautolog = DefaultName(args.imagein, '.fits', '.sextractor_config_auto.nosim.log.txt', self.logdir)
        self.nosim_sexlog = DefaultName(self.nosim_catalogmeasured, '.fits', '.log.txt', self.logdir)
        self.simruleslog = DefaultName(args.imagein, '.fits', '.simulation_rules.log.txt', self.logdir)
        self.catruleslog = DefaultName(args.imagein, '.fits', '.simulationcat_rules.log.txt', self.logdir)

        CreateDir(args.outdir)
        CreateSubDir(self.imgdir)
        CreateSubDir(self.catdir)
        CreateSubDir(self.logdir)
        CreateSubDir(self.sexdir)


        self.outimageext = 0
        self.outweightext = 0
        if self.weightout==self.imageout:
            self.outweightext = self.outimageext + 1

        self.subsample = True
        if args.xmin==1 and args.ymin==1 and args.xmax==pyfits.open(args.imagein)[args.imageext].header['NAXIS1'] and args.ymax==pyfits.open(args.imagein)[args.imageext].header['NAXIS2']:
            self.subsample = False


class BalrogConfig():
    def __init__(self, cargs, dargs):
        cdict = vars(cargs)
        for key in cdict.keys():
            exec "self.%s = cdict['%s']" %(key, key)

        ddict = vars(dargs)
        for key in ddict.keys():
            exec "self.%s = ddict['%s']" %(key, key)


def GetOpts():
    parser = argparse.ArgumentParser()
    DefaultArgs(parser)
    CustomArgs(parser) 

    cmdline_args = parser.parse_args()
    ParseDefaultArgs(cmdline_args)
    return cmdline_args


def ConfigureBalrog(cmdline_opts):
    cmdline_opts, derived_opts = GetMoreOpts(cmdline_opts)
    BalrogSetup = BalrogConfig(cmdline_opts, derived_opts)
    simulation_rules, ExtraSexConfig = UserDefinitions(cmdline_opts, BalrogSetup)
    return BalrogSetup, simulation_rules, ExtraSexConfig


def GetMoreOpts(cmdline_args):
    derived_args = DerivedArgs(cmdline_args)
    return [cmdline_args, derived_args]


def LogSimRules(catalog, BalrogSetup):
    out = open(BalrogSetup.catruleslog, 'w')
    for key in catalog.galaxyrule.keys():
        out.write('%s %s %s\n' %(key, catalog.galaxyrule[key].type, str(catalog.galaxyrule[key].param)) )

    out.write('\n')
    for i in range(len(catalog.rule)):
        for key in catalog.rule[i].keys():
            out.write('%s %s %s %s\n' %(str(i), key, catalog.rule[i][key].type, str(catalog.rule[i][key].param)) )


def LogDerivedOpts(cmdline_args, BalrogSetup):
    out = open(BalrogSetup.derivedlog, 'w')
    ArgsDict = vars(BalrogSetup)
    VetoDict = vars(cmdline_args)
    for key in ArgsDict:
        if key not in VetoDict:
            out.write('%s %s\n' %(key, ArgsDict[key]) )
    out.close()


def LogExtraSexConfig(ExtraSexConfig, BalrogSetup):
    out = open(BalrogSetup.extrasexlog, 'w')
    veto = NoOverride()
    for key in ExtraSexConfig:
        if key not in veto:
            out.write('%s %s\n' %(key, ExtraSexConfig[key]) )
    out.close()
  
    
def NoOverride():
    keys = ['IMAGE',
            'WEIGHT_IMAGE',
            'CATALOG_NAME',
            'c',
            'PARAMETERS_NAME',
            'STARNNW_NAME',
            'FILTER_NAME',
            'MAG_ZEROPOINT',
            'PSF_NAME',
            'ASSOC_NAME',
            'ASSOC_PARAMS',
            'ASSOC_DATA',
            'ASSOC_RADIUS',
            'ASSOC_TYPE',
            'ASSOCSELEC_TYPE']
    return keys



def LogCmdlineOpts(cmdline_args, cmdline_args_copy, BalrogSetup):
    out = open(BalrogSetup.cmdlinelog, 'w')

    ArgsDict = vars(cmdline_args)
    ordered = CmdlineListOrdered()
    for key in ordered:
        out.write('%s %s\n' %(key, ArgsDict[key]) )

    out.write('\n')
    ArgsDict = vars(cmdline_args_copy)
    for key in ArgsDict.keys():
        if key not in ordered:
            out.write('%s %s\n' %(key, ArgsDict[key]) ) 

    out.close()

   
def CmdlineListOrdered():
    args = ["imagein", "imageext", "weightin", "weightext", "psfin","outdir", "clean",
            "xmin", "xmax", "ymin", "ymax","ngal", "seed", "gain", "zeropoint","fluxthresh", 
            "sexpath", "sexconfig", "sexparam", "sexnnw", "sexconv", "noempty", "sexemptyparam", "noassoc"]
    return args


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
        args.xmax = pyfits.open(args.imagein)[args.imageext].header['NAXIS1']
    if args.ymax==-1:
        args.ymax = pyfits.open(args.imagein)[args.imageext].header['NAXIS2']

    try:
        args.gain = float(args.gain)
    except:
        try:
            args.gain = pyfits.open(args.imagein)[args.imageext].header[args.gain]
        except:
            args.gain = 1.0

    try:
        args.zeropoint = float(args.zeropoint)
    except:
        try:
            args.zeropoint = pyfits.open(args.imagein)[args.imageext].header[args.zeropoint]
        except:
            args.zeropoint = 30.0

    return args



def DefaultArgs(parser):
    # Input and (temporary) output Images
    parser.add_argument( "-od", "--outdir", help="Directory where to put output. Output files will be named based on the input file name.", default=None)
    parser.add_argument( "-ii", "--imagein", help="Input image to draw simulated galaxies into.", type=str, default=None)
    parser.add_argument( "-ie", "--imageext", help="FITS extension where the image lives in the input file.", type=int, default=0)
    parser.add_argument( "-wi", "--weightin", help="Weight map of input image.", type=str, default=None)
    parser.add_argument( "-we", "--weightext", help="FITS extension where the weight map lives in the input weight file.", type=int, default=None)
    parser.add_argument( "-pi", "--psfin", help="PSF of thin input image, to be convolved with simulated galaxies.", type=str, default=None)

    # Properties you want your simulated image to have
    parser.add_argument( "-xmin", "--xmin", help="Minimum column of extracted subimage, indexing ranges from (1,NumPixelsX).", type=int, default=1)
    parser.add_argument( "-xmax", "--xmax", help="Maximum column of extracted subimage, indexing ranges from (1,NumPixelsX).", type=int, default=-1)
    parser.add_argument( "-ymin", "--ymin", help="Minimum row of extracted subimage, indexing ranges from (1,NumPixelsY).", type=int, default=1)
    parser.add_argument( "-ymax", "--ymax", help="Maximum row of extracted subimage, indexing ranges from (1,NumPixelsY).", type=int, default=-1)
    parser.add_argument( "-ngal", "--ngal", help="Number of simulated galaxies", type=int, default=50)
    parser.add_argument( "-gain", "--gain", help="Gain, needed for adding noise. Can be a float or a keyword from the image header. (Default reads image header keyword 'GAIN'. If that fails, default is set to 1)", default='GAIN')
    parser.add_argument( "-zp", "--zeropoint", help="Zeropoint used to convert simulated magnitude to flux. Sextractor runs with this zeropoint. Can be a float or a keyword from the image header. (Default looks for keyword 'SEXMGZPT'. If given keyword is not found, zeropoint defaults to 30.)", default='SEXMGZPT')
    parser.add_argument( "-s", "--seed", help="Seed for random number generation when simulating galaxies. This does not apply to noise realizations, which are always random.", type=int, default=None)
    parser.add_argument( "-ft", "--fluxthresh", help="Flux value where to cutoff the postage stamp", type=float, default=0.01)
    parser.add_argument( "-c", "--clean", help="Delete output image files", action="store_true")

    # How to run sextractor
    parser.add_argument( "-spp", "--sexpath", help='Path for sextractor binary', type=str, default='sex')
    parser.add_argument( "-sc", "--sexconfig", help='Sextractor config file', type=str, default=None)
    parser.add_argument( "-sp", "--sexparam", help='Sextractor param file', type=str, default=None)
    parser.add_argument( "-sn", "--sexnnw", help='Sextractor neural network S/G file', type=str, default=None)
    parser.add_argument( "-sv", "--sexconv", help='Sextractor filter convolution file', type=str, default=None)
    parser.add_argument( "-na", "--noassoc", help="Don't do association mode matching in sextractor. Association mode is sextractor lingo for only look for sources at certain positions; in Balrog, the simulated galaxy positions. Using association mode is much faster.", action="store_true")
    parser.add_argument( "-ne", "--noempty", help="Skip sextractor run over original image, prior to any simulation. One usage for such a run is to identify cases where a galaxy is simulated in the same position as something originally there. Depending on how the objects' properties conspire, Sextractor may not know any blending happened, ", action="store_true")
    parser.add_argument( "-sep", "--sexemptyparam", help="Sextractor param file for run over original image, prior to any simulation, If only interested in the run for 'deblending' issues, the file's contents are mostly irrelevant. The default file does not do model fitting to be faster.", type=str, default=None)




if __name__ == "__main__":

    # Parse command line options and user configurations to setup Balrog.
    cmdline_opts = GetOpts()
    BalrogSetup, simulation_rules, extra_sex_config = ConfigureBalrog(cmdline_opts)


    # Take the the user's configurations and build the simulated truth catalog out of them.
    catalog = GetSimulatedGalaxies(BalrogSetup, simulation_rules)
   

    # Get the subsampled flux and weightmap images, along with the PSF model and WCS.
    bigImage, subWeight, psfmodel, wcs = ReadImages(BalrogSetup)


    # If desired, run sextractor over the image prior to inserting any simulated galaxies.
    if not BalrogSetup.noempty:
        NosimRunSextractor(BalrogSetup, bigImage, subWeight, extra_sex_config, catalog)


    # Insert simulated galaxies.
    bigImage = InsertSimulatedGalaxies(bigImage, catalog, psfmodel, BalrogSetup, wcs)
    WriteImages(BalrogSetup, bigImage, subWeight)
    WriteCatalog(catalog, BalrogSetup, txt=None, fits=True)


    # Run sextractor over the simulated image.
    RunSextractor(BalrogSetup, extra_sex_config, catalog)


    # If chosen, clean up image files you don't need anymore
    if BalrogSetup.clean:
        Cleanup(BalrogSetup)


    # Log some  extra stuff Balrog used along the way
    LogDerivedOpts(cmdline_opts, BalrogSetup)
