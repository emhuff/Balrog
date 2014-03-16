#!/usr/bin/env python

import copy
import datetime
import os
import sys
import subprocess
import argparse
import logging
import traceback
import numpy as np
import astropy.io.fits as pyfits
import galsim
import galsim.des
import sextractor_engine
from model_class import *
from config import *
from balrogexcept import *



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
        if os.path.lexists(BalrogSetup.catalogtruth):
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
    if os.path.lexists(psfout):
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
        if os.path.lexists(file):
            subprocess.call(['rm',file])


def UserDefinitions(cmdline_args, BalrogSetup):
    rules = SimRules(BalrogSetup.ngal)
    ExtraSexConfig = {}
    results = Results()
    cmdline_args_copy = copy.copy(cmdline_args)

    CustomParseArgs(cmdline_args_copy)
    SimulationRules(cmdline_args_copy,rules,results)
    SextractorConfigs(cmdline_args_copy, ExtraSexConfig)

    LogCmdlineOpts(cmdline_args, cmdline_args_copy, BalrogSetup)
    LogExtraSexConfig(ExtraSexConfig, BalrogSetup)
    return rules, ExtraSexConfig


def GetSimulatedGalaxies(BalrogSetup, simgals):
    simgals.Sample(BalrogSetup)
    LogSimRules(simgals, BalrogSetup)
    return simgals


def CompError(name, i):
    if name=='flux':
        name = 'magnitude'
    raise RulesAssignmentError(303, 'component %i of %s' %(i, name))

def GalError(name):
    raise RulesAssignmentError(303, name)


def InitializeSersic(rules, sampled, nProfiles=1):
    rules.InitializeSersic(nProfiles=nProfiles)
    sampled.InitializeSersic(nProfiles=nProfiles)


class CompRules(object):
    def __init__(self, nProfiles, name):
        super(CompRules, self).__setattr__('rules', [None]*nProfiles)
        super(CompRules, self).__setattr__('name', name)
        super(CompRules, self).__setattr__('nProfiles', nProfiles)

    def __setattr__(self, name, value):
        raise RulesComponentAttributeError(305)

    def __getattr__(self, name):
        raise RulesComponentAttributeError(305)
  
    def __len__(self):
        return self.nProfiles

    def __getitem__(self, index):
        if index >= self.nProfiles:
            raise RulesIndexOutOfRange(304, self.name, self.nProfiles)
        return self.rules[index]

    def __setitem__(self, index, value):
        if index >= self.nProfiles:
            raise RulesIndexOutOfRange(304, self.name, self.nProfiles)
        rule = self._CheckRule(value, index)
        self.rules[index] = rule


    def _CheckRule(self, rule, i):
        if type(rule).__name__!='Rule':
            if rule==None:
                pass
            elif type(rule)==float or type(rule)==int:
                rule = Value(float(rule))
            else:
                try:
                    arr = np.array(rule)
                    if arr.ndim==1 and arr.size==self.ngal:
                        rule = Array(arr)
                    else:
                        CompError(self.name, i)
                except:
                    CompError(self.name, i)
        return rule


class SimRules(object):
    def __init__(self, ngal):
        super(SimRules, self).__setattr__('ngal', ngal)
        for name in self._GetGalaxy():
            super(SimRules, self).__setattr__(name, None)
        self.InitializeSersic()

    def InitializeSersic(self, nProfiles=1):
        super(SimRules, self).__setattr__('nProfiles', nProfiles)
        for c in self._GetComponent():
            super(SimRules, self).__setattr__(c, CompRules(nProfiles,c))

    def _GetGalaxy(self):
        return ['x','y','g1','g2','magnification']

    def _GetComponent(self):
        return ['axisratio','beta','halflightradius','magnitude','sersicindex']

    def __getitem__(self, index):
        raise RulesIndexingError(302)

    def __setitem__(self, index, value):
        raise RulesIndexingError(302)

    def __getattr__(self, name):
        raise RulesAttributeError(301, name)

    def __setattr__(self, name, value):
        if name=='ngal':
            raise RulesNgalError(-1)

        elif name=='nProfiles':
            raise RulesnProfilesError(-2, name)

        elif name in ['x', 'y', 'g1', 'g2', 'magnification']:
            value = self._CheckRule(name, value, 'galaxy')
            super(SimRules, self).__setattr__(name, value)

        elif name in ['axisratio', 'beta', 'halflightradius', 'magnitude', 'sersicindex']:
            try:
                size = len(value)
            except:
                raise RulesAssignmentNoArrayError(306)
            if size!=self.nProfiles:
                raise RulesAssignmentNoArrayError(306)

            for i in range(size):
                value[i] = self._CheckRule(name, value[i], 'component', i=i)
                exec "self.%s[%i] = value[%i]" %(name, i, i)

        else:
            raise RulesAttributeError(301,name)



    def _CheckRule(self, name, rule, kind, i=None):
        if type(rule).__name__!='Rule':
            if rule==None:
                pass
            elif type(rule)==float or type(rule)==int:
                rule = Value(float(rule))
            else:
                try:
                    arr = np.array(rule)
                    if arr.ndim==1 and arr.size==self.ngal:
                        rule = Array(arr)
                    else:
                        if kind=='galaxy':
                            GalError(name)
                        else:
                            CompError(name, i)
                except:
                    if kind=='galaxy':
                        GalError(name)
                    else:
                        CompError(name, i)
        return rule


class CompResult(object):
    def __init__(self, nProfiles, name):
        super(CompResult, self).__setattr__('name', name)
        super(CompResult, self).__setattr__('nProfiles', nProfiles)

    def __getitem__(self,index):
        if index >= self.nProfiles:
            raise SampledIndexOutOfRange(404, self.name, self.nProfiles)

        if self.name=='magnitude':
            return Same( (index,'flux') )
        else:
            return Same( (index,self.name) )

    def __setitem__(self, index, value):
        raise SampledAssignmentError(403, '%s[%i]'%(self.name,index))

    def __setattr__(self, name, value):
        raise SampledAssignmentError(403, '%s.%s'%(self.name,name))

    def __getattr__(self, name):
        raise SampledComponentAttributeError(405)


class Results(object):
    def __init__(self):
        self.InitializeSersic()

    def InitializeSersic(self, nProfiles=1):
        super(Results, self).__setattr__('nProfiles', nProfiles)
        for c in self._GetComponent():
            super(Results, self).__setattr__(c, CompResult(nProfiles,c))

    def _GetGalaxy(self):
        return ['x','y','g1','g2','magnification']

    def _GetComponent(self):
        return ['axisratio','beta','halflightradius','magnitude','sersicindex']

    def __getattr__(self, name):
        if name not in self._GetGalaxy():
            raise SampledAttributeError(401, name)
        else:
            return Same(name)
    
    def __getitem__(self, index):
        raise SampledIndexingError(402)

    def __setitem__(self, index, value):
        raise SampledIndexingError(402)

    def __setattr__(self, name, value):
        raise SampledAssignmentError(403, name)
        
    

class DerivedArgs():
    def __init__(self,args,log):
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
        self.balroglog = DefaultName(args.imagein, '.fits', '.balrog.log.txt', self.logdir)

        CreateSubDir(self.imgdir)
        CreateSubDir(self.catdir)
        CreateSubDir(self.sexdir)

        thisdir = os.path.dirname( os.path.realpath(__file__) )
        config = os.path.join(thisdir, 'config.py')
        self.CopyFile(args.sexconfig, self.sexdir)
        self.CopyFile(args.sexparam, self.sexdir)
        self.CopyFile(args.sexnnw, self.sexdir)
        self.CopyFile(args.sexconv, self.sexdir)
        self.CopyFile(args.sexemptyparam, self.sexdir)
        self.CopyFile(config, self.logdir)

        self.outimageext = 0
        self.outweightext = 0
        if self.weightout==self.imageout:
            self.outweightext = self.outimageext + 1

        self.subsample = True
        if args.xmin==1 and args.ymin==1 and args.xmax==pyfits.open(args.imagein)[args.imageext].header['NAXIS1'] and args.ymax==pyfits.open(args.imagein)[args.imageext].header['NAXIS2']:
            self.subsample = False
         
        self.logger = log


    def CopyFile(self, file, dir):
        basename = os.path.basename(file)
        outname = os.path.join(dir,basename)
        if os.path.exists(outname):
            subprocess.call(['rm', outname])
        subprocess.call(['cp', file, outname])

        
class BalrogConfig():
    def __init__(self, cargs, dargs):
        cdict = vars(cargs)
        for key in cdict.keys():
            exec "self.%s = cdict['%s']" %(key, key)

        ddict = vars(dargs)
        for key in ddict.keys():
            exec "self.%s = ddict['%s']" %(key, key)


def SetLevel(h, choice):
    if choice=='q':
        h.setLevel(logging.ERROR)
    elif choice=='n':
        h.setLevel(logging.WARNING)
    elif choice=='v':
        h.setLevel(logging.INFO)
    elif choice=='vv':
        h.setLevel(logging.DEBUG)
    return h


def SetupLogger(cmdline_args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    outdir = os.path.join(defdir, 'output')
    if cmdline_args.outdir==None:
        cmdline_args.outdir = outdir

    CreateDir(cmdline_args.outdir)
    logdir = os.path.join(cmdline_args.outdir, 'balrog_log')
    CreateSubDir(logdir)
    
    log = logging.getLogger()
    log.setLevel(logging.NOTSET)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    file = os.path.join(logdir, 'run.log.txt')
    fh = logging.FileHandler(file, mode='w')
    fh.setFormatter(formatter)
    SetLevel(fh, cmdline_args.logverbosity)
    log.addHandler(fh)

    ph = logging.StreamHandler(stream=sys.stderr)
    ph.setFormatter(formatter)
    ph = SetLevel(ph, cmdline_args.stdverbosity)
    log.addHandler(ph)
    return log


def GetOpts():
    parser = argparse.ArgumentParser()
    DefaultArgs(parser)
    CustomArgs(parser) 
    cmdline_args = parser.parse_args()
    log = SetupLogger(cmdline_args)
    ParseDefaultArgs(cmdline_args,log)
    return cmdline_args, log


def ConfigureBalrog(cmdline_opts, log):
    cmdline_opts, derived_opts = GetMoreOpts(cmdline_opts, log)
    BalrogSetup = BalrogConfig(cmdline_opts, derived_opts)
    return BalrogSetup


def GetMoreOpts(cmdline_args, log):
    derived_args = DerivedArgs(cmdline_args, log)
    return [cmdline_args, derived_args]


def LogSimRules(catalog, BalrogSetup):
    out = open(BalrogSetup.catruleslog, 'w')
    for key in catalog.galaxyrule.keys():
        out.write('%s %s %s\n' %(key, catalog.galaxyrule[key].type, str(catalog.galaxyrule[key].param)) )

    out.write('\n')
    for i in range(len(catalog.componentrule)):
        for key in catalog.componentrule[i].keys():
            out.write('%s %s %s %s\n' %(str(i), key, catalog.componentrule[i][key].type, str(catalog.componentrule[i][key].param)) )


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
    err = 0
    if not os.path.lexists(subdir):
        err = subprocess.call(['mkdir', subdir], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    if err!=0:
        raise SubdirWriteError(202, subdir)


def CreateDir(dir):
    gdir = copy.copy(dir)
    full = False
    while dir[0]=='/':
        dir = dir[1:]
        full = True
    while dir[-1]=='/':
        dir = dir[:-1]
    dirs = dir.strip().split('/')
    if full:
        subdir = '/'
    else:
        subdir = './'

    err = 0
    for dir in dirs:
        subdir = os.path.join(subdir,dir)
        if not os.path.lexists(subdir):
            err = subprocess.call( ['mkdir', subdir], stderr=subprocess.PIPE, stdout=subprocess.PIPE )
            if err!=0:
                raise OutdirWriteError(201, gdir, subdir)


def ImagesExistenceCheck(args, log):
    if not os.path.lexists(args.imagein):
        raise ImageInputError(101, 'image', 'imagein', args.imagein)
    if not os.path.lexists(args.weightin):
        raise ImageInputError(102, 'weight', 'weightin', args.weightin)
    if not os.path.lexists(args.psfin):
        raise ImageInputError(103, 'PSF', 'psfin', args.psfin)


def ImagesAreFITS(args, log):
    try:
        pyfits.open(args.imagein)
    except:
        raise FitsFileError(111, 'image', 'imagein', args.imagein) 
    try:
        pyfits.open(args.weightin)
    except:
        raise FitsFileError(112, 'weight', 'weightin', args.weightin)
    try:
        pyfits.open(args.psfin)
    except:
        raise FitsFileError(113, 'PSF', 'psfin', args.psfin)


def ExtsExist(args, log):
    try: 
        hdu = pyfits.open(args.imagein)[args.imageext]
    except:
        raise FitsExtError(121, 'image', 'imageext', args.imageext, 'imagein', args.imagein)
    try:
        hdu = pyfits.open(args.weightin)[args.weightext]
    except:
        raise FitsExtError(122, 'weight', 'weightext', args.weightext, 'weightin', args.weightin)


def SizesOK(args, log):
    try:
        ix = pyfits.open(args.imagein)[args.imageext].header['NAXIS1']
    except:
        raise FitsHeaderError(131,'image','NAXIS1', 'imagein',args.imagein, 'imageext',args.imageext)
    try:
        iy = pyfits.open(args.imagein)[args.imageext].header['NAXIS2']
    except:
        raise FitsHeaderError(132,'image','NAXIS2', 'imagein',args.imagein, 'imageext',args.imageext)
    try:
        wx = pyfits.open(args.weightin)[args.weightext].header['NAXIS1']
    except:
        raise FitsHeaderError(133,'weight','NAXIS1', 'weightin',args.weightin, 'weightext',args.weightext)
    try:
        wy = pyfits.open(args.weightin)[args.weightext].header['NAXIS2']
    except:
        raise FitsHeaderError(134,'weight','NAXIS2', 'weightin',args.weightin, 'weightext',args.weightext)

    if ix!=wx or iy!=wy:
        raise SizeMismatchError(135, ix,iy,wx,wy)

    if args.xmax==-1:
        args.xmax = ix
    if args.ymax==-1:
        args.ymax = iy

    if args.xmin > args.xmax:
        raise SizeError(136, 'x', args.xmin, args.xmax)
    if args.ymin > args.ymax:
        raise SizeError(137, 'y', args.ymin, args.ymax)

    if args.xmin < 1:
        log.warning('Given --xmin = %i is less than 1. Setting --xmin = 1' %(args.xmin))
        args.xmin = 1
    if args.ymin < 1:
        log.warning('Given --ymin = %i is less than 1. Setting --ymin = 1' %(args.ymin))
        args.ymin = 1

    if args.xmax > ix:
        log.warning('Given --xmax = %i exceeds the image bounds. Setting --xmax to the right edge of image: --xmax = %i' %(args.xmax,ix))
        args.xmax = ix
    if args.ymax > iy:
        log.warning('Given --ymax = %i exceeds the image bounds. Setting --ymax to the upper edge of image: --ymax = %i' %(args.ymax,iy))
        args.ymax = iy



def FindImages(args, log, indir):
    default = os.path.join(indir, 'example.fits')
    if args.imagein!=None and os.path.abspath(args.imagein)!=default and args.psfin==None:
        raise PsfInputError(104, args.imagein)

    if args.imagein==None:
        args.imagein = default
        log.info('No --imagein explicitly given. Will use the default. Setting --imagein = %s' %args.imagein)
    if args.weightin==None:
        log.info('No --weightin explicitly given. Assuming it lives in same file as the image. Setting --weightin = %s' %args.imagein)
        args.weightin = args.imagein
    if args.weightext==None and os.path.abspath(args.weightin)==args.imagein:
        args.weightin = os.path.abspath(args.weightin)
        log.info('No --weightext explicitly given. Assuming it lives in the extension after the image. Setting --weightext = %s' %(args.imageext+1))
        args.weightext = args.imageext + 1
    if args.weightext==None:
        log.info('No --weightext explicitly given. Assuming it lives in the 0th extension. Setting --weightext = 0')
        args.weightext = 0
    if args.psfin==None:
        args.psfin = os.path.join(indir, 'example.psf')
        log.info('No --psfin explicitly given. Will use the default. Setting --psfin = %s' %args.psfin)
    

def ParseImages(args, log, indir):
    FindImages(args, log, indir)
    ImagesExistenceCheck(args, log)
    ImagesAreFITS(args, log)
    ExtsExist(args, log)
    SizesOK(args, log)


def ParseFloatKeyword(args, log):
    try:
        args.gain = float(args.gain)
        log.debug('--gain given as a float value: %f' %(args.gain))
    except:
        try:
            g = args.gain
            args.gain = pyfits.open(args.imagein)[args.imageext].header[args.gain]
            log.debug('--gain given as image header keyword: %s = %f' %(g,args.gain))
        except:
            if args.gain=='GAIN':
                log.info('Did not find keyword GAIN in image header. Setting --gain = 1.0')
            else:
                log.warning('Did not find keyword %s in image header. Setting --gain = 1.0' %(args.gain))
            args.gain = 1.0

    try:
        args.zeropoint = float(args.zeropoint)
        log.debug('--zeropoint given as a float value: %f' %(args.zeropoint))
    except:
        try:
            z = args.zeropoint
            args.zeropoint = pyfits.open(args.imagein)[args.imageext].header[args.zeropoint]
            log.debug('--zeropoint given as image header keyword: %s = %f' %(z,args.zeropoint))
        except:
            if args.zeropoint=='SEXMGZPT':
                log.info('Did not find keyword SEXMGZPT in image header. Setting --zeropoint = 30.0')
            else:
                log.warning('Did not find keyword %s in image header. Setting --zeropoint = 30.0' %(args.zeropoint))
            args.zeropoint = 30.0


def FindSexFile(arg, log, configdir, default, label):
    default = os.path.join(configdir, default)
    if arg==None:
        arg = default
        log.info('No --%s explicitly given. Will use the default. Setting --%s = %s' %(label,label,arg))
    if not os.path.lexists(arg):
        log.warning('--%s %s does not exist. Assuming a default file. Setting --%s = %s' %(label, arg, label, default))
        arg = default
    return arg


def ParseSex(args, log, configdir):
    args.sexconfig = FindSexFile(args.sexconfig, log, configdir, 'sex.config', 'sexconfig')
    args.sexparam = FindSexFile(args.sexparam, log, configdir, 'bulge.param', 'sexparam')
    args.sexemptyparam = FindSexFile(args.sexemptyparam, log, configdir, 'sex.param', 'sexemptyparam')
    args.sexnnw = FindSexFile(args.sexnnw, log, configdir, 'sex.nnw', 'sexnnw')
    args.sexconv = FindSexFile(args.sexconv, log, configdir, 'sex.conv', 'sexconv')


def ParseDefaultArgs(args,log):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    indir = os.path.join(defdir, 'input')
    configdir = os.path.join(thisdir, 'astro_config')
    
    ParseImages(args, log, indir)
    ParseFloatKeyword(args, log) 
    ParseSex(args, log, configdir)

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

    # Other Balrog stuff
    parser.add_argument( "-c", "--clean", help="Delete output image files", action="store_true")
    parser.add_argument( "-sv", "--stdverbosity", help="Verbosity level of stdout", type=str, default='n', choices=['n','v','vv','q'])
    parser.add_argument( "-lv", "--logverbosity", help="Verbosity level of log file", type=str, default='n', choices=['n','v','vv'])

    # How to run sextractor
    parser.add_argument( "-spp", "--sexpath", help='Path for sextractor binary', type=str, default='sex')
    parser.add_argument( "-sc", "--sexconfig", help='Sextractor config file', type=str, default=None)
    parser.add_argument( "-sp", "--sexparam", help='Sextractor param file', type=str, default=None)
    parser.add_argument( "-sn", "--sexnnw", help='Sextractor neural network S/G file', type=str, default=None)
    parser.add_argument( "-sf", "--sexconv", help='Sextractor filter convolution file', type=str, default=None)
    parser.add_argument( "-na", "--noassoc", help="Don't do association mode matching in sextractor. Association mode is sextractor lingo for only look for sources at certain positions; in Balrog, the simulated galaxy positions. Using association mode is much faster.", action="store_true")
    parser.add_argument( "-ne", "--noempty", help="Skip sextractor run over original image, prior to any simulation. One usage for such a run is to identify cases where a galaxy is simulated in the same position as something originally there. Depending on how the objects' properties conspire, Sextractor may not know any blending happened, ", action="store_true")
    parser.add_argument( "-sep", "--sexemptyparam", help="Sextractor param file for run over original image, prior to any simulation, If only interested in the run for 'deblending' issues, the file's contents are mostly irrelevant. The default file does not do model fitting to be faster.", type=str, default=None)


def RaiseException():
    exc_info = sys.exc_info()
    config_errs = []
    err_list = traceback.extract_tb(exc_info[2])
    for err in err_list: 
        file = err[0]
        if file.find('config.py')!=-1:
            config_errs.append(err)
        '''
        else:
            config_errs.append(err)
        '''
    keep = traceback.format_list(config_errs)
    keep_tb = ''.join(keep)
    logging.error('Run error caused Balrog to exit.\n%s' %(keep_tb), exc_info=(exc_info[0], exc_info[1], None))
    sys.exit()


def GetNativeOptions():
    parser = argparse.ArgumentParser()
    DefaultArgs(parser)
    return parser


def AddCustomOptions(parser):
    CustomArgs(parser) 


def NativeParse(parser):
    cmdline_opts = parser.parse_args()
    log = SetupLogger(cmdline_opts)
    ParseDefaultArgs(cmdline_opts,log)
    BalrogSetup = ConfigureBalrog(cmdline_opts, log)
    return cmdline_opts, BalrogSetup


def CustomParse(cmdline_opts, BalrogSetup):
    rules, extra_sex_config = UserDefinitions(cmdline_opts, BalrogSetup)
    rules = DefineRules(BalrogSetup, x=rules.x, y=rules.y, g1=rules.g1, g2=rules.g2, magnification=rules.magnification, nProfiles=rules.nProfiles, axisratio=rules.axisratio, beta=rules.beta, halflightradius=rules.halflightradius, magnitude=rules.magnitude, sersicindex=rules.sersicindex)
    return rules, extra_sex_config


def RunBalrog():
    # Get Native argument parser
    parser = GetNativeOptions()

    # Add the user's command line options
    AddCustomOptions(parser)

    # Parse the command line agruments and interpret the user's settings for the simulation
    cmdline_opts, BalrogSetup = NativeParse(parser)
    rules, extra_sex_config = CustomParse(cmdline_opts, BalrogSetup)

    # Take the the user's configurations and build the simulated truth catalog out of them.
    catalog = GetSimulatedGalaxies(BalrogSetup, rules)
   
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



if __name__ == "__main__":
    try:
        RunBalrog()
    except:
        RaiseException()

