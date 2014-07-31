#!/usr/bin/env python

import time
import imp
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

def FindInCat(rule, name):
    cat = rule.param[0]
    ext = rule.param[1]
    n = rule.param[2]

    hdu = pyfits.open(cat)[ext]
    #cut = np.in1d(hdu.columns.names,np.array([name]))
    cut = np.in1d(hdu.columns.names,np.array([n]))

    return cut, hdu

def WriteCatalog(sample, BalrogSetup, txt=None, fits=False, TruthCatExtra=None, extracatalog=None):
    columns = []
    for key in sample.galaxy.keys():
        name = '%s' %(key)
        arr = sample.galaxy[key]
        if key=='x':
            arr = arr - BalrogSetup.xmin + 1
            unit = 'pix'
        elif key=='y':
            arr = arr - BalrogSetup.ymin + 1
            unit = 'pix'
        elif key.find('flux')!=-1:
            unit = 'ADU'
        else:
            unit = 'dimensionless'
        col = pyfits.Column(name=name, array=arr,format='E', unit=unit)
        columns.append(col)

    if TruthCatExtra!=None:
        for rule,name,fmt,unit in zip(TruthCatExtra.rules,TruthCatExtra.names,TruthCatExtra.fmts, TruthCatExtra.units):
            if unit==None:
                if rule.type!='catalog':
                    unit = 'unspecified'
                    BalrogSetup.runlogger.info('No unit was specificed for additional truth catalog column %s' %(name))
                else:
                    cut,hdu = FindInCat(rule, name)
                    unit = np.array(hdu.columns.units)[cut][0]
            if unit.strip()=='':
                unit = 'unspecified'

            if fmt==None:
                if rule.type!='catalog':
                    if len(extracatalog.galaxy[name])==0:
                        #fmt = 'E'
                        raise TableUnknownType(401, name)
                    else:
                        arrtype = str(type(extracatalog.galaxy[name][0]))
                        if arrtype.find('str')!=-1:
                            strlen = len(max(extracatalog.galaxy[name]))
                            fmt = '%iA' %(strlen)
                        elif arrtype.find('int')!=-1:
                            fmt = 'J'
                        else:
                            fmt = 'E'
                else:
                    cut,hdu = FindInCat(rule, name)
                    fmt = np.array(hdu.columns.formats)[cut][0]
                
            col = pyfits.Column(array=extracatalog.galaxy[name], name=name, format=fmt, unit=unit)
            columns.append(col)


    for i in range(len(sample.component)):
        for key in sample.component[i].keys():
            name = '%s_%i' %(key,i)
            if key.find('halflightradius')!=-1:
                col = pyfits.Column(name=name, array=sample.component[i][key]/np.sqrt(sample.component[i]['axisratio']), format='E', unit='arcsec')
            else:
                if key.find('sersicindex')!=-1:
                    unit = 'dimensionless'
                if key.find('flux')!=-1:
                    unit = 'ADU'
                if key.find('beta')!=-1:
                    unit = 'deg'
                if key.find('axisratio')!=-1:
                    unit = 'dimensionless'
                col = pyfits.Column(name=name, array=sample.component[i][key],format='E', unit=unit)
            columns.append(col)

    try:
        tbhdu = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(columns))
    except:
        tbhdu = pyfits.new_table(pyfits.ColDefs(columns))
    tbhdu.header['XSTART'] = BalrogSetup.xmin
    tbhdu.header['XEND'] = BalrogSetup.xmax
    tbhdu.header['YSTART'] = BalrogSetup.ymin
    tbhdu.header['YEND'] = BalrogSetup.ymax
    tbhdu.header['NSIM'] = BalrogSetup.ngal
    tbhdu.header['ZP'] = BalrogSetup.zeropoint

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
    mhead = mhdus[BalrogSetup.catext].header
    for i in range(len(BalrogSetup.assocnames)):
        mhead[ 'V%i'%i ] = BalrogSetup.assocnames[i]
        if BalrogSetup.assocnames[i] in ['x','y']:
            unit = 'pix'
        elif BalrogSetup.assocnames[i].find('halflightradius')!=-1:
            unit = 'arcsec'
        elif BalrogSetup.assocnames[i].find('beta')!=-1:
            unit = 'deg'
        elif BalrogSetup.assocnames[i].find('flux')!=-1:
            unit = 'ADU'
        else:
            unit = 'dimensionless'
        mhead[ 'VUNIT%i'%i ] = unit
    mhdus.close() 



def ReadImages(BalrogSetup):
    image = galsim.fits.read(BalrogSetup.image, hdu=BalrogSetup.imageext)
    weight = galsim.fits.read(BalrogSetup.weight, hdu=BalrogSetup.weightext)
    if image.wcs==galsim.PixelScale(1):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'fiducialwcs.fits')
        image.wcs = galsim.GSFitsWCS(file)
        weight.wcs = image.wcs
        BalrogSetup.wcshead = file
        BalrogSetup.runlogger.warning('No WCS was found in the input image header. Using default pixel scale of 0.263 arcsec/pixel.')
    wcs = image.wcs

    subBounds = galsim.BoundsI(BalrogSetup.xmin,BalrogSetup.xmax,BalrogSetup.ymin,BalrogSetup.ymax)
    image = image[subBounds]
    weight = weight[subBounds]
    psfmodel = galsim.des.DES_PSFEx(BalrogSetup.psf, BalrogSetup.wcshead)

    return image, weight, psfmodel, wcs


def WriteImages(BalrogSetup, image, weight, nosim=False):
    weightout = BalrogSetup.weightout
    if nosim:
        imageout = BalrogSetup.nosim_imageout
    else:
        imageout = BalrogSetup.imageout

    if BalrogSetup.weight==BalrogSetup.image:
        galsim.fits.writeMulti(image_list=[image,weight], file_name=imageout)
    else:
        galsim.fits.write(image=image, file_name=imageout)
        if BalrogSetup.nonosim or nosim:
            galsim.fits.write(image=weight, file_name=weightout)

    if not BalrogSetup.psf_written:
        WritePsf(BalrogSetup, BalrogSetup.psf, BalrogSetup.psfout)
        BalrogSetup.psf_written = True


def WritePsf(BalrogSetup, psfin, psfout):
    psfhdus = pyfits.open(psfin)
    psfhdus[1].header['POLZERO1'] = psfhdus[1].header['POLZERO1'] - (BalrogSetup.xmin - 1)
    psfhdus[1].header['POLZERO2'] = psfhdus[1].header['POLZERO2'] - (BalrogSetup.ymin - 1)
    if os.path.lexists(psfout):
        subprocess.call(['rm', psfout])
    psfhdus.writeto(psfout)



def InsertSimulatedGalaxies(bigImage, simulatedgals, psfmodel, BalrogSetup, wcs, gspcatalog):
    #psizes, athresh = simulatedgals.GetPSizes(BalrogSetup, wcs)

    t0 = datetime.datetime.now()
    rt = long( t0.microsecond )

    '''
    simulatedgals.galaxy['flux_noised'] = np.copy(simulatedgals.galaxy['x'])
    simulatedgals.galaxy['flux_noiseless'] = np.copy(simulatedgals.galaxy['x'])
    '''
    simulatedgals.galaxy['flux_noised'] = np.zeros( len(simulatedgals.galaxy['x']) )
    simulatedgals.galaxy['flux_noiseless'] = np.zeros( len(simulatedgals.galaxy['x']) )
    simulatedgals.galaxy['not_drawn'] = np.array( [0]*len(simulatedgals.galaxy['x']) )

    for i in range(BalrogSetup.ngal):
        start = datetime.datetime.now()

        #postageStampSize = int(psizes[i])
        d = {}
        for key in gspcatalog.galaxy.keys():
            if gspcatalog.galaxy[key]!=None:
                if key in ['minimum_fft_size', 'maximum_fft_size', 'range_for_extrema']:
                    mod = gspcatalog.galaxy[key][i] % 1
                    if mod > 0.0:
                        d[key] = gspcatalog.galaxy[key][i]
                    else:
                        d[key] = int(gspcatalog.galaxy[key][i])
                else:
                    d[key] = gspcatalog.galaxy[key][i]
        gsparams = galsim.GSParams(**d)

        try:
            combinedObjConv = simulatedgals.GetConvolved(psfmodel, i, wcs, gsparams, BalrogSetup)
        except:
            simulatedgals.galaxy['not_drawn'][i] = 1
            print simulatedgals.component[0]['sersicindex'][i],simulatedgals.component[0]['halflightradius'][i],simulatedgals.component[0]['flux'][i],simulatedgals.component[0]['axisratio'][i],simulatedgals.component[0]['beta'][i]; sys.stdout.flush()
            continue

        ix = int(simulatedgals.galaxy['x'][i])
        iy = int(simulatedgals.galaxy['y'][i])
   
        '''
        smallImage = galsim.Image(postageStampSize,postageStampSize)
        smallImage.setCenter(ix,iy)
        smallImage.wcs = bigImage.wcs
        smallImage = combinedObjConv.draw(image=smallImage)
        '''
       
        pos = galsim.PositionD(simulatedgals.galaxy['x'][i], simulatedgals.galaxy['y'][i])
        local = wcs.local(image_pos=pos)
        localscale = np.sqrt(local.dudx * local.dvdy)
        #smallImage = combinedObjConv.draw(scale=localscale)

        try:
            smallImage = combinedObjConv.draw(scale=localscale, use_true_center=False)
        except:
            simulatedgals.galaxy['not_drawn'][i] = 1
            print simulatedgals.component[0]['sersicindex'][i],simulatedgals.component[0]['halflightradius'][i],simulatedgals.component[0]['flux'][i],simulatedgals.component[0]['axisratio'][i],simulatedgals.component[0]['beta'][i]; sys.stdout.flush()
            continue

        smallImage.setCenter(ix,iy)

        t1 = datetime.datetime.now()
        dt = t1 - t0
        micro = long( (dt.days*24*60*60 + dt.seconds)*1.0e6 + dt.microseconds ) + rt + i

        bounds = smallImage.bounds & bigImage.bounds
        simulatedgals.galaxy['flux_noiseless'][i] = smallImage.added_flux 
        smallImage.addNoise(galsim.CCDNoise(gain=BalrogSetup.gain,read_noise=0,rng=galsim.BaseDeviate(micro)))
        flux_noised = np.sum(smallImage.array.flatten())
        simulatedgals.galaxy['flux_noised'][i] = flux_noised
        bounds = smallImage.bounds & bigImage.bounds
        bigImage[bounds] += smallImage[bounds]
        
        end = datetime.datetime.now()
        #print (end - start).total_seconds(), simulatedgals.component[0]['sersicindex'][i], simulatedgals.component[0]['halflightradius'][i], simulatedgals.component[0]['axisratio'][i], simulatedgals.component[0]['flux'][i]; sys.stdout.flush()

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
        param_file = BalrogSetup.nosimsexparam

    pfile = DefaultName(catalogmeasured, '.fits', '.sex.params', BalrogSetup.sexdir)
    txt = ParamTxtWithoutAssoc(param_file)
    if not BalrogSetup.noassoc: #and not nosim:
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


def AutoConfig(BalrogSetup, imageout, weightout, catalogmeasured, config_file, param_file, afile, eng, nosim):
    eng.Path(BalrogSetup.sexpath)
    eng.config['IMAGE'] = '%s[%i],%s[%s]' %(imageout,BalrogSetup.outimageext,imageout,BalrogSetup.outimageext)
    eng.config['WEIGHT_IMAGE'] = '%s[%i],%s[%i]' %(weightout,BalrogSetup.outweightext,weightout,BalrogSetup.outweightext)
    eng.config['CATALOG_NAME'] = catalogmeasured
    eng.config['c'] = config_file
    eng.config['PARAMETERS_NAME'] = param_file
    eng.config['STARNNW_NAME'] = BalrogSetup.sexnnw
    eng.config['FILTER_NAME'] = BalrogSetup.sexconv
    eng.config['MAG_ZEROPOINT'] = BalrogSetup.zeropoint
    eng.config['PSF_NAME'] = '%s,%s' %(BalrogSetup.psfout,BalrogSetup.psfout)
    eng.config['CATALOG_TYPE'] = '%s' %(BalrogSetup.catfitstype)

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
        eng.config['ASSOC_PARAMS'] = '%i,%i' %(x,y)

        eng.config['ASSOC_DATA'] = ','.join(inds)
        '''
        if not nosim:
            eng.config['ASSOC_DATA'] = ','.join(inds)
        else:
            eng.config['ASSOC_DATA'] = '%i,%i' %(x,y)
        '''

        eng.config['ASSOC_RADIUS'] = '2.0'
        eng.config['ASSOC_TYPE'] = 'NEAREST'
        eng.config['ASSOCSELEC_TYPE'] = 'MATCHED'
   

def RunSextractor(BalrogSetup, ExtraSexConfig, catalog, nosim=False):
    afile = None
    weightout = BalrogSetup.weightout
    if nosim:
        catalogmeasured = BalrogSetup.nosim_catalogmeasured
        imageout = BalrogSetup.nosim_imageout
        if not BalrogSetup.noassoc:
            afile = BalrogSetup.assoc_nosimfile
    else:
        catalogmeasured = BalrogSetup.catalogmeasured
        imageout = BalrogSetup.imageout
        if not BalrogSetup.noassoc:
            afile = BalrogSetup.assoc_simfile
        if BalrogSetup.image==BalrogSetup.weight:
            weightout = imageout

    if not BalrogSetup.noassoc:
        WriteCatalog(catalog, BalrogSetup, txt=afile, fits=False)

    param_file = WriteParamFile(BalrogSetup, catalogmeasured, nosim)
    config_file = BalrogSetup.sexconfig
    #if not BalrogSetup.noassoc:
    config_file = WriteConfigFile(BalrogSetup, config_file, catalogmeasured)

    eng = sextractor_engine.SextractorEngine()
    for key in ExtraSexConfig.keys():
        eng.config[key] = ExtraSexConfig[key]

    AutoConfig(BalrogSetup, imageout, weightout, catalogmeasured, config_file, param_file, afile, eng, nosim)
    if nosim:
        msg = '# Running sextractor prior to inserting simulated galaxies\n'
    else:
        msg = '# Running sextractor after inserting simulated galaxies\n'
    BalrogSetup.sexlogger.info(msg)
    eng.run(logger=BalrogSetup.sexlogger)

    if not BalrogSetup.noassoc: #and not nosim:
        CopyAssoc(BalrogSetup, catalogmeasured)


def NosimRunSextractor(BalrogSetup, bigImage, subweight, ExtraSexConfig, catalog):
    if BalrogSetup.subsample:
        WriteImages(BalrogSetup, bigImage, subweight, nosim=True)
    else:
        if os.path.lexists(BalrogSetup.nosim_imageout):
            subprocess.call( ['rm', BalrogSetup.nosim_imageout] )
        if os.path.lexists(BalrogSetup.weightout):
            subprocess.call( ['rm', BalrogSetup.weightout] )
        if os.path.lexists(BalrogSetup.psfout):
            subprocess.call( ['rm', BalrogSetup.psfout] )

        subprocess.call( ['ln', '-s', BalrogSetup.image, BalrogSetup.nosim_imageout] )
        subprocess.call( ['ln', '-s', BalrogSetup.psf, BalrogSetup.psfout] )
        BalrogSetup.psf_written = True
        
        if BalrogSetup.weight!=BalrogSetup.image:
            subprocess.call( ['ln', '-s', BalrogSetup.weight, BalrogSetup.weightout] )

    RunSextractor(BalrogSetup, ExtraSexConfig, catalog, nosim=True)


def Cleanup(BalrogSetup):
    files = [BalrogSetup.imageout, BalrogSetup.psfout, BalrogSetup.weightout, BalrogSetup.nosim_imageout]
    for file in files:
        if os.path.lexists(file):
            subprocess.call(['rm',file])


def UserDefinitions(cmdline_args, BalrogSetup, config, galkeys, compkeys):
    rules = SimRules(BalrogSetup.ngal, galkeys, compkeys)
    ExtraSexConfig = {}
    results = Results(galkeys, compkeys)
    cmdline_args_copy = copy.copy(cmdline_args)
    TruthCatExtra = TableColumns(BalrogSetup.ngal)

    if config!=None:
        if 'CustomParseArgs' not in dir(config):
            BalrogSetup.runlogger.warning('The function CustomParseArgs was not found in your Balrog python config file: %s. Will continue without parsing any custom command line arguments.' %BalrogSetup.pyconfig)
        else:
            config.CustomParseArgs(cmdline_args_copy)
        copy1 = copy.copy(cmdline_args_copy)
        copy2 = copy.copy(cmdline_args_copy)

        if 'SimulationRules' not in dir(config):
            BalrogSetup.runlogger.warning('The function SimulationRules was not found in your Balrog python config file: %s. All properties of the simulated galaxies will assume their defaults.' %BalrogSetup.pyconfig)
        else:
            config.SimulationRules(copy1,rules,results, TruthCatExtra)

        if 'SextractorConfigs' not in dir(config):
            BalrogSetup.runlogger.info('The function SextractorConfigs  was not found in your Balrog python config file: %s. Add this function to manually override configurations in the sextractor config file.' %BalrogSetup.pyconfig)
        else:
            config.SextractorConfigs(copy2, ExtraSexConfig)

    LogCmdlineOpts(cmdline_args, cmdline_args_copy, BalrogSetup.arglogger, '\n# Final parsed values for each command line option')
    return rules, ExtraSexConfig, cmdline_args_copy, TruthCatExtra


class Result4GSP(object):
    def __init__(self, cat):
        super(Result4GSP, self).__setattr__('galkeys', cat.galaxy.keys())
        for key in self.galkeys:
            super(Result4GSP, self).__setattr__(key, cat.galaxy[key])
    
        super(Result4GSP, self).__setattr__('compkeys', cat.component[0].keys())

        if len(cat.component)>1:
            for key in self.compkeys:
                super(Result4GSP, self).__setattr__(key, [None]*len(cat.component))
            for i in range(len(cat.component)):
                for key in self.compkeys:
                    exec 'self.%s[%i] = cat.component[%i]["%s"]' %(key,i,i,key)
        else:
            for key in self.compkeys:
                super(Result4GSP, self).__setattr__(key, cat.component[0][key])


    def __getattr__(self, name):
        raise SampledAttributeError(401, name, 'galaxies')
    
    def __getitem__(self, index):
        raise SampledIndexingError(402, 'galaxies')

    def __setitem__(self, index, value):
        raise SampledIndexingError(402, 'galaxies')

    def __setattr__(self, name, value):
        raise SampledAssignmentError(403, name, 'galaxies')


def GetSimulatedGalaxies(BalrogSetup, simgals, config, cmdline_opts_copy, TruthCatExtra):
    used = simgals.Sample(BalrogSetup)
    
    gkeys = ['minimum_fft_size','maximum_fft_size','alias_threshold','stepk_minimum_hlr','maxk_threshold','kvalue_accuracy','xvalue_accuracy','table_spacing','realspace_relerr','realspace_abserr','integration_relerr','integration_abserr']
    gsp = SimRules(BalrogSetup.ngal, gkeys, [])
    if config!=None:
        if 'GalsimParams' not in dir(config):
            BalrogSetup.runlogger.warning('The function GalsimParams  was not found in your Balrog python config file: %s. Add this function to manually override the Galsim GSParams.' %BalrogSetup.pyconfig)
        else:
            s = Result4GSP(simgals)
            config.GalsimParams(cmdline_opts_copy, gsp, s)
    grules =  [gsp.minimum_fft_size, gsp.maximum_fft_size, gsp.alias_threshold, gsp.stepk_minimum_hlr, gsp.maxk_threshold, gsp.kvalue_accuracy, gsp.xvalue_accuracy, gsp.table_spacing, gsp.realspace_relerr, gsp.realspace_abserr, gsp.integration_relerr, gsp.integration_abserr]
    gsprules = DefineRules(BalrogSetup.ngal, gkeys, grules, [], [], 0)
    used = gsprules.SimpleSample(BalrogSetup, used)

    if len(TruthCatExtra.rules)==0:
        ExtraTruthCat = None
        TruthRules = None
    else:
        ncomp = len(simgals.component)
        ExtraTruthRules = TruthCatExtra.rules
        TruthRules = DefineRules(BalrogSetup.ngal, TruthCatExtra.names, ExtraTruthRules, [], [], ncomp)

        for gal in simgals.galaxy.keys():
            TruthRules.galaxy[gal] = simgals.galaxy[gal]
        
        for i in range(ncomp):
            for comp in simgals.component[i].keys():
                TruthRules.component[i][comp] = simgals.component[i][comp]

        used = TruthRules.SimpleSample(BalrogSetup, used)

    simgals.galaxy['balrog_index'] = BalrogSetup.indexstart + np.arange(0, BalrogSetup.ngal)

    return simgals, gsprules, TruthRules, TruthCatExtra


def CompError(name, i):
    if name=='flux':
        name = 'magnitude'
    raise RulesAssignmentError(303, 'component %i of %s' %(i, name))

def GalError(name):
    raise RulesAssignmentError(303, name)



class TableColumns(object):
    def __init__(self, ngal):
        super(TableColumns, self).__setattr__('ngal', ngal)
        super(TableColumns, self).__setattr__('rules', [])
        super(TableColumns, self).__setattr__('names', [])
        super(TableColumns, self).__setattr__('fmts', [])
        super(TableColumns, self).__setattr__('units', [])
        self.InitializeSersic()

    def InitializeSersic(self, nProfiles=1):
        super(TableColumns, self).__setattr__('nProfiles', nProfiles)

    def AddColumn(self, rule=None, name=None, fmt=None, unit=None):
        rule = self._CheckRule(rule, name)
        if rule.type != 'catalog':
            if name==None:
                raise ColumnNameError(703)
        else:
            if name==None:
                name = rule.param[2]

        self.rules.append(rule)
        self.names.append(name)
        self.fmts.append(fmt)
        self.units.append(unit)


    def _CheckRule(self, rule, name):
        if type(rule).__name__!='Rule':
            if type(rule).__name__=='CompResult':
                if rule.nProfiles==1:
                    return self._CheckRule(rule[0], name)
                else:
                    raise ColumnArrayError(705, name)
            elif type(rule)==float or type(rule)==int or type(rule)==str:
                rule = Value(rule)
            else:
                try:
                    arr = np.array(rule)
                    if arr.ndim==1 and arr.size==self.ngal:
                        rule = Array(arr)
                    else:
                        raise ColumnSizeError(701, name, len(arr), self.ngal)
                except:
                    raise ColumnDefinitionError(702, name)

        return rule

    def __setattr__(self, name, value):
        raise ColumnAttributeError(706, name)


## Class used with the {sersicindex, halflightradius, magnitude, axisratio, beta} components of rules.
#  Since a simulated galaxy can have as many Sersic components as desired, the basic object of the class is an array called @p rules.
#  Index/attribute get/set methods are overwritten to put restrictions on how users can interact with the class
#  and then handle errors if they do something bad.
class CompRules(object):

    ## Initialize the rules, setting each rule to None, which will equate to using Balrog's defaults if the rule is not reassigned.
    #  @param nProfiles Integer number of Sersic profiles which make up the simulated galaxy
    #  @param name String name for the Sersic parameter, e.g. @p halflightradius
    def __init__(self, nProfiles, name):
        super(CompRules, self).__setattr__('rules', [None]*nProfiles)
        super(CompRules, self).__setattr__('name', name)
        super(CompRules, self).__setattr__('nProfiles', nProfiles)

    ## Throw an error if the user tries to define a new attribute.
    #  e.g. @p rules.beta.nonsense = 100
    #  @param name Attempted attribute name
    #  @param value Attempted attribute value
    def __setattr__(self, name, value):
        raise RulesComponentAttributeError(305)

    ## Throw an error if the user asks for an attribute which does not exist.
    #  The only attributes which exist are {rules, name, nProfiles}.
    #  However, unless they dig through the code, users will not know those three exist anyway.
    #  I would have preferred if python let me forbid access to these as well, but since
    #  they cannot be reassigned access to them does not hurt.
    #  @param name Attempted attribute name
    def __getattr__(self, name):
        raise RulesComponentAttributeError(305)
 
    ## Return the length of the @p rules array
    def __len__(self):
        return self.nProfiles

    ## Throw an error if the requested index is out of range; otherwise get element @p index of the rules array.
    #  @param index Attempted integer array position
    def __getitem__(self, index):
        if index >= self.nProfiles:
            raise RulesIndexOutOfRange(304, self.name, self.nProfiles)
        return self.rules[index]

    ## Throw an error if the requested index is out of range; otherwise check to make sure the rule given is valid before assigning the new rule or raising an exception.
    #  @param index Attempted integer array position
    #  @param value Attempted rule
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


## Class which defines @p rules.
#  The {sersicindex, halflightradius, magnitude, axisratio, beta} attributes are of type @p CompRules.
#  No set methods are allowed to be called by the user.
class SimRules(object):

    ## Initialize the simulation parameters to their balrog defaults.
    #  @param ngal Integer number of galaxies simulated
    def __init__(self, ngal, galkeys, compkeys):
        super(SimRules, self).__setattr__('ngal', ngal)
        super(SimRules, self).__setattr__('galkeys', galkeys)
        super(SimRules, self).__setattr__('compkeys', compkeys)
        for name in self._GetGalaxy():
            super(SimRules, self).__setattr__(name, None)
        self.InitializeSersic()

    ## Setup the attributes of type @p CompRules: {sersicindex, halflightradius, magnitude, axisratio, beta}.
    #  This function allocates the proper size array and will reset all rules to their defaults.
    #  @param nProfiles Integer number of Sersic profiles the simulated galaxies are composed of
    def InitializeSersic(self, nProfiles=1):
        super(SimRules, self).__setattr__('nProfiles', nProfiles)
        for c in self._GetComponent():
            super(SimRules, self).__setattr__(c, CompRules(nProfiles,c))

    def _GetGalaxy(self):
        #return ['x','y','g1','g2','magnification']
        return self.galkeys

    def _GetComponent(self):
        #return ['axisratio','beta','halflightradius','magnitude','sersicindex']
        return self.compkeys

    ## Throw an error if the user tries to index @p rules
    #  @param index Attempted integer array element position
    def __getitem__(self, index):
        raise RulesIndexingError(302)

    ## Throw an error if the user tries to index @p rules
    #  @param index Attempted integer array element position
    #  @param value Attempted assignment value
    def __setitem__(self, index, value):
        raise RulesIndexingError(302)

    ## Throw an error if the user asks for an attribute that does not exist.
    #  @param name Attempted attribute name
    def __getattr__(self, name):
        raise RulesAttributeError(301, name)

    ## Set a rule.
    #  Before new rules are assigned they are check and if found to be invalid an exception occurs.
    #  Attributes @p ngal and @p nProfiles cannot be reassigned.
    #  @param name Attempted attribute to reassign
    #  @param value Attempted assignment value
    def __setattr__(self, name, value):
        if name=='ngal':
            raise RulesHiddenError(-1, name)

        elif name=='nProfiles':
            raise RulesnProfilesError(-2, name)

        elif name=='galkeys':
            raise RulesHiddenError(-1, name)

        elif name=='compkeys':
            raise RulesHiddenError(-1, name)

        elif name in self._GetGalaxy():
            value = self._CheckRule(name, value, 'galaxy')
            super(SimRules, self).__setattr__(name, value)

        elif name in self._GetComponent():
            try:
                size = len(value)
            except:
                if self.nProfiles!=1:
                    raise RulesAssignmentNoArrayError(306)
                else:
                    size = 1
                    value = [value]
                    #print type(value[0]).__name__

            if size!=self.nProfiles:
                raise RulesAssignmentNoArrayError(306)

            for i in range(size):
                val = self._CheckRule(name, value[i], 'component', i=i)
                exec "self.%s[%i] = val" %(name, i)

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


## Class for use
#  Testing more stuff
class CompResult(object):
    def __init__(self, nProfiles, name):
        super(CompResult, self).__setattr__('name', name)
        super(CompResult, self).__setattr__('nProfiles', nProfiles)

    def __len__(self):
        return self.nProfiles

    def __getitem__(self,index):
        if index >= self.nProfiles:
            raise SampledIndexOutOfRange(404, self.name, self.nProfiles)

        if self.name=='magnitude':
            return Same( (index,'flux') )
        else:
            return Same( (index,self.name) )

    def __setitem__(self, index, value):
        raise SampledAssignmentError(403, '%s[%i]'%(self.name,index), 'sampled')

    def __setattr__(self, name, value):
        raise SampledAssignmentError(403, '%s.%s'%(self.name,name), sampled)

    def __getattr__(self, name):
        raise SampledComponentAttributeError(405)


class Results(object):
    def __init__(self, galkeys, compkeys):
        super(Results, self).__setattr__('galkeys', galkeys)
        super(Results, self).__setattr__('compkeys', compkeys)
        self.InitializeSersic()

    def InitializeSersic(self, nProfiles=1):
        super(Results, self).__setattr__('nProfiles', nProfiles)
        for c in self._GetComponent():
            super(Results, self).__setattr__(c, CompResult(nProfiles,c))

    def _GetGalaxy(self):
        #return ['x','y','g1','g2','magnification']
        return self.galkeys

    def _GetComponent(self):
        #return ['axisratio','beta','halflightradius','magnitude','sersicindex']
        return self.compkeys

    def __getattr__(self, name):
        if name not in self._GetGalaxy():
            raise SampledAttributeError(401, name, 'sampled')
        else:
            return Same(name)
    
    def __getitem__(self, index):
        raise SampledIndexingError(402, 'sampled')

    def __setitem__(self, index, value):
        raise SampledIndexingError(402, 'sampled')

    def __setattr__(self, name, value):
        raise SampledAssignmentError(403, name, 'sampled')
        
    

class DerivedArgs():
    def __init__(self,args, known):
        self.imgdir = os.path.join(args.outdir, 'balrog_image')
        self.catdir = os.path.join(args.outdir, 'balrog_cat')
        #self.logdir = os.path.join(args.outdir, 'balrog_log')
        self.sexdir = os.path.join(args.outdir, 'balrog_sexconfig')

        self.subsample = True
        if args.xmin==1 and args.ymin==1 and args.xmax==pyfits.open(args.image)[args.imageext].header['NAXIS1'] and args.ymax==pyfits.open(args.image)[args.imageext].header['NAXIS2']:
            self.subsample = False

        length = len('.sim.fits')
        self.outimageext = 0
        self.outweightext = 1
        self.imageout = DefaultName(args.image, '.fits', '.sim.fits', self.imgdir)

        self.psfout = DefaultName(args.psf, '.psf', '.psf', self.imgdir)
        self.catalogtruth = DefaultName(args.image, '.fits', '.truthcat.sim.fits', self.catdir)
        self.catalogmeasured = DefaultName(args.image, '.fits', '.measuredcat.sim.fits', self.catdir)
        if not args.noassoc:
            self.assoc_simfile = DefaultName(args.image, '.fits', '.assoc.sim.txt', self.sexdir)
            self.assoc_nosimfile = DefaultName(args.image, '.fits', '.assoc.nosim.txt', self.sexdir)

        self.psf_written = False
        self.wcshead = args.image
        ext = '.nosim.fits'
        self.nosim_imageout = '%s%s' %(self.imageout[:-length],ext)
        self.nosim_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)


        if args.image==args.weight:
            if args.nonosim:
                self.weightout = self.imageout
            else:
                self.weightout = self.nosim_imageout
        else:
            self.weightout = DefaultName(args.weight, '.fits', '.weight.fits', self.imgdir)
            self.outweightext = 0

        CreateSubDir(self.imgdir)
        CreateSubDir(self.catdir)
        CreateSubDir(self.sexdir)

        thisdir = os.path.dirname( os.path.realpath(__file__) )
        #config = os.path.join(thisdir, 'config.py')
        self.CopyFile(args.sexconfig, self.sexdir)
        self.CopyFile(args.sexparam, self.sexdir)
        self.CopyFile(args.sexnnw, self.sexdir)
        self.CopyFile(args.sexconv, self.sexdir)
        self.CopyFile(args.nosimsexparam, self.sexdir)
        if os.path.lexists(args.pyconfig):
            self.CopyFile(args.pyconfig, known.logdir)

        if args.catfitstype=='FITS_LDAC':
            self.catext = 2
        elif args.catfitstype=='FITS_1.0':
            self.catext = 1

        self.runlogger = known.logs[0]
        self.runlog = known.runlogfile
        self.sexlogger = known.logs[1]
        self.sexlog = known.sexlogfile
        self.arglogger = known.logs[2]
        self.arglog = known.arglogfile


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


def SetupLogger(known):
    
    log = logging.getLogger()
    log.setLevel(logging.NOTSET)

    runlog = logging.getLogger('run')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    known.runlogfile = os.path.join(known.logdir, 'run.log.txt')
    fh = logging.FileHandler(known.runlogfile, mode='w')
    fh.setFormatter(formatter)
    SetLevel(fh, known.logverbosity)
    runlog.addHandler(fh)
    ph = logging.StreamHandler(stream=sys.stderr)
    ph.setFormatter(formatter)
    ph = SetLevel(ph, known.stdverbosity)
    runlog.addHandler(ph)

    sexlog = logging.getLogger('sex')
    sexlog.setLevel(logging.INFO)
    known.sexlogfile = os.path.join(known.logdir, 'sex.log.txt')
    fh = logging.FileHandler(known.sexlogfile, mode='w')
    SetLevel(fh, logging.INFO)
    sexlog.addHandler(fh)

    arglog = logging.getLogger('arg')
    arglog.setLevel(logging.INFO)
    known.arglogfile = os.path.join(known.logdir, 'args.log.txt')
    fh = logging.FileHandler(known.arglogfile, mode='w')
    SetLevel(fh, logging.INFO)
    arglog.addHandler(fh)

    return [runlog, sexlog, arglog]


def ConfigureBalrog(cmdline_opts, known):
    derived_opts = DerivedArgs(cmdline_opts, known)
    BalrogSetup = BalrogConfig(cmdline_opts, derived_opts)
    return BalrogSetup


def LogDerivedOpts(cmdline_args, BalrogSetup, desc):
    ArgsDict = vars(BalrogSetup)
    VetoDict = vars(cmdline_args)
    BalrogSetup.arglogger.info(desc)
    extra_veto = ['runlogger','sexlogger','arglogger','assocnames']
    for key in ArgsDict.keys():
        if (key not in VetoDict.keys()) and (key not in extra_veto):
            BalrogSetup.arglogger.info('%s %s' %(key, ArgsDict[key]) )


def LogCmdlineOpts(cmdline_args, cmdline_args_copy, logger, desc):
    logger.info('%s' %desc)
    ArgsDict = vars(cmdline_args)
    ordered = CmdlineListOrdered()
    logger.info('# Native args')
    for key in ordered:
        logger.info('%s %s' %(key, ArgsDict[key]) )

    ArgsDict = vars(cmdline_args_copy)
    logger.info('# User-defined args')
    for key in ArgsDict.keys():
        if key not in ordered:
            logger.info('%s %s' %(key, ArgsDict[key]) ) 


def CmdlineListOrdered():
    args = ["outdir","clean",
            "pyconfig",
            "image", "imageext", "weight", "weightext", "psf", 
            "xmin", "xmax", "ymin", "ymax",
            "ngal", "seed", "gain", "zeropoint", "indexstart",
            "stdverbosity", "logverbosity", "fulltraceback", 
            "sexpath", "sexconfig", "sexparam", "sexnnw", "sexconv", "nonosim", "nosimsexparam", "noassoc", "catfitstype"]
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
    '''
    err = 0
    if not os.path.lexists(subdir):
        err = subprocess.call(['mkdir', subdir], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    if err!=0:
    '''
    if not os.path.lexists(subdir):
        try:
            os.makedirs(subdir)
        except:
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
    if not os.path.lexists(args.image):
        raise ImageInputError(101, 'image', 'image', args.image)
    if not os.path.lexists(args.weight):
        raise ImageInputError(102, 'weight', 'weight', args.weight)
    if not os.path.lexists(args.psf):
        raise ImageInputError(103, 'PSF', 'psf', args.psf)

    args.image = os.path.abspath(args.image)
    args.weight = os.path.abspath(args.weight)
    args.psf = os.path.abspath(args.psf)


def ImagesAreFITS(args, log):
    try:
        pyfits.open(args.image)
    except:
        raise FitsFileError(111, 'image', 'image', args.image) 
    try:
        pyfits.open(args.weight)
    except:
        raise FitsFileError(112, 'weight', 'weight', args.weight)
    try:
        pyfits.open(args.psf)
    except:
        raise FitsFileError(113, 'PSF', 'psf', args.psf)


def ExtsExist(args, log):
    try: 
        hdu = pyfits.open(args.image)[args.imageext]
    except:
        raise FitsExtError(121, 'image', 'imageext', args.imageext, 'image', args.image)
    try:
        hdu = pyfits.open(args.weight)[args.weightext]
    except:
        raise FitsExtError(122, 'weight', 'weightext', args.weightext, 'weight', args.weight)


def SizesOK(args, log):
    try:
        ix = pyfits.open(args.image)[args.imageext].header['NAXIS1']
    except:
        raise FitsHeaderError(131,'image','NAXIS1', 'image',args.image, 'imageext',args.imageext)
    try:
        iy = pyfits.open(args.image)[args.imageext].header['NAXIS2']
    except:
        raise FitsHeaderError(132,'image','NAXIS2', 'image',args.image, 'imageext',args.imageext)
    try:
        wx = pyfits.open(args.weight)[args.weightext].header['NAXIS1']
    except:
        raise FitsHeaderError(133,'weight','NAXIS1', 'weight',args.weight, 'weightext',args.weightext)
    try:
        wy = pyfits.open(args.weight)[args.weightext].header['NAXIS2']
    except:
        raise FitsHeaderError(134,'weight','NAXIS2', 'weight',args.weight, 'weightext',args.weightext)

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
    if args.image!=None and os.path.abspath(args.image)!=default and args.psf==None:
        raise PsfInputError(104, args.image)

    if args.image==None:
        args.image = default
        log.info('No --image explicitly given. Will use the default. Setting --image = %s' %args.image)
    if args.weight==None:
        log.info('No --weight explicitly given. Assuming it lives in same file as the image. Setting --weight = %s' %args.image)
        args.weight = args.image
    if args.weightext==None and os.path.abspath(args.weight)==args.image:
        args.weight = os.path.abspath(args.weight)
        log.info('No --weightext explicitly given. Assuming it lives in the extension after the image. Setting --weightext = %s' %(args.imageext+1))
        args.weightext = args.imageext + 1
    if args.weightext==None:
        log.info('No --weightext explicitly given. Assuming it lives in the 0th extension. Setting --weightext = 0')
        args.weightext = 0
    if args.psf==None:
        args.psf = os.path.join(indir, 'example.psf')
        log.info('No --psf explicitly given. Will use the default. Setting --psf = %s' %args.psf)
    

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
            args.gain = pyfits.open(args.image)[args.imageext].header[args.gain]
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
            args.zeropoint = pyfits.open(args.image)[args.imageext].header[args.zeropoint]
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
    args.nosimsexparam = FindSexFile(args.nosimsexparam, log, configdir, 'sex.param', 'nosimsexparam')
    args.sexnnw = FindSexFile(args.sexnnw, log, configdir, 'sex.nnw', 'sexnnw')
    args.sexconv = FindSexFile(args.sexconv, log, configdir, 'sex.conv', 'sexconv')
    args.catfitstype = 'FITS_%s' %(args.catfitstype.upper())

    try:
        sex = subprocess.check_output(['which', args.sexpath])
    except:
        raise SextractorPathError(140, args.sexpath)


def ParseDefaultArgs(args,known):
    args.pyconfig = known.pyconfig
    args.outdir = known.outdir

    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    indir = os.path.join(defdir, 'input')
    configdir = os.path.join(thisdir, 'astro_config')
    
    ParseImages(args, known.logs[0], indir)
    ParseFloatKeyword(args, known.logs[0]) 
    ParseSex(args, known.logs[0], configdir)

    return args


def DefaultArgs(parser):
    # Input and (temporary) output Images
    parser.add_argument( "-o", "--outdir", help="Directory where to put output. Output files will be named based on the input file name.", default=None)
    parser.add_argument( "-i", "--image", help="Input image to draw simulated galaxies into.", type=str, default=None)
    parser.add_argument( "-ie", "--imageext", help="FITS extension where the image lives in the input file.", type=int, default=0)
    parser.add_argument( "-w", "--weight", help="Weight map of input image.", type=str, default=None)
    parser.add_argument( "-we", "--weightext", help="FITS extension where the weight map lives in the input weight file.", type=int, default=None)
    parser.add_argument( "-p", "--psf", help="PSF of the input image.", type=str, default=None)

    # Properties you want your simulated image to have
    parser.add_argument( "-x1", "--xmin", help="Minimum column index for subsampling, >=1.", type=int, default=1)
    parser.add_argument( "-x2", "--xmax", help="Maximum column index for subsampling, <=N_pix,x.", type=int, default=-1)
    parser.add_argument( "-y1", "--ymin", help="Minimum row index for subsampling, >=1.", type=int, default=1)
    parser.add_argument( "-y2", "--ymax", help="Maximum row index for subsampling, <=N_pix,y.", type=int, default=-1)
    parser.add_argument( "-n", "--ngal", help="Number of simulated galaxies", type=int, default=50)
    parser.add_argument( "-g", "--gain", help="Gain, needed for adding noise. Can be a float or a keyword from the image header.", default='GAIN')
    parser.add_argument( "-z", "--zeropoint", help="Zeropoint used to convert simulated magnitude to flux. Sextractor runs with this zeropoint. Can be a float or a keyword from the image header.", default='SEXMGZPT')
    parser.add_argument( "-s", "--seed", help="Seed for random number generation when simulating galaxies. This does not apply to noise realizations, which are always random.", type=int, default=None)

    # Other Balrog stuff
    parser.add_argument( "-c", "--clean", help="Delete output image files", action="store_true")
    parser.add_argument( "-sv", "--stdverbosity", help="Verbosity level of stdout/stderr", type=str, default='n', choices=['n','v','vv','q'])
    parser.add_argument( "-lv", "--logverbosity", help="Verbosity level of log file", type=str, default='n', choices=['n','v','vv'])
    parser.add_argument( "-ft", "--fulltraceback", help="Full traceback is written out", action="store_true")
    parser.add_argument( "-pc", "--pyconfig", help="Balrog python config file", type=str, default=None)
    parser.add_argument( "-is", "--indexstart", help="Index for first simulated galaxy", type=int, default=0)


    # How to run sextractor
    parser.add_argument( "-sex", "--sexpath", help='Path for sextractor binary', type=str, default='sex')
    parser.add_argument( "-sc", "--sexconfig", help='Sextractor config file', type=str, default=None)
    parser.add_argument( "-sp", "--sexparam", help='Sextractor param file', type=str, default=None)
    parser.add_argument( "-sn", "--sexnnw", help='Sextractor neural network star-galaxy file', type=str, default=None)
    parser.add_argument( "-sf", "--sexconv", help='Sextractor filter convolution file', type=str, default=None)
    parser.add_argument( "-na", "--noassoc", help="Don't do association mode matching in sextractor.", action="store_true")
    parser.add_argument( "-nn", "--nonosim", help="Skip sextractor run over original image, prior to any simulation.", action="store_true")
    parser.add_argument( "-nsp", "--nosimsexparam", help="Sextractor param file for run over original image, prior to any simulation.", type=str, default=None)
    parser.add_argument( "-ct", "--catfitstype", help="Type of FITS file for sextractor to write out.", type=str, default='ldac', choices=['ldac','1.0'])


def RaiseException(log, fulltraceback=False):

    if not fulltraceback:
        exc_info = sys.exc_info()
        config_errs = []
        err_list = traceback.extract_tb(exc_info[2])
        for err in err_list: 
            file = err[0]
            if file.find('balrogexcept.py')!=-1:
                pass
            elif file.find('balrog.py')!=-1:
                pass
            elif file.find('model_class.py')!=-1:
                pass
            elif file.find('sextractor_engine.py')!=-1:
                pass
            else:
                config_errs.append(err)

        keep = traceback.format_list(config_errs)
        keep_tb = ''.join(keep)
        log.error('Run error caused Balrog to exit.\n%s' %(keep_tb), exc_info=(exc_info[0], exc_info[1], None))
        #logging.error('Run error caused Balrog to exit.\n%s' %(keep_tb), exc_info=(exc_info[0], exc_info[1], None))

    else:
        log.exception('Run error caused Balrog to exit.')
        #logging.exception('Run error caused Balrog to exit.')

    sys.exit()


def GetNativeOptions():
    parser = argparse.ArgumentParser()
    DefaultArgs(parser)
    return parser


def AddCustomOptions(parser, config, log):
    if config!=None:
        if 'CustomArgs' not in dir(config):
            log.warning('The function CustomArgs was not found in your Balrog python config file. Will continue without adding any custom command line arguments.')
        else:
            config.CustomArgs(parser) 


def NativeParse(parser, known):
    cmdline_opts = parser.parse_args()
    known.logs[2].info('# Exact command call')
    known.logs[2].info(' '.join(sys.argv))
    LogCmdlineOpts(cmdline_opts, cmdline_opts, known.logs[2], '\n# Values received for each possible command line option, filling with defaults if necessary')
    
    ParseDefaultArgs(cmdline_opts, known)
    BalrogSetup = ConfigureBalrog(cmdline_opts, known)
    return cmdline_opts, BalrogSetup


def CustomParse(cmdline_opts, BalrogSetup, config):
    galkeys = ['x', 'y', 'g1', 'g2', 'magnification']
    compkeys = ['sersicindex', 'halflightradius', 'magnitude', 'axisratio', 'beta']

    rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra = UserDefinitions(cmdline_opts, BalrogSetup, config, galkeys, compkeys)
    compkeys[2] = 'flux'

    galrules = [rules.x, rules.y, rules.g1, rules.g2, rules.magnification]
    comprules = [rules.sersicindex, rules.halflightradius, rules.magnitude, rules.axisratio, rules.beta]
    rules = DefineRules(BalrogSetup.ngal, galkeys, galrules, compkeys, comprules, rules.nProfiles )

    return rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra


def GetKnown(parser):
    known, unknown = parser.parse_known_args(sys.argv)

    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    outdir = os.path.join(defdir, 'output')
    if known.outdir==None:
        known.outdir = outdir
    #CreateDir(known.outdir)
    known.logdir = os.path.join(known.outdir, 'balrog_log')
    CreateSubDir(known.logdir)
    
    known.logs = SetupLogger(known)
    return known


def GetConfig(known):

    if known.pyconfig==None:
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        known.pyconfig = os.path.join(thisdir, 'config.py')

    if not os.path.lexists(known.pyconfig):
        #raise ConfigFileNotFound(150, known.config)
        known.logs[0].warning('Path to Balrog python config file not found: %s. All properties of the simulated galaxies will assume their defaults.' %known.pyconfig)
        config = None
    else:
        try:
            config = imp.load_source('config', known.pyconfig)
        except:
            #raise ConfigImportError(151, known.pyconfig)
            known.logs[0].warning('Python could not import your Balrog config file: %s. This means it has the python equivalent of a compiling error, apart from any possible runtime errors. First check for an error global in scope, such as an import. Continuing by assigning all properties of the simulated galaxies to their defaults.' %known.pyconfig)
            config = None

    return config

## @namespace balrog
#  Balrog stuff is in this namespace.

## The main function called to run Balrog.
#  @param parser Command line parser object made by argparse.ArgumentParser() which has had the native Balrog arguments added to it.
#  @param known  Contains a few parsed arguments, namely those needed for logging which have already been set up.
#
def RunBalrog(parser, known):

    # Find the user's config file
    config = GetConfig(known)

    # Add the user's command line options
    AddCustomOptions(parser, config, known.logs[0])

    # Parse the command line agruments and interpret the user's settings for the simulation
    cmdline_opts, BalrogSetup = NativeParse(parser, known)
    rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra = CustomParse(cmdline_opts, BalrogSetup, config)

    # Take the the user's configurations and build the simulated truth catalog out of them.
    catalog, gspcatalog, extracatalog, TruthCatExtra = GetSimulatedGalaxies(BalrogSetup, rules, config, cmdline_opts_copy, TruthCatExtra)
   
    # Get the subsampled flux and weightmap images, along with the PSF model and WCS.
    bigImage, subWeight, psfmodel, wcs = ReadImages(BalrogSetup)

    # If desired, run sextractor over the image prior to inserting any simulated galaxies.
    if not BalrogSetup.nonosim:
        NosimRunSextractor(BalrogSetup, bigImage, subWeight, extra_sex_config, catalog)

    # Insert simulated galaxies.
    bigImage = InsertSimulatedGalaxies(bigImage, catalog, psfmodel, BalrogSetup, wcs, gspcatalog)
    WriteImages(BalrogSetup, bigImage, subWeight)
    WriteCatalog(catalog, BalrogSetup, txt=None, fits=True, TruthCatExtra=TruthCatExtra, extracatalog=extracatalog)

    # Run sextractor over the simulated image.
    RunSextractor(BalrogSetup, extra_sex_config, catalog)

    # If chosen, clean up image files you don't need anymore
    if BalrogSetup.clean:
        Cleanup(BalrogSetup)

    # Log some  extra stuff Balrog used along the way
    LogDerivedOpts(cmdline_opts, BalrogSetup, '\n#Psuedo-args. Other values derived from the command line arguments.')



if __name__ == "__main__":
   
    # First get the needed info to setup up the logger, which allows everything to be logged even if things fail at very early stages.
    # Only a writable outdir is required to be able to get the output log file.
    parser = GetNativeOptions()
    known = GetKnown(parser)
    try:
        RunBalrog(parser, known)
    except:
        RaiseException(known.logs[0], fulltraceback=known.fulltraceback)

