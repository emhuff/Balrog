#!/usr/bin/env python

import shutil
import distutils.spawn

import signal
import warnings
import re
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

def WriteCatalog(sample, BalrogSetup, txt=None, fits=False, TruthCatExtra=None, extracatalog=None, setup=None):
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

        #col = pyfits.Column(name=name, array=arr,format='E', unit=unit)
        col = pyfits.Column(name=name, array=arr,format='D', unit=unit)
        columns.append(col)

    if TruthCatExtra is not None:
        for rule,name,fmt,unit in zip(TruthCatExtra.rules,TruthCatExtra.names,TruthCatExtra.fmts, TruthCatExtra.units):
            if unit is None:
                if rule.type!='catalog':
                    unit = 'unspecified'
                    BalrogSetup.runlogger.info('No unit was specificed for additional truth catalog column %s' %(name))
                else:
                    cut,hdu = FindInCat(rule, name)
                    unit = np.array(hdu.columns.units)[cut][0]
            if unit.strip()=='':
                unit = 'unspecified'

            if fmt is None:
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
            os.remove(BalrogSetup.catalogtruth)
        hdus.writeto(BalrogSetup.catalogtruth)

    if txt is not None:
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
    if BalrogSetup.nodraw and (not BalrogSetup.subsample):
        return [None]*4
    
    matchwcs = False
    image = galsim.fits.read(BalrogSetup.image, hdu=BalrogSetup.imageext)
    if image.wcs==galsim.PixelScale(1):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'fiducialwcs.fits')
        image.wcs = galsim.GSFitsWCS(file)
        matchwcs = True
        BalrogSetup.wcshead = file
        BalrogSetup.runlogger.warning('No WCS was found in the input image header. Using default pixel scale of 0.263 arcsec/pixel.')
    wcs = image.wcs

    subBounds = galsim.BoundsI(BalrogSetup.xmin,BalrogSetup.xmax,BalrogSetup.ymin,BalrogSetup.ymax)
    image = image[subBounds]
    psfmodel = galsim.des.DES_PSFEx(BalrogSetup.psf, BalrogSetup.wcshead)
    
    weight = None
    if not BalrogSetup.noweightread:
        weight = galsim.fits.read(BalrogSetup.weight, hdu=BalrogSetup.weightext)
        if matchwcs:
            weight.wcs = image.wcs
        weight = weight[subBounds]


    return image, weight, psfmodel, wcs


def WriteImages(BalrogSetup, image, weight, nosim=False, setup=None):
    weightout = BalrogSetup.weightout
    if nosim:
        imageout = BalrogSetup.nosim_imageout
    else:
        imageout = BalrogSetup.imageout


    if (BalrogSetup.nodraw) and (not BalrogSetup.subsample):
        rm_link(imageout)
        os.symlink(BalrogSetup.image, imageout)
        if (BalrogSetup.weight!=BalrogSetup.image) and (not BalrogSetup.noweightread):
            rm_link(weightout)
            os.symlink(BalrogSetup.weight, weightout)
    else:
        if BalrogSetup.weight==BalrogSetup.image:
            if not BalrogSetup.noweightread:
                galsim.fits.writeMulti(image_list=[image,weight], file_name=imageout)
            else:
                galsim.fits.writeMulti(image_list=[image], file_name=imageout)
        else:
            galsim.fits.write(image=image, file_name=imageout)
            if not BalrogSeutp.noweightread:
                if (BalrogSetup.nonosim) or (nosim):
                    if BalrogSetup.subsample:
                        galsim.fits.write(image=weight, file_name=weightout)
                    else:
                        rm_link(weightout)
                        os.symlink(BalrogSetup.weight, weightout)

    if not BalrogSetup.psf_written:
        WritePsf(BalrogSetup, BalrogSetup.psf, BalrogSetup.psfout, setup=setup)
        if BalrogSetup.detpsf!=BalrogSetup.psf:
            WritePsf(BalrogSetup, BalrogSetup.detpsf, BalrogSetup.detpsfout, setup=setup)
        BalrogSetup.psf_written = True


def WritePsf(BalrogSetup, psfin, psfout, setup=None):
    psfhdus = pyfits.open(psfin)
    psfhdus[1].header['POLZERO1'] = psfhdus[1].header['POLZERO1'] - (BalrogSetup.xmin - 1)
    psfhdus[1].header['POLZERO2'] = psfhdus[1].header['POLZERO2'] - (BalrogSetup.ymin - 1)
    if os.path.lexists(psfout):
        os.remove(psfout)
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
            if not IsNone(gspcatalog.galaxy[key]):
                if key in ['minimum_fft_size', 'maximum_fft_size', 'range_for_extrema']:
                    mod = gspcatalog.galaxy[key][i] % 1
                    if mod > 0.0:
                        d[key] = gspcatalog.galaxy[key][i]
                    else:
                        d[key] = int(gspcatalog.galaxy[key][i])
                else:
                    d[key] = gspcatalog.galaxy[key][i]

        with warnings.catch_warnings():
            try:
                warnings.simplefilter("ignore", category=galsim.GalSimDeprecationWarning)
            except:
                pass

            gsparams = galsim.GSParams(**d)

            try:
                combinedObjConv = simulatedgals.GetConvolved(psfmodel, i, wcs, gsparams, BalrogSetup)
            except:
                simulatedgals.galaxy['not_drawn'][i] = 1
                #print simulatedgals.component[0]['sersicindex'][i],simulatedgals.component[0]['halflightradius'][i],simulatedgals.component[0]['flux'][i],simulatedgals.component[0]['axisratio'][i],simulatedgals.component[0]['beta'][i], simulatedgals.galaxy['magnification'][i]; sys.stdout.flush()
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

        with warnings.catch_warnings():
            try:
                warnings.simplefilter("ignore", category=galsim.GalSimDeprecationWarning)
            except:
                pass

            try:
                smallImage = combinedObjConv.draw(scale=localscale, use_true_center=False)
            except:
                simulatedgals.galaxy['not_drawn'][i] = 1
                #print simulatedgals.component[0]['sersicindex'][i],simulatedgals.component[0]['halflightradius'][i],simulatedgals.component[0]['flux'][i],simulatedgals.component[0]['axisratio'][i],simulatedgals.component[0]['beta'][i], simulatedgals.galaxy['magnification'][i]; sys.stdout.flush()
                continue

        smallImage.setCenter(ix,iy)

        t1 = datetime.datetime.now()
        dt = t1 - t0

        if BalrogSetup.noiseseed is None:
            micro = long( (dt.days*24*60*60 + dt.seconds)*1.0e6 + dt.microseconds ) + rt + i
        else:
            micro = BalrogSetup.noiseseed + i

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


def AutoConfig(BalrogSetup, detimageout, imageout, detweightout, weightout, catalogmeasured, config_file, param_file, afile, eng, nosim):
    eng.Path(BalrogSetup.sexpath)
    eng.config['IMAGE'] = '%s[%i],%s[%s]' %(detimageout,BalrogSetup.outdetimageext, imageout,BalrogSetup.outimageext)
    eng.config['WEIGHT_IMAGE'] = '%s[%i],%s[%i]' %(detweightout,BalrogSetup.outdetweightext, weightout,BalrogSetup.outweightext)
    eng.config['CATALOG_NAME'] = catalogmeasured
    eng.config['c'] = config_file
    eng.config['PARAMETERS_NAME'] = param_file
    eng.config['STARNNW_NAME'] = BalrogSetup.sexnnw
    eng.config['FILTER_NAME'] = BalrogSetup.sexconv
    eng.config['MAG_ZEROPOINT'] = BalrogSetup.zeropoint
    eng.config['PSF_NAME'] = '%s,%s' %(BalrogSetup.detpsfout, BalrogSetup.psfout)
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
   

def RunSextractor(BalrogSetup, ExtraSexConfig, catalog, nosim=False, sim_noassoc_seg=False, setup=None):
    afile = None
    weightout = BalrogSetup.weightout
    detweightout = BalrogSetup.detweightout

    if nosim:
        catalogmeasured = BalrogSetup.nosim_catalogmeasured
        imageout = BalrogSetup.nosim_imageout
        detimageout = BalrogSetup.nosim_detimageout
        if not BalrogSetup.noassoc:
            afile = BalrogSetup.assoc_nosimfile
    elif not sim_noassoc_seg:
        catalogmeasured = BalrogSetup.catalogmeasured
        imageout = BalrogSetup.imageout
        detimageout = BalrogSetup.detimageout
        if not BalrogSetup.noassoc:
            afile = BalrogSetup.assoc_simfile
        if (BalrogSetup.image==BalrogSetup.weight) and (not BalrogSetup.noweightread):
            weightout = imageout
        if (BalrogSetup.detimage==BalrogSetup.detweight) and (not BalrogSetup.noweightread):
            detweightout = detimageout
    else:
        catalogmeasured = BalrogSetup.sim_noassoc_catalogmeasured
        imageout = BalrogSetup.imageout
        detimageout = BalrogSetup.detimageout
        if not BalrogSetup.noassoc:
            afile = BalrogSetup.assoc_simfile
        if (BalrogSetup.image==BalrogSetup.weight) and (not BalrogSetup.noweightread):
            weightout = imageout
        if (BalrogSetup.detimage==BalrogSetup.detweight) and (not BalrogSetup.noweightread):
            detweightout = detimageout

    if not BalrogSetup.noassoc:
        WriteCatalog(catalog, BalrogSetup, txt=afile, fits=False, setup=setup)
    if sim_noassoc_seg and BalrogSetup.sim_noassoc_seg_param_file:
        param_file = BalrogSetup.sim_noassoc_seg_param_file
    else:
        param_file = WriteParamFile(BalrogSetup, catalogmeasured, nosim)
    config_file = BalrogSetup.sexconfig
    #if not BalrogSetup.noassoc:
    config_file = WriteConfigFile(BalrogSetup, config_file, catalogmeasured)


    '''
    if setup.redirect is not None:
        logtosend = setup.redirect
    '''
    if setup.kind=='system':
        logtosend = BalrogSetup.sexlog
    elif setup.kind=='popen':
        logtosend = BalrogSetup.sexlogger
    sexlogsetup = setup.Copy(redirect=logtosend)

    eng = sextractor_engine.SextractorEngine(setup=sexlogsetup)

    for key in ExtraSexConfig.keys():
        eng.config[key] = ExtraSexConfig[key]

    if BalrogSetup.nonosim:
        if (BalrogSetup.nodraw) and (not BalrogSetup.subsample):
            #detweightout = BalrogSetup.detweight
            #weightout = BalrogSetup.weight
            detimageout = BalrogSetup.detimage
            imageout = BalrogSetup.image


    AutoConfig(BalrogSetup, detimageout, imageout, detweightout, weightout, catalogmeasured, config_file, param_file, afile, eng, nosim)

    if nosim:
        msg = '# Running sextractor prior to inserting simulated galaxies\n'
    else:
        msg = '\n\n# Running sextractor after inserting simulated galaxies\n'

    eng.run(msg=msg)

    if not BalrogSetup.noassoc: #and not nosim:
        CopyAssoc(BalrogSetup, catalogmeasured)


def rm_link(attr):
    if os.path.lexists(attr):
        os.remove(attr)


def NosimRunSextractor(BalrogSetup, bigImage, subweight, ExtraSexConfig, catalog, setup=None):
    if BalrogSetup.subsample:
        WriteImages(BalrogSetup, bigImage, subweight, nosim=True, setup=setup)
    else:
        rm_link(BalrogSetup.nosim_imageout)
        rm_link(BalrogSetup.psfout)
        
        os.symlink(BalrogSetup.psf, BalrogSetup.psfout)
        os.symlink(BalrogSetup.image, BalrogSetup.nosim_imageout)

        if BalrogSetup.psf!=BalrogSetup.detpsf:
            rm_link(BalrogSetup.detpsfout)
            os.symlink(BalrogSetup.detpsf, BalrogSetup.detpsfout)

        if (BalrogSetup.weight!=BalrogSetup.image) and (not BalrogSetup.noweightread):
            rm_link(BalrogSetup.weightout)
            os.symlink(BalrogSetup.weight, BalrogSetup.weightout)

        if (BalrogSetup.detweight!=BalrogSetup.detimage) and (not BalrogSetup.noweightread):
            if BalrogSetup.detweightout!=BalrogSetup.weightout:
                rm_link(BalrogSetup.detweightout)
                os.symlink(BalrogSetup.detweight, BalrogSetup.detweightout)

        if BalrogSetup.nosim_detimageout!=BalrogSetup.detimagein:
            rm_link(BalrogSetup.nosim_detimageout)
            os.symlink(BalrogSetup.detimagein, BalrogSetup.nosim_detimageout)

        BalrogSetup.psf_written = True

    RunSextractor(BalrogSetup, ExtraSexConfig, catalog, nosim=True, setup=setup)


def Cleanup(BalrogSetup, setup=None):
    files = [BalrogSetup.imageout, BalrogSetup.psfout, BalrogSetup.weightout, BalrogSetup.nosim_imageout]
    for file in files:
        if os.path.lexists(file):
            os.remove(file)


def UserDefinitions(cmdline_args, BalrogSetup, config, galkeys, compkeys):
    rules = SimRules(BalrogSetup.ngal, galkeys, compkeys)
    ExtraSexConfig = {}
    results = Results(galkeys, compkeys)
    cmdline_args_copy = copy.copy(cmdline_args)
    TruthCatExtra = TableColumns(BalrogSetup.ngal)

    if config is not None:
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
    if config is not None:
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
            if name is None:
                raise ColumnNameError(703)
        else:
            if name is None:
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
            if rule is None:
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

            if size!=self.nProfiles:
                raise RulesAssignmentNoArrayError(306)

            for i in range(size):
                val = self._CheckRule(name, value[i], 'component', i=i)
                exec "self.%s[%i] = val" %(name, i)

        else:
            raise RulesAttributeError(301,name)



    def _CheckRule(self, name, rule, kind, i=None):
        if type(rule).__name__!='Rule':

            if IsNone(rule):
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
    def __init__(self,args, known, setup=None):
        self.imgdir = os.path.join(args.outdir, 'balrog_image')
        self.catdir = os.path.join(args.outdir, 'balrog_cat')
        #self.logdir = os.path.join(args.outdir, 'balrog_log')
        self.sexdir = os.path.join(args.outdir, 'balrog_sexconfig')

        self.subsample = True
        if args.xmin==1 and args.ymin==1 and args.xmax==pyfits.open(args.image)[args.imageext].header['NAXIS1'] and args.ymax==pyfits.open(args.image)[args.imageext].header['NAXIS2']:
            self.subsample = False
        if self.subsample and args.noweightread:
            raise Exception("--noweightread and subsampling isn't possible.")

        length = len('.sim.fits')
        #dlength = len('.sim.det.fits')
        dlength = len('.fits')

        self.outimageext = 0
        self.outdetimageext = 0
        self.outweightext = 1
        self.outdetweightext = 1

        self.imageout = DefaultName(args.image, '.fits', '.sim.fits', self.imgdir)
        #self.detimageout = DefaultName(args.image, '.fits', '.sim.det.fits', self.imgdir)

        ext = '.nosim.fits'
        dext = '.nosim.det.fits'
        self.nosim_imageout = '%s%s' %(self.imageout[:-length],ext)
        #self.nosim_detimageout = '%s%s' %(self.detimageout[:-dlength],dext)
        #self.nosim_detimageout = DefaultName(args.detimage, '.fits', '.nosim.det.fits', self.imgdir)

        if args.detimage==args.image:
            self.detimagein = args.image
            self.detimageout = self.imageout
            #self.nosim_detimageout = self.detimagein
            self.nosim_detimageout = self.nosim_imageout
        else:
            self.detimagein = args.detimage
            self.detimageout = args.detimage
            self.nosim_detimageout = DefaultName(args.detimage, '.fits', '.nosim.det.fits', self.imgdir)
            if self.subsample:
                raise Exception('subsampling with different measurement and deteciton images is not supported currently, though may be some day.')

        if args.weight!=args.detweight:
            if self.subsample:
                raise Exception('subsampling with different measurement and deteciton images is not supported currently, though may be some day.')


        self.psfout = DefaultName(args.psf, '.psf', '.psf', self.imgdir)
        if args.psf!=args.detpsf:
            self.detpsfout = DefaultName(args.psf, '.psf', '.det.psf', self.imgdir)
        else:
            args.detpsfout = args.psf


        self.catalogtruth = DefaultName(args.image, '.fits', '.truthcat.sim.fits', self.catdir)
        self.catalogmeasured = DefaultName(args.image, '.fits', '.measuredcat.sim.fits', self.catdir)
        if not args.noassoc:
            self.assoc_simfile = DefaultName(args.image, '.fits', '.assoc.sim.txt', self.sexdir)
            self.assoc_nosimfile = DefaultName(args.image, '.fits', '.assoc.nosim.txt', self.sexdir)            
        self.psf_written = False
        self.wcshead = args.image



        self.nosim_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)
        if args.sim_noassoc_seg is not None:
            ext = '.sim_noassoc.fits'
            self.sim_noassoc_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)

        if args.image==args.weight:
            if args.nonosim:
                self.weightout = self.imageout
            else:
                self.weightout = self.nosim_imageout
        else:
            self.weightout = DefaultName(args.weight, '.fits', '.weight.fits', self.imgdir)
            self.outweightext = 0

        if args.detimage==args.detweight:
            if args.nonosim:
                self.detweightout = self.detimageout
            else:
                self.detweightout = self.nosim_detimageout
        else:
            #self.detweightout = DefaultName(args.detweight, '.fits', '.weight.det.fits', self.imgdir)
            self.detweightout = args.detweight
            self.outdetweightext = 0


        if args.noweightread:
            self.weightout = args.weight
            self.outweightext = args.weightext
            self.detweightout = args.detweight
            self.outdetweightext = args.detweightext

        CreateSubDir(self.imgdir)
        CreateSubDir(self.catdir)
        CreateSubDir(self.sexdir)

        thisdir = os.path.dirname( os.path.realpath(__file__) )
        #config = os.path.join(thisdir, 'config.py')
        self.CopyFile(args.sexconfig, self.sexdir, setup=setup)
        self.CopyFile(args.sexparam, self.sexdir, setup=setup)
        self.CopyFile(args.sexnnw, self.sexdir, setup=setup)
        self.CopyFile(args.sexconv, self.sexdir, setup=setup)
        self.CopyFile(args.nosimsexparam, self.sexdir, setup=setup)
        if os.path.lexists(args.pyconfig):
            self.CopyFile(args.pyconfig, known.logdir, setup=setup)

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


    def CopyFile(self, file, dir, setup=None):
        basename = os.path.basename(file)
        outname = os.path.join(dir,basename)
        if os.path.exists(outname):
            os.remove(outname)
        shutil.copy(file, outname)

        
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

    #r = np.random.random()
    r = time.time()
    i = known.indexstart
    runlog = logging.getLogger('run_%f_%i'%(r,i))

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

    sexlog = logging.getLogger('sex_%f'%(r))
    sexlog.setLevel(logging.INFO)
    known.sexlogfile = os.path.join(known.logdir, 'sex.log.txt')
    fh = logging.FileHandler(known.sexlogfile, mode='w')
    SetLevel(fh, logging.INFO)
    sexlog.addHandler(fh)

    '''
    if os.path.exists(known.sexlogfile):
        os.remove(known.sexlogfile)
    '''

    arglog = logging.getLogger('arg_%f'%(r))
    arglog.setLevel(logging.INFO)
    known.arglogfile = os.path.join(known.logdir, 'args.log.txt')
    fh = logging.FileHandler(known.arglogfile, mode='w')
    SetLevel(fh, logging.INFO)
    arglog.addHandler(fh)

    return [runlog, sexlog, arglog]


def ConfigureBalrog(cmdline_opts, known, setup=None):
    derived_opts = DerivedArgs(cmdline_opts, known, setup=setup)
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
    args = ["outdir",
            "image", "imageext", "detimage", "detimageext",
            "weight", "weightext", "detweight", "detweightext",
            "psf","detpsf", 
            "xmin", "xmax", "ymin", "ymax","ngal", "gain", "zeropoint", "seed", "noiseseed",
            "imageonly","nodraw","clean","stdverbosity", "logverbosity", "fulltraceback", "pyconfig","indexstart",
            "sexpath", "sexconfig", "sexparam", "sexnnw", "sexconv", "noassoc", "nonosim", "nosimsexparam",  "catfitstype",
            "sim_noassoc_seg", "sim_noassoc_seg_param_file",
            "systemcmd", "retrycmd", "useshell"]

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
    if not os.path.lexists(subdir):
        try:
            os.makedirs(subdir)
        except:
            raise SubdirWriteError(202, subdir)


def get_check_array_ext(args):
    attrs = [args.image, args.weight, args.psf, args.detimage, args.detweight, args.detpsf]
    ss = ['image', 'weight', 'psf', 'detimage', 'detweight','detpsf']
    eattrs = [args.imageext, args.weightext, 0, args.detimageext, args.detweightext, 0]
    ess = ['imageext', 'weightext', 'psfext', 'detimageext', 'detweightext', 'detpsfext']
    codes = np.array( [1,2,3, 4,5,6] )
    return attrs, ss, eattrs, ess, codes

def images_existence_check(attr, code, s):
    if not os.path.lexists(attr):
        raise ImageInputError(code, s, s, attr)

def images_are_fits(attr, code, s):
    try:
        pyfits.open(attr)
    except:
        raise FitsFileError(code, s, s, attr) 

def exts_exist(attr, sattr, ext, sext, code):
    try: 
        hdu = pyfits.open(attr)[ext]
    except:
        raise FitsExtError(code, sattr, sext, ext, sattr, attr)


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
    if args.image is not None and os.path.abspath(args.image)!=default and args.psf is None:
        raise PsfInputError(104, args.image)

    if args.image is None:
        args.image = default
        log.info('No --image explicitly given. Will use the default. Setting --image = %s' %args.image)
    if args.detimage is None:
        args.detimage = args.image
        log.info('No --detimage explicitly given. Will use the image itself. Setting --detimage = %s' %args.detimage)
    if args.weight is None:
        args.weight = args.image
        log.info('No --weight explicitly given. Assuming it lives in same file as the image. Setting --weight = %s' %args.weight)
    if args.detweight is None:
        if os.path.abspath(args.image)!=os.path.abspath(args.weight):
            args.detweight = args.weight
        else:
            args.detweight = args.detimage
        log.info('No --detweight explicitly given. Assuming it lives in same file as the detimage. Setting --detweight = %s' %args.detimage)

    if args.weightext is None and os.path.abspath(args.weight)==os.path.abspath(args.image):
        args.weight = os.path.abspath(args.weight)
        args.weightext = args.imageext + 1
        log.info('No --weightext explicitly given. Assuming it lives in the extension after the image. Setting --weightext = %s' %(args.weightext))
    if args.weightext is None:
        args.weightext = 0
        log.info('No --weightext explicitly given. Assuming it lives in the 0th extension. Setting --weightext = 0')
    if args.detweightext is None and os.path.abspath(args.detweight)==os.path.abspath(args.detimage):
        args.detweight = os.path.abspath(args.detweight)
        args.detweightext = args.detimageext + 1
        log.info('No --detweightext explicitly given. Assuming it lives in the extension after the detimage. Setting --detweightext = %s' %(args.detweightext))
    if args.detweightext is None:
        args.detweightext = 0
        log.info('No --detweightext explicitly given. Assuming it lives in the 0th extension. Setting --detweightext = 0')

    if args.psf is None:
        args.psf = os.path.join(indir, 'example.psf')
        log.info('No --psf explicitly given. Will use the default. Setting --psf = %s' %args.psf)
    if args.detpsf is None:
        args.detpsf = args.psf
        log.info('No --detpsf explicitly given. Assuming it is the same file as the image psf. Setting --detpsf = %s' %args.detpsf)
    

def ParseImages(args, log, indir):
    FindImages(args, log, indir)

    args.image = os.path.abspath(args.image)
    args.weight = os.path.abspath(args.weight)
    args.psf = os.path.abspath(args.psf)
    args.detimage = os.path.abspath(args.detimage)
    args.detweight = os.path.abspath(args.detweight)
    args.detpsf = os.path.abspath(args.detpsf)

    attrs, ss, eattrs, ess, codes = get_check_array_ext(args)
    for attr,s, eattr,es, code in zip(attrs,ss, eattrs,ess, codes):
        images_existence_check(attr, code, s)
        images_are_fits(attr, code+10, s)
        #if s!='psf':
        exts_exist(attr, s, eattr, es, code+20)

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
    if arg is None:
        arg = default
        log.info('No --%s explicitly given. Will use the default. Setting --%s = %s' %(label,label,arg))
    if not os.path.lexists(arg):
        log.warning('--%s %s does not exist. Assuming a default file. Setting --%s = %s' %(label, arg, label, default))
        arg = default
    return arg


def SysInfoPrint(setup, msg, level='warning', exc_info=None, skip=False):

    if (not skip) and (setup.redirect is not None):
        if setup.kind=='system':
            time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
            with open(setup.redirect, 'a') as f:
                f.write('%s: %s\n'%(time, msg))
        elif setup.kind=='popen':
            if level=='warning':
                setup.redirect.warning(msg)
            elif level=='info':
                setup.redirect.info(msg)
            elif level=='error':
                #setup.redirect.error(msg, exc_info=exc_info)
                setup.redirect.error(msg, exc_info=exc_info)


def SystemCall(oscmd, setup=None, delfiles=[], keeps=[], timeout=None):
    if setup is None:
        raise Exception('system call exception')


    if setup.kind=='popen':
        kwargs = {}
        if setup.useshell:
            kwargs = {'shell':True}
        if setup.redirect is not None:
            kwargs['stdout'] = subprocess.PIPE
            kwargs['stderr'] = subprocess.PIPE

    if (setup.kind=='system') or (setup.useshell):
        syscmd = subprocess.list2cmdline(oscmd)
        if (setup.redirect is not None) and (setup.kind=='system'):
            syscmd = '%s >> %s 2>&1'%(syscmd, setup.redirect)
    else:
        syscmd = oscmd


    done = False
    attempt = 0
    while not done:
        
        '''
        if attempt==1:
            SysInfoPrint(setup, 'First retry: Doing system command:\n%s\n'%(syscmd), level='info')
        '''

        if (attempt%1000==0):
            skip = False
        else:
            skip = True
        SysInfoPrint(setup, 'Attempt %i: Doing system command:\n%s\n'%(attempt, syscmd), level='info', skip=skip)

        if setup.kind=='system': 
            rcode = os.system(syscmd)

        elif setup.kind=='popen':
            p = subprocess.Popen(syscmd, preexec_fn=os.setsid, **kwargs)

            if timeout is not None:
                dt = 0
                t1 = time.time()
                while True:
                    status = p.poll()
                    if status is not None:
                        SysInfoPrint(setup, 'Timed command exited after %.2f seconds'%(dt), level='info', skip=skip)
                        break
                    if dt > timeout:
                        os.killpg(os.getpgid(p.pid), signal.SIGTERM)
                        SysInfoPrint(setup, 'Command timed out after %.2f seconds'%(dt), level='info', skip=skip)
                        break
                    dt = time.time() - t1

            stdout, stderr = p.communicate()
            rcode = p.returncode
            SysInfoPrint(setup, 'stdout:\n%s\n'%(stdout), level='info', skip=skip)
            SysInfoPrint(setup, 'stderr:\n%s\n'%(stderr), level='info', skip=skip)
        SysInfoPrint(setup, 'returncode: %s'%(str(rcode)), level='info', skip=skip)

        if rcode==0:
            done = True
        else:
            SysInfoPrint(setup, 'System call failed', skip=skip)
            if not setup.retry:
                done = True
            else:
                SysInfoPrint(setup, 'Retrying the command', skip=skip)

            for file in delfiles:
                dir = os.path.dirname(file)
                base = os.path.basename(file)
                files = os.listdir(dir)
                for f in files:
                    fullpath = os.path.join(dir,f)
                    if fullpath in keeps:
                        continue
                    if f.startswith(base):
                        os.remove(fullpath)

        if done and skip:
            SysInfoPrint(setup, 'Finished on attempt %i\n'%(attempt), level='info')
            if setup.kind=='popen':
                SysInfoPrint(setup, 'stdout:\n%s\n'%(stdout), level='info')
                SysInfoPrint(setup, 'stderr:\n%s\n'%(stderr), level='info')
                SysInfoPrint(setup, 'returncode: %s'%(str(rcode)), level='info')
                
        attempt += 1

    return rcode


def ParseSex(args, log, configdir, setup=None):
    args.sexconfig = FindSexFile(args.sexconfig, log, configdir, 'sex.config', 'sexconfig')
    args.sexparam = FindSexFile(args.sexparam, log, configdir, 'bulge.param', 'sexparam')
    args.nosimsexparam = FindSexFile(args.nosimsexparam, log, configdir, 'sex.param', 'nosimsexparam')
    args.sexnnw = FindSexFile(args.sexnnw, log, configdir, 'sex.nnw', 'sexnnw')
    args.sexconv = FindSexFile(args.sexconv, log, configdir, 'sex.conv', 'sexconv')
    args.catfitstype = 'FITS_%s' %(args.catfitstype.upper())

    which = distutils.spawn.find_executable(args.sexpath)
    if which is None:
        raise SextractorPathError(140, args.sexpath)


def ParseDefaultArgs(args,known, setup=None):
    args.pyconfig = known.pyconfig
    args.outdir = known.outdir

    thisdir = os.path.dirname( os.path.realpath(__file__) )
    defdir = os.path.join(thisdir, 'default_example')
    indir = os.path.join(defdir, 'input')
    configdir = os.path.join(thisdir, 'astro_config')
    
    ParseImages(args, known.logs[0], indir)
    ParseFloatKeyword(args, known.logs[0]) 
    ParseSex(args, known.logs[0], configdir, setup=setup)

    args.syslog = setup

    return args


def DefaultArgs(parser):
    # Input and (temporary) output Images
    parser.add_argument( "-o", "--outdir", help="Directory where to put output. Output files will be named based on the input file name.", default=None)

    parser.add_argument( "-i", "--image", help="Input image to draw simulated galaxies into.", type=str, default=None)
    parser.add_argument( "-ie", "--imageext", help="FITS extension where the image lives in the input file.", type=int, default=0)
    parser.add_argument( "-di", "--detimage", help="Input detection image ", type=str, default=None)
    parser.add_argument( "-die", "--detimageext", help="FITS extension where the detection image lives in the input file.", type=int, default=0)

    parser.add_argument( "-w", "--weight", help="Weight map of input image.", type=str, default=None)
    parser.add_argument( "-we", "--weightext", help="FITS extension where the weight map lives in the input weight file.", type=int, default=None)
    parser.add_argument( "-dw", "--detweight", help="Weight map of input detection image.", type=str, default=None)
    parser.add_argument( "-dwe", "--detweightext", help="FITS extension where the detection weight map lives in the input weight file.", type=int, default=None)
    parser.add_argument( "-nwr", "--noweightread", help="Skip any weight map reading/writing, because you don't actually need to. But as a consequence you get less nicely organized image output.", action="store_true")

    parser.add_argument( "-p", "--psf", help="PSF of the input image.", type=str, default=None)
    parser.add_argument( "-dp", "--detpsf", help="PSF of the input detimage.", type=str, default=None)

    # Properties you want your simulated image to have
    parser.add_argument( "-x1", "--xmin", help="Minimum column index for subsampling, >=1.", type=int, default=1)
    parser.add_argument( "-x2", "--xmax", help="Maximum column index for subsampling, <=N_pix,x.", type=int, default=-1)
    parser.add_argument( "-y1", "--ymin", help="Minimum row index for subsampling, >=1.", type=int, default=1)
    parser.add_argument( "-y2", "--ymax", help="Maximum row index for subsampling, <=N_pix,y.", type=int, default=-1)
    parser.add_argument( "-n", "--ngal", help="Number of simulated galaxies", type=int, default=50)
    parser.add_argument( "-g", "--gain", help="Gain, needed for adding noise. Can be a float or a keyword from the image header.", default='GAIN')
    parser.add_argument( "-z", "--zeropoint", help="Zeropoint used to convert simulated magnitude to flux. Sextractor runs with this zeropoint. Can be a float or a keyword from the image header.", default='SEXMGZPT')
    parser.add_argument( "-s", "--seed", help="Seed for random number generation when simulating galaxies. This does not apply to noise realizations.", type=int, default=None)
    parser.add_argument( "-nse", "--noiseseed", help="Seed for random number generation noise", type=int, default=None)

    # Other Balrog stuff
    parser.add_argument( "-io", "--imageonly", help="Only write the image, don't run sextractor", action="store_true")
    parser.add_argument( "-nd", "--nodraw", help="Don't draw generated galaxies into image", action="store_true")
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
    parser.add_argument( "-nas", "--sim_noassoc_seg", help="If set, run non assoctation mode sextractor over simulated image, as well as assoc mode,\
        and save resulting seg map to outdir/sim_noassoc_seg", type=str, default=None)
    parser.add_argument( "--sim_noassoc_seg_param_file", help="param file for nas run (since only need it for seg map, this can be minimal)", default=None)
    parser.add_argument( "-nn", "--nonosim", help="Skip sextractor run over original image, prior to any simulation.", action="store_true")
    parser.add_argument( "-nsp", "--nosimsexparam", help="Sextractor param file for run over original image, prior to any simulation.", type=str, default=None)
    parser.add_argument( "-ct", "--catfitstype", help="Type of FITS file for sextractor to write out.", type=str, default='ldac', choices=['ldac','1.0'])


    # System call type stuff
    parser.add_argument( "-syscmd", "--systemcmd", help="How to do system calls", type=str, default='popen', choices=['popen','system'])
    parser.add_argument( "-rt", "--retrycmd", help="Retry system commands is they fail", action='store_true')
    parser.add_argument( "-ush", "--useshell", help="Do Popen with shell=True", action='store_true')


def RaiseException(log, fulltraceback=False, sendto=None, message=None):
    exc_info = sys.exc_info()
    err_list = traceback.extract_tb(exc_info[2])

    if not fulltraceback:
        config_errs = []
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
        if keep_tb!='':
            keep_tb = '\n%s'%(keep_tb)

        if message is None:
            msg = 'Run error caused Balrog to exit.%s'%(keep_tb)
        else:
            msg = '%s%s'%(message,keep_tb)

        exc_info = (exc_info[0], exc_info[1], None)
        if sendto is not None:
            m = '%s' %(msg)
            if sendto.kind=='system':
                s = logging.Formatter().formatException(exc_info)
                m = '%s\n%s'%(m, s)
            SysInfoPrint(sendto, m, level='error', exc_info=exc_info)
        if log is not None:
            log.error(msg, exc_info=exc_info)

    else:
        if message is None:
            msg = 'Run error caused Balrog to exit.'
        else:
            msg = message

        if sendto is not None:
            m = '%s'%(msg)
            if sendto.kind=='system':
                s = logging.Formatter().formatException(exc_info)
                m = '%s\n%s'%(m, s)
            SysInfoPrint(sendto, m, level='error', exc_info=exc_info)
        if log is not None:
            log.exception(msg)


def GetNativeOptions():
    parser = argparse.ArgumentParser()
    DefaultArgs(parser)
    return parser


def AddCustomOptions(parser, config, log):
    if config is not None:
        if 'CustomArgs' not in dir(config):
            log.warning('The function CustomArgs was not found in your Balrog python config file. Will continue without adding any custom command line arguments.')
        else:
            config.CustomArgs(parser) 


def NativeParse(parser, known, setup=None):
    cmdline_opts = parser.parse_args()
    known.logs[2].info('# Exact command call')
    known.logs[2].info(' '.join(sys.argv))
    LogCmdlineOpts(cmdline_opts, cmdline_opts, known.logs[2], '\n# Values received for each possible command line option, filling with defaults if necessary')
    
    ParseDefaultArgs(cmdline_opts, known, setup=setup)
    BalrogSetup = ConfigureBalrog(cmdline_opts, known, setup=setup)
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
    if known.outdir is None:
        known.outdir = outdir
    known.logdir = os.path.join(known.outdir, 'balrog_log')
    CreateSubDir(known.logdir)
    
    known.logs = SetupLogger(known)
    return known


def GetConfig(known):

    if known.pyconfig is None:
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


def ImgOpenIfNeeded(setup, BalrogSetup, extra_sex_config, catalog, gspcatalog, TruthCatExtra, extracatalog):
    SysInfoPrint(setup, 'Reading images', level='info')
    # Get the subsampled flux and weightmap images, along with the PSF model and WCS.
    bigImage, subWeight, psfmodel, wcs = ReadImages(BalrogSetup)


    SysInfoPrint(setup, 'Nosim sex', level='info')
    # Get the subsampled flux and weightmap images, along with the PSF model and WCS.
    if not BalrogSetup.imageonly:
        # If desired, run sextractor over the image prior to inserting any simulated galaxies.
        if not BalrogSetup.nonosim:
            NosimRunSextractor(BalrogSetup, bigImage, subWeight, extra_sex_config, catalog, setup=setup)

    # Insert simulated galaxies.
    if not BalrogSetup.nodraw:
        SysInfoPrint(setup, 'Insert objects', level='info')
        bigImage = InsertSimulatedGalaxies(bigImage, catalog, psfmodel, BalrogSetup, wcs, gspcatalog)
        SysInfoPrint(setup, 'Write out', level='info')
        WriteImages(BalrogSetup, bigImage, subWeight, setup=setup)
        SysInfoPrint(setup, 'Write catalog', level='info')
        WriteCatalog(catalog, BalrogSetup, txt=None, fits=True, TruthCatExtra=TruthCatExtra, extracatalog=extracatalog, setup=setup)
    else:
        SysInfoPrint(setup, 'Write out', level='info')
        WriteImages(BalrogSetup, bigImage, subWeight, setup=setup)

def RunBalrog(parser, known, setup):
   
    SysInfoPrint(setup, 'Starting RunBalrog', level='info')
    # Find the user's config file
    config = GetConfig(known)

    SysInfoPrint(setup, 'Reading user options', level='info')
    # Add the user's command line options
    AddCustomOptions(parser, config, known.logs[0])

    SysInfoPrint(setup, 'Parsing options', level='info')
    # Parse the command line agruments and interpret the user's settings for the simulation
    cmdline_opts, BalrogSetup = NativeParse(parser, known, setup=setup)
    rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra = CustomParse(cmdline_opts, BalrogSetup, config)

    SysInfoPrint(setup, 'Getting simulated parameters', level='info')
    # Take the the user's configurations and build the simulated truth catalog out of them.
    catalog, gspcatalog, extracatalog, TruthCatExtra = GetSimulatedGalaxies(BalrogSetup, rules, config, cmdline_opts_copy, TruthCatExtra)
  
    # Do it like this so the image doesn't need to be open when sim sextractor is running
    ImgOpenIfNeeded(setup, BalrogSetup, extra_sex_config, catalog, gspcatalog, TruthCatExtra, extracatalog)

    SysInfoPrint(setup, 'Sim sex', level='info')
    if not BalrogSetup.imageonly:

        # Run sextractor over the simulated image.
        RunSextractor(BalrogSetup, extra_sex_config, catalog, setup=setup)

        # Run sextractor over simulated image in noassoc mode, only saving seg mask
        if BalrogSetup.sim_noassoc_seg:
            BalrogSetup.noassoc = True
            extra_sex_config['CHECKIMAGE_TYPE']='SEGMENTATION'
            extra_sex_config['CHECKIMAGE_NAME']=BalrogSetup.sim_noassoc_seg
            RunSextractor(BalrogSetup, extra_sex_config, catalog, sim_noassoc_seg=True, setup=setup)

    SysInfoPrint(setup, 'Cleanup', level='info')
    # If chosen, clean up image files you don't need anymore
    if BalrogSetup.clean:
        Cleanup(BalrogSetup, setup=setup)

    SysInfoPrint(setup, 'Writing some logs', level='info')
    # Log some  extra stuff Balrog used along the way
    LogDerivedOpts(cmdline_opts, BalrogSetup, '\n#Psuedo-args. Other values derived from the command line arguments.')



class SystemCallSetup(object):

    def Copy(self, redirect=None):
        if redirect is None:
            raise Exception('Must give redirect')
        return SystemCallSetup(retry=self.retry, redirect=redirect, kind=self.kind, useshell=self.useshell)

    def __init__(self, retry=False, redirect=None, kind='popen', useshell=False):
        self.retry = retry
        self.redirect = redirect
        self.kind = kind
        self.useshell = useshell



def BalrogFunction(args=None, syslog=None):

    if args is not None:
        argv = [os.path.realpath(__file__)]
        for arg in args:
            argv.append(arg)
        sys.argv = argv

    # First get the needed info to setup up the logger, which allows everything to be logged even if things fail at very early stages.
    # Only a writable outdir is required to be able to get the output log file.
    
    parser = GetNativeOptions()
    known = GetKnown(parser)

    SystemSetup = SystemCallSetup(retry=known.retrycmd, redirect=syslog, kind=known.systemcmd, useshell=known.useshell)
    try:
        RunBalrog(parser, known, SystemSetup)
        ret = 0
    except Exception as e:
        RaiseException(known.logs[0], fulltraceback=known.fulltraceback, sendto=SystemSetup)
        ret = 1

    return ret


if __name__ == "__main__":
   
    ret = BalrogFunction()
    sys.exit(ret)
