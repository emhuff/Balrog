from balrog import *
import os
import argparse
import sys
import shutil
import astropy.table
import astropy.io.fits as pyfits
import multi_epoch_tools.coords as coords
import multi_epoch_tools.file_tools as file_tools
import multi_epoch_tools.des_meds as des_meds
#import pylab
#import psutil
#import functools
import copy
import subprocess
from mpi4py import MPI

MAGZP_REF=30.
BALROG_DESDATA=os.environ['BALROG_DESDATA']
thisdir=os.path.dirname(os.path.realpath(__file__))
RunConfig={'swarp_config_default':os.path.join(thisdir,"astro_config/sva1/default.swarp"),
'swarp':"swarp"}

def Mkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError:
        if not os.path.isdir(dirname):
            raise

def combined_seg(seg_orig_file, seg_new_file):
    """We normally run SExtractor in association mode on image simulated with balrog
    The downside to this is that the segmentation map produced does not inculde any of
    the objects in the original image. So combine the original seg map with the new one,
    by setting all pixels in the new seg map to -1, if there is an object detected there
    in the original seg map. If a pixel in the new seg map had a detection in the old seg map
    do not set to -1, but create a new blended_with_orig flag

    returns:
    seg - new combiend seg map
    blended_with_orig - list of simulated object ids which are blended with old objects
    """
    #Read original seg map
    seg=galsim.fits.read(seg_orig_file,hdu=1)
    #Read new seg map
    #seg_new=galsim.fits.read(BalrogSetup.extra_sex_config['CHECKIMAGE_NAME'])
    seg_new=galsim.fits.read(seg_new_file)
    assert seg.array.shape==seg_new.array.shape
    new_obj_ids=np.unique(seg_new.array)
    #First find new objects which are blended with old ones
    blended_with_orig_pixels=(seg_new.array>0)*(seg.array>0)
    blended_with_orig=np.unique(seg_new.array[blended_with_orig_pixels]) #list of new ids which are blended with original objects

    #Start with original seg map, setting all detections to -1
    seg.array[seg.array>0]=-1
    #Then simply overwrite pixels where there are new detections
    seg.array[seg_new.array>0]=seg_new.array[seg_new.array>0]
    return seg, blended_with_orig

class SEBalrogSetup(dict):
    def __init__(self,CoaddSetup,srclist_entry,funpack=True):
        #Use same logging files as CoaddSetup        
        self.runlogger = CoaddSetup.runlogger
        self.sexlogger = CoaddSetup.sexlogger
        self.arglogger = CoaddSetup.arglogger
        self.gspcatalog = CoaddSetup.gspcatalog
        self.ngal = CoaddSetup.ngal

        #Now overwrite the relevant ones
        img_path=srclist_entry['red_image']
        self['sky_path']=srclist_entry['red_bkg']
        if funpack:
            img_path=img_path.strip(".fz")
            if not os.path.isfile(img_path):
                file_tools.funpack(srclist_entry['red_image'])
                file_tools.funpack(srclist_entry['red_bkg'])
            self['sky_path']=self['sky_path'].strip(".fz")            
            img_suffix=".fits"
            self['skyext']=0
            self['imageext']=0
            self['weightext']=2
            self['compression']=None
        else:
            img_suffix='.fits.fz'
            self['imageext']=1
            self['weightext']=3
            self['compression']='rice'
            self['skyext']=1

        imgdir_orig=os.path.dirname(img_path)
        img_basename=os.path.basename(img_path)
        imgdir=imgdir_orig.replace(os.environ['DESDATA'],BALROG_DESDATA)
        try:
            os.makedirs(imgdir)
        except OSError:
            if not os.path.isdir(imgdir):
                raise
        sexconfigdir=os.path.join(imgdir,'_sexconfig')
        try:
            os.mkdir(sexconfigdir)
        except OSError:
            if not os.path.isdir(sexconfigdir):
                raise    
        self['detimageext']=self['imageext']
        self['imgdir']=imgdir
        self['zeropoint']=srclist_entry['magzp']
        self['weight']=img_path
        self['image']=img_path
        self['weightout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.nosim.fits'))
        self['imageout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.sim.fits'))
        self['detimageout']=self['imageout']
        self['psfout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.psf'))
        self['detpsf']= img_path.replace(img_suffix,'_psfcat.psf')
        self['catalogmeasured']=os.path.join(imgdir,img_basename.replace(img_suffix,'.measuredcat.sim.fits'))
        self['detimagein']=img_path
        self['catalogtruth']=os.path.join(imgdir,img_basename.replace(img_suffix,'.truthcat.sim.fits'))
        self['assoc_simfile']=os.path.join(imgdir,img_basename.replace(img_suffix,'.assoc.sim.txt'))
        self['detweight']=img_path
        self['assoc_nosimfile']=os.path.join(imgdir,img_basename.replace(img_suffix,'.assoc.nosim.txt'))
        self['psf']= img_path.replace(img_suffix,'_psfcat.psf')
        if not CoaddSetup.no_noise:
            self['gain'] = srclist_entry['gaina']
        else:
            self['gain'] = None
        self['sexdir']= sexconfigdir
        self['scampheader']=srclist_entry['astro_refine']
        self['detimage']= self['image']
        self['detpsfout']= self['psfout']
        self['nosim_detimageout']= self['image']
        self['nosim_catalogmeasured']=os.path.join(imgdir,img_basename.replace(img_suffix,'.measuredcat.nosim.fits'))
        self['nosim_imageout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.nosim.fits'))
        self['xmin']=1
        self['ymin']=1
        self['xmax']=pyfits.open(self['image'])[self['imageext']].header['NAXIS1']
        self['ymax']=pyfits.open(self['image'])[self['imageext']].header['NAXIS2']
        for key,val in self.iteritems():
            setattr(self,key,val)
        #if self['catfitstype']=='FITS_LDAC':
        #    self['catext'] = 2
        #elif self['catfitstype']=='FITS_1.0':
        #    self['catext'] = 1

        #Set keys as attributes:
        for key,val in self.iteritems():
            setattr(self,key,val)

    def ReadImages(self):
        _,self.input_hdulist,_=galsim.fits.readFile(self.image)
        image = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.imageext,compression=self.compression)
        weight = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.weightext,compression=self.compression) 
        badpix = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.imageext+1,compression=self.compression)
        #If reading from a compressed image, the badpix mask has, for some reason 32768 added to it (or whatever 'BZERO' is in 
        #the header). Also galsim doesn't seem to know its an int...
        if self.compression is not None:
            bzero=pyfits.getheader(self.image,self.imageext+1)['BZERO']
            badpix = galsim.ImageI(badpix.array.astype(int) - bzero)
        sky = galsim.fits.read(file_name=self.sky_path, hdu=self.skyext, compression=self.compression)
        if self.scampheader is not None:
            h=pyfits.Header.fromtextfile(self.scampheader)
            wcs=galsim.wcs.readFromFitsHeader(h)[0]
        else:
            wcs=image.wcs
        image.wcs = wcs
        weight.wcs = image.wcs

        subBounds = galsim.BoundsI(self.xmin,self.xmax,self.ymin,self.ymax)
        image = image[subBounds]
        weight = weight[subBounds]
        badpix = badpix[subBounds]
        psfmodel = galsim.des.DES_PSFEx(self.psf, wcs=wcs)
        return image, weight, badpix, sky, psfmodel, wcs

    def WriteImages(self, image, nosim=False):
        if nosim:
            imageout = self.nosim_imageout
        else:
            imageout = self.imageout

        #Replace image hdu data
        self.input_hdulist[self.imageext].data=np.zeros_like(image.array)
        #Write to imageout
        galsim.fits.writeFile(imageout,self.input_hdulist)

        if not self.psf_written:
            WritePsf(self, self.psf, self.psfout)
            if self.detpsf!=self.psf:
                WritePsf(self, self.detpsf, self.detpsfout)
            self.psf_written = True

def GetNativeOptions():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    TileArgs(parser)
    DefaultArgs(parser)
    return parser

def GetKnown(parser):
    known, unknown = parser.parse_known_args(sys.argv[1:])
    return known


def NativeParse(parser, known):
    cmdline_opts = parser.parse_args()
    known.logs[2].info('# Exact command call')
    known.logs[2].info(' '.join(sys.argv))
    LogCmdlineOpts(cmdline_opts, cmdline_opts, known.logs[2], '\n# Values received for each possible command line option, filling with defaults if necessary')
    #SetupCoadd(cmdline_opts,known.logs[0], known)
    MiscSetup(cmdline_opts, known.logs[0], known)
    #BalrogSetup = ConfigureBalrog(cmdline_opts, known)
    return cmdline_opts

def MiscSetup(args, log, known):
    args.pyconfig = known.pyconfig
    args.outdir = known.outdir
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    configdir = os.path.join(thisdir, 'astro_config')
    #ParseFloatKeyword(args, known.logs[0]) 
    ParseSex(args, known.logs[0], configdir)

class CoaddSetup(dict):
    def __init__(self, cmdline_opts, known):
        opt_dict=vars(cmdline_opts)
        known_dict=vars(known)
        for key, val in known_dict.iteritems():
            self[key]=val
        for key, val in opt_dict.iteritems():
            self[key]=val
        for key,val in self.iteritems():
            setattr(self, key, val)
        self.init_logs()

    def init_logs(self):
        self.runlogger = self.logs[0]
        self.runlog = self.runlogfile
        self.sexlogger = self.logs[1]
        self.sexlog = self.sexlogfile
        self.arglogger = self.logs[2]
        self.arglog = self.arglogfile
            
    def init_coadd(self, coadd):

        #Setup coadd image 
        try:
            self.image = coadd['image_url_funpacked']
        except KeyError:
            coadd['image_url_funpacked']=(coadd['image_url']).rsplit('.',1)[0]
            self.image = coadd['image_url_funpacked']
        #self.detimage = self.image
        self.imageext = 0
        #self.detimageext = 0
        self.weight = self.image
        self.weightext = 1
        #self.detweight = self.image
        #self.detweightext = 1
        self.psf = coadd['psf_path']
        self.detpsf = self.psf

        #Also set path to segmentation map
        imgdir = os.path.dirname(coadd['image_url'])
        imgbase = os.path.basename(coadd['image_url'])
        self.seg = os.path.join(os.path.dirname(imgdir),'QA','segmap',imgbase.replace('.fits.fz','_seg.fits.fz'))

        """
        attrs, ss, eattrs, ess, codes = get_check_array_ext(self)
        for attr,s, eattr,es, code in zip(attrs,ss, eattrs,ess, codes):
            print attr,s, eattr,es, code
            images_existence_check(attr, code, s)
            images_are_fits(attr, code+10, s)
            #if s!='psf':
            exts_exist(attr, s, eattr, es, code+20)
        """
        #Get zeropoint and gain from header
        self.zeropoint=pyfits.open(self.image)[self.imageext].header['SEXMGZPT']
        self.gain=pyfits.open(self.image)[self.imageext].header['GAIN']

        SizesOK(self, self.logs[0])
        self.DerivedArgs()

    def DerivedArgs(self):
        self.imgdir = os.path.join(self.outdir, 'balrog_image')
        self.catdir = os.path.join(self.outdir, 'balrog_cat')
        #self.logdir = os.path.join(self.outdir, 'balrog_log')
        self.sexdir = os.path.join(self.outdir, 'balrog_sexconfig')

        self.subsample = True
        if self.xmin==1 and self.ymin==1 and self.xmax==pyfits.open(self.image)[self.imageext].header['NAXIS1'] and self.ymax==pyfits.open(self.image)[self.imageext].header['NAXIS2']:
            self.subsample = False

        length = len('.sim.fits')
        #dlength = len('.sim.det.fits')
        dlength = len('.fits')

        self.outimageext = 0
        #self.outdetimageext = 0
        self.outweightext = 1
        #self.outdetweightext = 1

        self.imageout = DefaultName(self.image, '.fits', '.sim.fits', self.imgdir)
        self.segout=self.imageout.replace('.fits','_seg.fits')
        #self.detimageout = DefaultName(self.image, '.fits', '.sim.det.fits', self.imgdir)
        self.psfout = DefaultName(self.psf, '.psf', '.psf', self.imgdir)

        self.catalogtruth = DefaultName(self.image, '.fits', '.truthcat.sim.fits', self.catdir)
        self.catalogmeasured = DefaultName(self.image, '.fits', '.measuredcat.sim.fits', self.catdir)
        if not self.noassoc:
            self.assoc_simfile = DefaultName(self.image, '.fits', '.assoc.sim.txt', self.sexdir)
            self.assoc_nosimfile = DefaultName(self.image, '.fits', '.assoc.nosim.txt', self.sexdir)            
        self.psf_written = False
        self.wcshead = self.image

        ext = '.nosim.fits'
        self.nosim_imageout = '%s%s' %(self.imageout[:-length],ext)
        #self.nosim_detimageout = '%s%s' %(self.detimageout[:-dlength],dext)
        #self.nosim_detimageout = DefaultName(self.detimage, '.fits', '.nosim.det.fits', self.imgdir)


        self.nosim_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)
        if self.sim_noassoc_seg is not None:
            ext = '.sim_noassoc.fits'
            self.sim_noassoc_catalogmeasured = '%s%s' %(self.catalogmeasured[:-length],ext)

        if self.image==self.weight:
            if self.nonosim:
                self.weightout = self.imageout
            else:
                self.weightout = self.nosim_imageout
        else:
            self.weightout = DefaultName(self.weight, '.fits', '.weight.fits', self.imgdir)
            self.outweightext = 0

        CreateSubDir(self.imgdir)
        CreateSubDir(self.catdir)
        CreateSubDir(self.sexdir)

        thisdir = os.path.dirname( os.path.realpath(__file__) )
        #config = os.path.join(thisdir, 'config.py')
        shutil.copy(self.sexconfig, self.sexdir)
        shutil.copy(self.sexparam, self.sexdir)
        shutil.copy(self.sexnnw, self.sexdir)
        shutil.copy(self.sexconv, self.sexdir)
        shutil.copy(self.nosimsexparam, self.sexdir)
        if os.path.lexists(self.pyconfig):
            shutil.copy(self.pyconfig, self.logdir)

        if self.catfitstype=='FITS_LDAC':
            self.catext = 2
        elif self.catfitstype=='FITS_1.0':
            self.catext = 1


    def init_det(self, detimage, detweight, detimageext=None, detweightext=None):
        self.detimage=detimage
        self.detweight=detweight
        if detimageext is None:
            self.detimageext = 0
        if detweightext is None:
            if self.detweight==self.detimage:
                self.detweightext = 1
            else:
                self.detweightext = 0
        self.outdetimageext = self.detimageext
        self.outdetweightext = self.detweightext

        if self.detimage==self.image:
            self.detimagein = self.image
            self.detimageout = self.imageout
            self.nosim_detimageout = self.detimagein
        else:
            self.detimagein = self.detimage
            self.detimageout = self.detimage
            self.nosim_detimageout = DefaultName(self.detimage, '.fits', '.nosim.det.fits', self.imgdir)

        if self.psf!=self.detpsf:
            self.detpsfout = DefaultName(self.psf, '.psf', '.det.psf', self.imgdir)
        else:
            self.detpsfout = self.psf

        if self.detimage==self.detweight:
            if self.nonosim:
                self.detweightout = self.detimageout
            else:
                self.detweightout = self.nosim_detimageout
        else:
            self.detweightout = self.detweight # DefaultName(self.detweight, '.fits', '.weight.det.fits', self.imgdir)
            self.outdetweightext = 0

        attrs, ss, eattrs, ess, codes = get_check_array_ext(self)
        for attr,s, eattr,es, code in zip(attrs,ss, eattrs,ess, codes):
            #print attr,s, eattr,es, code
            images_existence_check(attr, code, s)
            images_are_fits(attr, code+10, s)
            #if s!='psf':
            exts_exist(attr, s, eattr, es, code+20)

def CustomParse(cmdline_opts, BalrogSetup, config, bands):
    galkeys = ['x', 'y', 'g1', 'g2', 'magnification']
    compkeys = ['sersicindex', 'halflightradius', 'axisratio', 'beta', 'flux']
    for band in bands:
    	compkeys.append('magnitude_%s'%band)

    rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra = UserDefinitions(cmdline_opts, BalrogSetup, config, galkeys, compkeys)
    rules.flux=0
    #compkeys[4] = 'flux'

    galrules = [rules.x, rules.y, rules.g1, rules.g2, rules.magnification]
    comprules = [rules.sersicindex, rules.halflightradius, rules.axisratio, rules.beta, rules.flux]
    for band in bands:
        comprules.append(getattr(rules, 'magnitude_%s'%band))
    rules = DefineRules(BalrogSetup.ngal, galkeys, galrules, compkeys, comprules, rules.nProfiles )

    return rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra

def mag2flux(mag, zp):
    return np.power(10.0, (zp - mag) / 2.5)

def SwarpConfig(imgs, RunConfig, BalrogSetup, iext=0, wext=1):
    config = {'RESAMPLE': 'N',
              'COMBINE': 'Y',
              'COMBINE_TYPE': 'CHI-MEAN',
              'SUBTRACT_BACK': 'N',
              'DELETE_TMPFILES': 'Y',
              'WEIGHT_TYPE': 'MAP_WEIGHT',
              'PIXELSCALE_TYPE': 'MANUAL',
              'PIXEL_SCALE': str(0.270),
              'CENTER_TYPE': 'MANUAL',
              'HEADER_ONLY': 'N',
              'WRITE_XML': 'N'}

    header = pyfits.open(imgs[0])[iext].header
    xsize = header['NAXIS1']
    ysize = header['NAXIS2']
    config['IMAGE_SIZE'] = '%i,%i' %(xsize,ysize)
    xc = header['CRVAL1']
    yc = header['CRVAL2']
    config['CENTER'] = '%f,%f' %(xc,yc)

    ims = []
    ws = []
    for i in range(len(imgs)):
        ims.append( '%s[%i]' %(imgs[i],iext) )
        ws.append( '%s[%i]' %(imgs[i],wext) )
    ims = ','.join(ims)
    ws = ','.join(ws)

    dir = os.path.join(BalrogSetup['outdir'], 'det')
    Mkdir(dir)
    imbase=os.path.basename(imgs[0])
    imout = os.path.join(dir, imbase.replace(BalrogSetup.detbands[0]+'.sim','det.sim'))
    wout = imout.replace('.fits', '_weight.fits')
    config['IMAGEOUT_NAME'] = imout
    config['WEIGHTOUT_NAME'] = wout
    
    call = [RunConfig['swarp'], ims, '-c', RunConfig['swarp_config_default'], '-WEIGHT_IMAGE', ws]
    for key in config:
        call.append( '-%s'%(key) )
        call.append( config[key] )
    #call = ' '.join(call)
    return call, imout, wout

def get_stamp(image,x_pos,y_pos,stamp_size,offset=(0,0),scale=1.,int_stamp=False):
    "Get postage stamp to run im3shape on"
    "Expects galsim image, SExtractor X_IMAGE,Y_IMAGE and stamp size as input"
    "Returns stamp, x0 and y0 - coordinates of object in stamp (starting from 0)"
    "and galsim bounds object"
    half=stamp_size/2
    x_center, y_center = int(x_pos),int(y_pos)
    #x_center,y_center = int(x_pos+0.5+offset[0]),int(y_pos+0.5+offset[1])
    x_min,x_max,y_min,y_max = x_center-half,x_center+half,y_center-half,y_center+half
    out_of_bounds = (x_pos<image.bounds.xmin or x_pos>image.bounds.xmax or y_pos<image.bounds.ymin or y_pos>image.bounds.ymax)
    if out_of_bounds: raise OutOfBoundsError(0,0,x_pos,y_pos)
    subBounds = galsim.BoundsI(x_min,x_max-1,y_min,y_max-1)
    if not int_stamp:
        stamp = galsim.ImageD(stamp_size,stamp_size)
    else:
        stamp = galsim.ImageI(stamp_size,stamp_size)
    stamp.setCenter(x_center,y_center)
    bounds=stamp.bounds & image.bounds
    stamp[bounds]=image[bounds] * scale
    x0,y0 = x_pos - x_min, y_pos - y_min    
    return stamp, x0, y0, subBounds

def coadd_extract(coaddsetup, catalog, sim_img, sim_weight, sim_wcs_list):

    coaddsetup.extra_sex_config['CHECKIMAGE_TYPE']='SEGMENTATION'
    coaddsetup.extra_sex_config['CHECKIMAGE_NAME']=coaddsetup.segout
    RunSextractor(coaddsetup, coaddsetup.extra_sex_config, catalog)

    seg, blended_with_orig=combined_seg(coaddsetup.seg, coaddsetup.extra_sex_config['CHECKIMAGE_NAME'])

    sex_data=astropy.table.Table(pyfits.getdata(coaddsetup.catalogmeasured, coaddsetup.catext))
    #replace vector_assoc with balrog_index
    balrog_index_det=sex_data['VECTOR_ASSOC'][:,coaddsetup.assocnames.index('balrog_index')].astype(int)
    col=astropy.table.Column(balrog_index_det, name='balrog_index')
    sex_data.remove_column('VECTOR_ASSOC')
    sex_data.add_column(col)
    coaddsetup.runlogger.info("Ran SExtractor")

    #Make new column specifying whether object is blended with original object (using seg map)
    blended_with_orig=list(blended_with_orig)
    new_sex_ids=np.array(sex_data['NUMBER'])
    col_data=np.zeros_like(new_sex_ids)
    for i,new_id in enumerate(new_sex_ids):
        if new_id in blended_with_orig:
            col_data[i]=1        
    blended_with_orig_col=astropy.table.Column(col_data, name='blended_with_orig')
    sex_data.add_column(blended_with_orig_col)
    #Overwrite sextractor catalog with this aumented data
    sex_data.write(coaddsetup.catalogmeasured, overwrite=True)

    Cleanup(coaddsetup)

    Mkdir(coaddsetup.sim_outdir)
    for cat in [coaddsetup.catalogtruth,coaddsetup.catalogmeasured]:
        shutil.move(cat, coaddsetup.sim_outdir)
    coaddsetup.catalogtruth=os.path.join(coaddsetup.sim_outdir,os.path.basename(coaddsetup.catalogtruth))
    coaddsetup.catalogmeasured=os.path.join(coaddsetup.sim_outdir,os.path.basename(coaddsetup.catalogmeasured))
    """
    try:
        shutil.copytree(coaddsetup.outdir,coaddsetup.sim_outdir)
    except OSError:
        shutil.rmtree(coaddsetup.sim_outdir)
        shutil.copytree(coaddsetup.outdir,coaddsetup.sim_outdir)
    """
    coaddsetup.runlogger.info("Moved output files to sim ouput dir")

    #Now make meds inputs
    #Set boxsize to be same for all objects for now
    try: 
        box_size=known.box_size
    except AttributeError:
        box_size=[48]*known.ngal
    catalog.galaxy['box_size']=box_size

    #If making meds on the fly, need to get image, weight and badpix cutouts
    #Also need to give each meds entry the correct 'number' - this the SExtractor object number,
    #and should correspond to the value of the detection pixels in the seg map
    balrog_index_list=list(balrog_index_det)
    sex_numbers=np.array(sex_data['NUMBER'].astype(int))

    meds_inputs=[]
    for i in range(coaddsetup.ngal):
        meds_inputs.append(des_meds.MedsInput(id=coaddsetup.indexstart+i))

    if not coaddsetup.no_meds:
        for i,(x,y) in enumerate(zip(catalog.galaxy['x'],catalog.galaxy['y'])):
            try:
                sex_cat_ind=balrog_index_list.index(catalog.galaxy['balrog_index'][i])
                number=sex_numbers[sex_cat_ind]
            except ValueError:
                number=0
            meds_inputs[i]['number']=number
            if catalog.galaxy['not_drawn'][i]:
                continue
            try:
                image_stamp,x0,y0,bounds=get_stamp(sim_img,x,y,box_size[i])
            except OutOfBoundsError:
                continue
            weight_stamp,_,_,_=get_stamp(sim_weight,x,y,box_size[i])
            seg_stamp,_,_,_=get_stamp(seg,x,y,box_size[i])
            badpix_stamp=galsim.ImageI(np.zeros_like(image_stamp.array))
            meds_inputs[i].add(image_stamp,weight_stamp,badpix_stamp,seg_stamp,sim_wcs_list[i],y-1,x-1,bounds.ymin-1,bounds.xmin-1,y0,x0,0)
    return sex_data, meds_inputs, coaddsetup.catalogtruth, coaddsetup.catalogmeasured

def SERun(catalog, se_info, BalrogSetup):

    x,y=coords.get_pcoords(catalog.galaxy['ra'],catalog.galaxy['dec'],BalrogSetup.wcs)
    catalog.galaxy['x'],catalog.galaxy['y']=x,y
    # Insert simulated galaxies.
    catalog.galaxy['not_drawn']=np.zeros_like(catalog.galaxy['not_drawn'])
    #print catalog.component[0]['flux']
    catalog.component[0]['flux'] = mag2flux(catalog.component[0]['magnitude_%s'%se_info['band']], BalrogSetup.zeropoint)
    #print catalog.component[0]['flux']
    bigImage, wcs_list = InsertSimulatedGalaxies(BalrogSetup.bigImage.copy(), 
                                                 catalog, BalrogSetup.psfmodel, 
                                                 BalrogSetup, BalrogSetup.wcs, BalrogSetup.gspcatalog, 
                                                 return_wcs_list=True, out_of_bounds_threshold=48)
    #If making meds on the fly, need to get image, weight and badpix cutouts
    #First calculate scaling factor to set to magzp_ref
    scale=10**(0.4*(MAGZP_REF-se_info['magzp']))
    se_objects=[]
    for i_obj,(x,y) in enumerate(zip(catalog.galaxy['x'],catalog.galaxy['y'])):
        if catalog.galaxy['not_drawn'][i_obj]:
            continue
        try:
            image_stamp,x0,y0,bounds=get_stamp(bigImage,x,y,catalog.galaxy['box_size'][i_obj])
        except OutOfBoundsError:
            BalrogSetup.runlogger.info("object %d OutOfBoundsError, x=%f,y=%f"%(i_obj,x,y))
            continue
        weight_stamp,_,_,_=get_stamp(BalrogSetup.weight,x,y,catalog.galaxy['box_size'][i_obj])
        if known.no_noise:
            w=1000**2/(np.sum(image_stamp.array**2))
            weight_stamp=galsim.ImageD(catalog.galaxy['box_size'][i_obj], catalog.galaxy['box_size'][i_obj], init_value=np.float(w))
        badpix_stamp,_,_,_=get_stamp(BalrogSetup.badpix,x,y,catalog.galaxy['box_size'][i_obj], int_stamp=True)
        sky_stamp,_,_,_=get_stamp(BalrogSetup.sky,x,y,catalog.galaxy['box_size'][i_obj])
        #sky subtract, rescale image stamp and set badpix to zero
        if not known.no_noise:
            image_stamp-=sky_stamp
        image_stamp*=scale
        weight_stamp/=(scale*scale)
        #print 'i_obj',i_obj
        image_stamp_array=image_stamp.array
        #print image_stamp_array
        badpix_array=badpix_stamp.array
        #print badpix_array
        image_stamp_array[badpix_array!=0]=0.
        #print 'i_obj',i_obj
        #print image_stamp_array
        image_stamp=galsim.ImageD(image_stamp_array)
        seg_stamp=None
        se_object=des_meds.SEObject(image_stamp,weight_stamp,badpix_stamp,seg_stamp,wcs_list[i_obj],y-1,x-1,
                               bounds.ymin-1,bounds.xmin-1,y0,x0,se_info['file_id'])
        se_objects.append((i_obj, se_object))

    return se_objects

def MultiBandCoaddRun(coaddsetup, sim_inds=[0], slave=False, comm=None, master=0, shears=[[0.,0.]]):
    #First read in image etc. - if doing more than 1 sim, don't want to do this more than once
    ImageOrig=[]
    for i_band,band in enumerate(coaddsetup.bands):
        coaddsetup.init_coadd(coaddsetup.coadds[i_band])
        bigImage, subWeight, psfmodel, wcs = ReadImages(coaddsetup)
        #Copy big image since this will get overwritten by InsertSimulatedGalaxies
        ImageOrig.append((bigImage.copy(), subWeight, psfmodel, wcs))

    catalogs_allsims, sex_data_allbands_allsims, meds_inputs_allbands_allsims = [],[],[]
    for i_sim in sim_inds:
        coaddsetup.runlogger.info('sim %d',i_sim)
        coaddsetup.indexstart = i_sim * coaddsetup.ngal
        first_shear=True
        catalogs=[]
        for i_shear,shear in enumerate(shears):
            print 'shear',shear
            print 'catalogs',catalogs
            det_imgouts=[]
            meds_inputs_allbands={}
            sex_data_allbands={}
            coadd_images=[]
            coadd_weights=[]
            coadd_wcs_lists=[]

            #For each band, get catalog of simulated galaxies (only the flux column changes currently), inserting
            #simulated galaxies, and saving images.
            #band_setups=[]
            first_band=True
            coaddsetup.runlogger.info('Getting simulated galaxy catalogs, and making simulated coadd images')
            for i_band,band in enumerate(coaddsetup.bands):
                print 'band',band
                coaddsetup.runlogger.info('band: %s'%band)
                coaddsetup.init_coadd(known.coadds[i_band])
                if first_shear:
                    if first_band:
                        coaddsetup.rules, coaddsetup.extra_sex_config, coaddsetup.cmdline_opts_copy, coaddsetup.TruthCatExtra = CustomParse(coaddsetup, coaddsetup, coaddsetup.config, coaddsetup.bands)
                        coaddsetup.magcol='magnitude_%s'%band
                        catalog, gspcatalog, extracatalog, TruthCatExtra = GetSimulatedGalaxies(coaddsetup, coaddsetup.rules, 
                                                                                        coaddsetup.config, 
                                                                                        coaddsetup.cmdline_opts_copy, 
                                                                                        coaddsetup.TruthCatExtra)
                        coaddsetup.gspcatalog=gspcatalog
                        first_band
                    else:
                        catalog.component[0]['flux']=mag2flux(catalog.component[0]['magnitude_%s'%band],coaddsetup.zeropoint)

                else:
                    catalog=catalogs[i_band]
                
                catalog.galaxy['g1'],catalog.galaxy['g2']=shear[0]*np.ones(coaddsetup.ngal),shear[1]*np.ones(coaddsetup.ngal)

                #Insert galaxies
                print 'inserting galaxies'
                bigImage, subWeight, psfmodel, wcs = ImageOrig[i_band]
                if coaddsetup.no_noise:
                    bigImage *= 0 
                coadd_x,coadd_y=catalog.galaxy['x'],catalog.galaxy['y']
                catalog.galaxy['ra'],catalog.galaxy['dec']=coords.get_wcoords(coadd_x,coadd_y,wcs)
                bigImage,wcs_list_coadd = InsertSimulatedGalaxies(bigImage, catalog, psfmodel, 
                                                              coaddsetup, wcs, 
                                                              gspcatalog, return_wcs_list=True)
                #Save coadd image and truth catalog
                print 'writing image'
                coadd_images.append(bigImage)
                coadd_weights.append(subWeight)
                coadd_wcs_lists.append(wcs_list_coadd)
                WriteImages(coaddsetup, bigImage, subWeight)
                if band in coaddsetup.detbands:
                    det_imgouts.append(coaddsetup.imageout)
                WriteCatalog(catalog, coaddsetup, txt=None, fits=True, TruthCatExtra=coaddsetup.TruthCatExtra, extracatalog=extracatalog)      
                catalogs.append(copy.deepcopy(catalog))
                #band_setups.append(coaddsetup)

            first_shear=False

            #Now run swarp if more than one det band, else just set det image
            print det_imgouts
            coaddsetup.runlogger.info('det images: %s'%str(det_imgouts))
            if len(coaddsetup.detbands)>1:
                swarp_call, detimg, detweight = SwarpConfig(det_imgouts, RunConfig, coaddsetup)
                coaddsetup.runlogger.info('swarp call:')
                coaddsetup.runlogger.info(str(swarp_call))
                subprocess.call(swarp_call)
            else:
                detimg,detweight = det_imgouts[0],det_imgouts[0]

            coaddsetup.runlogger.info('Now SExtracting and initialising MEDS inputs')
            truthcat_files,measuredcat_files=[],[]
            for i_band,band in enumerate(coaddsetup.bands):
                print 'band',band
                print 'running coadd extract'
                catalog, sim_img, sim_weight, sim_wcs_list = catalogs[i_band], coadd_images[i_band], coadd_weights[i_band], coadd_wcs_lists[i_band]    
                coaddsetup.init_coadd(known.coadds[i_band])
                coaddsetup.init_det(detimg, detweight)
                coaddsetup.sim_outdir = os.path.join(coaddsetup.parent_output,coaddsetup.tilename+".%d.%d"%(i_sim,i_shear))
                sex_data, meds_inputs, truthcat_file, measuredcat_file = coadd_extract(coaddsetup, catalog, sim_img, sim_weight, sim_wcs_list)
                truthcat_files.append(truthcat_file)
                measuredcat_files.append(measuredcat_file)
                meds_inputs_allbands[band] = meds_inputs
                sex_data_allbands[band] = sex_data
                #coadd_not_drawn=np.copy(catalog.galaxy['not_drawn'])
                #coaddsetup.runlogger.info("%d objects not drawn"%(coadd_not_drawn).sum())
                #Cleanup(coaddsetup)

            #Also clean up det image and weight (if multiband)
            if len(coaddsetup.detbands) > 1:
                for file in [detimg,detweight]:
                    os.remove(file)

            #print meds_inputs_allbands

            if slave:
                #info=(i_sim,i_shear,,)
                info={'i_sim':i_sim,
                      'i_shear':i_shear,
                      'truthcat_files':truthcat_files,
                      'measuredcat_files':measuredcat_files}
                comm.send((copy.deepcopy(catalogs), sex_data_allbands, meds_inputs_allbands, info), dest=master, tag=i_sim)
            else:
                 catalogs_allsims.append(catalogs)
                 sex_data_allbands_allsims.append(sex_data_allbands)
                 meds_inputs_allbands_allsims.append(meds_inputs_allbands)

    if slave and i_sim==sim_inds[-1]:
        #print 'rank %d returning gspcatalog'%(rank)
        return gspcatalog
    else:
        return catalogs_allsims, sex_data_allbands_allsims, meds_inputs_allbands_allsims, gspcatalog

def Run(parser, known, comm=None):

    print 'rank %d starting Run()'%(comm.Get_rank())

    #First step is to read in the coadd and generate the truth catalog based on this:
    # Find the user's config file
    config = GetConfig(known)

    # Add the user's command line options
    AddCustomOptions(parser, config, known.logs[0])

    #Setup coadd object to get coadd path etc.
    if comm is not None:
        rank = comm.Get_rank()
        if rank==0:            
            known.coadds=[file_tools.setup_tile(known.tilename, band=band, sync=False, funpack=True) for band in known.bands]
        else:
            known.coadds=[file_tools.setup_tile(known.tilename, band=band, sync=False, funpack=False) for band in known.bands]
            print 'rank %d got coadd info'%(comm.Get_rank())
        comm.Barrier()
    else:
        known.coadds=[file_tools.setup_tile(known.tilename, band=band, sync=False, funpack=True) for band in known.bands]

    #Check shear_perturb settings
    known.shears=[[0.,0.]]
    if known.shear_perturb is not None:
        known.shears.append(known.shear_perturb)
    known.n_shears=len(known.shears)

    # Parse the command line arguments and interpret the user's settings for the simulation
    cmdline_opts = NativeParse(parser, known)
    coaddsetup_base = CoaddSetup(cmdline_opts, known)
    coaddsetup_base.config = config
    coaddsetup_base.runlogger=coaddsetup_base.logs[0]

    #Get number of processes. If no mpi, set rank to 0 and everything else *should* work
    master_slave=False
    master=False
    if comm is not None:
        rank = comm.Get_rank()
        size = comm.Get_size()
        if size>1:
            master_slave=True
            if rank==0:
                master=True
                print 'Doing master-worker, with %d workers'%(size-1)
    else:
        rank, size = 0,1    

    # Do coadd part, first known.coadd_nproc processes do this
    sim_inds=np.arange(known.n_sims)
    coadd_outputs=[]
    idle=False
    if master_slave:
        catalogs_allsims, sex_data_allbands_allsims, meds_inputs_allbands_allsims=[],[],[]
        if (known.coadd_nproc is None) or (known.coadd_nproc>size-1):
            coadd_size = size-1
        else:
            coadd_size = known.coadd_nproc
        gspcatalog=None
        if not master:
            catalogs_allsims,sim_nums,shear_nums=None,None,None
            #gspcatalog=None
            local_sim_inds=np.array_split(sim_inds,coadd_size)[rank-1]
            print 'rank,local_sim_inds',rank, local_sim_inds
            if len(local_sim_inds)>0:
                coaddsetup_base.runlogger.info('Doing simulations: '+str(local_sim_inds))
                gspcatalog=MultiBandCoaddRun(coaddsetup_base, sim_inds=local_sim_inds, slave=True, comm=comm, 
                                             shears=known.shears)
            else:
                idle=True
                coaddsetup_base.runlogger.info("Fewer sims than processes, I'm idle...")
        else:
            """
            while len(catalogs_allsims)<known.n_sims:
                status=MPI.Status()
                rec_data=comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                #print rec_data
                tag = status.Get_tag()
                coaddsetup_base.runlogger.info('reveived data from sim %d'%tag)
                catalogs_allsims.append(rec_data[0])
                sex_data_allbands_allsims.append(rec_data[1])
                meds_inputs_allbands_allsims.append(rec_data[2])
            """
            rec_datas={}
            sims_collected=0
            while sims_collected<known.n_sims*known.n_shears:
                status=MPI.Status()
                rec_data=comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                info=rec_data[-1]
                rec_i_sim,rec_i_shear=info['i_sim'],info['i_shear']
                #i_sim_rec = status.Get_tag()
                coaddsetup_base.runlogger.info('reveived data from sim %d, shear %d'%(rec_i_sim,rec_i_shear))
                rec_datas[rec_i_sim,rec_i_shear]=rec_data
                sims_collected+=1
                #catalogs_allsims.append(rec_data[0])
                #sex_data_allbands_allsims.append(rec_data[1])
                #meds_inputs_allbands_allsims.append(rec_data[2])

            coaddsetup_base.runlogger.info('received data from all sims')
            sim_nums,shear_nums=[],[]
            truthcats,measuredcats=[],[]
            for i_sim in range(known.n_sims):
                for i_shear in range(known.n_shears):
                    catalogs_allsims.append(rec_datas[i_sim,i_shear][0])
                    sex_data_allbands_allsims.append(rec_datas[i_sim,i_shear][1])
                    meds_inputs_allbands_allsims.append(rec_datas[i_sim,i_shear][2])
                    sim_nums.append(i_sim)
                    shear_nums.append(i_shear)
                    truthcats.append(rec_datas[i_sim,i_shear][3]['truthcat_files'])
                    measuredcats.append(rec_datas[i_sim,i_shear][3]['measuredcat_files'])

            #print 'catalogs_allsims, sex_data_allbands_allsims, meds_inputs_allbands_allsims'    
            #print catalogs_allsims, sex_data_allbands_allsims, meds_inputs_allbands_allsims

        print 'checkpoint 1',comm.Get_rank()
        catalogs_allsims=comm.bcast(catalogs_allsims, root=0)

        gspcatalog=comm.bcast(gspcatalog, root=1)
        if idle:
            print rank, 'received gspcatalog'
            coaddsetup_base.gspcatalog=gspcatalog

        comm.Barrier()
    else:
        catalogs_allsims, sex_data_allbands_allsims, meds_inputs_allbands_allsims,gspcatalog=MultiBandCoaddRun(coaddsetup_base, sim_inds=sim_inds, slave=False, comm=None)
    
    #Now do single-epoch stuff
    #Loop through bands, and then through single-epoch images
    coaddsetup_base.runlogger.info('Now for the single-epoch images...')
    for i_band,band in enumerate(known.bands):
        print rank,'Now doing SE, band:',band
        coaddsetup_base.runlogger.info('band: '+band)
        coadd = known.coadds[i_band]
        coadd.srclist.sort(key=lambda x: x['red_cat'])
        se_image_inds=np.arange(len(coadd.srclist))
        file_ids=1+np.arange(len(coadd.srclist))

        if (master_slave):
            if not master:
                local_se_image_inds=np.array_split(se_image_inds,size-1)[rank-1]
                coaddsetup_base.runlogger.info('doing se images: '+str(local_se_image_inds))
                local_srclist=[coadd.srclist[i] for i in local_se_image_inds]
                local_file_ids=file_ids[local_se_image_inds]
                se_objects=[]
                for i_se, se_info in enumerate(local_srclist):
                    known.logs[0].info('processing image %s'%se_info['red_image'])
                    file_id=local_file_ids[i_se]
                    se_info['file_id']=file_id
                    se_setup = SEBalrogSetup(coaddsetup_base, se_info, funpack=False)
                    se_setup.bigImage,se_setup.weight,se_setup.badpix,se_setup.sky,se_setup.psfmodel,se_setup.wcs = se_setup.ReadImages()
                    if known.no_noise:
                        se_setup.bigImage*=0
                    
                    i_cat=0
                    for i_sim in range(known.n_sims):
                        for i_shear in range(known.n_shears):
                            catalog = catalogs_allsims[i_cat][i_band]
                            #sex_data = sex_data_allbands_allsims[i_sim][band]
                            sim_se_objects = SERun(catalog, se_info, se_setup)
                            comm.send(sim_se_objects, dest=master, tag=i_cat)
                            i_cat+=1
            else:
                sim_se_objects=[]
                se_returned=0
                while se_returned<len(coadd.srclist*len(catalogs_allsims)):
                    status=MPI.Status()                
                    rec_data=comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                    i_cat=status.Get_tag()
                    for r in rec_data:
                        sim_se_objects.append((i_cat,r[0],r[1]))
                    se_returned+=1
                    known.logs[0].info('se_returned: %d/%d'%(se_returned,len(coadd.srclist*known.n_sims)))

                for se_obj in sim_se_objects:
                    i_cat,i_obj = se_obj[0], se_obj[1]
                    meds_inputs_allbands_allsims[i_cat][band][i_obj].add_SEObject(se_obj[2])

                #prepend coadd info to source list
                srclist=[{'id':coadd['image_id'], 'red_image':coadd['image_url']}]+coadd.srclist
                meds_files=[]
                for i_cat in range(len(catalogs_allsims)):
                    i_sim,i_shear=sim_nums[i_cat],shear_nums[i_cat]
                    MEobjs=[]
                    for i, meds_input in enumerate(meds_inputs_allbands_allsims[i_cat][band]):
                        MEobjs.append(meds_input.make_MEobj(dummy_segs=False))
                    meds_name = os.path.join(known.parent_output, "%s-%s-meds.%d.%d.fits")%(known.tilename,band,i_sim,i_shear)
                    meds_files.append(meds_name)
                    print MEobjs
                    des_meds.write_meds(meds_name, MEobjs, srclist=srclist, clobber=True)
                    known.logs[0].info("wrote %d objects to meds file %s"%(len(meds_inputs_allbands_allsims[i_cat][band]),meds_name))
                #Save catalogs
                with open(os.path.join(known.parent_output,'cat_list.%s.txt'%band),'w') as f:
                    for truth, measured, meds in zip(truthcats,measuredcats,meds_files):
                        outline="%s\t%s\t%s\n"%(truth[i_band],measured[i_band],meds)
                        f.write(outline)
                    
            print comm.Get_rank()
            comm.Barrier()


def TileArgs(parser):
    parser.add_argument("tilename",type=str)
    parser.add_argument("parent_output",type=str,help="the output directory (each simulation will be saved to a different directory within this one)")
    parser.add_argument("--bands",nargs='*',type=str,default=['i'])
    parser.add_argument("--detbands",nargs='*',type=str,default=['i'])
    parser.add_argument("--no_meds",action='store_true',default=False)
    parser.add_argument("--no_noise",action='store_true',default=False,help="make 'noiseless' images i.e. set background image to zero everywhere")
    parser.add_argument("--n_sims",type=int,default=1, help="number of simulations (of same tile)")
    parser.add_argument("--pickle_seobjects",action='store_true',help="Save each SE image's list SEObjects as a pickle, then have master process read them in and construct meds file")
    parser.add_argument("--single_meds",action='store_true',help="Save all sims to a single meds file")
    parser.add_argument("--no_mpi", action='store_true', default=False)
    parser.add_argument("--coadd_nproc", default=None)
    parser.add_argument("--shear_perturb", default=None, type=float, nargs=2, help="run each simulation twice, with same catalog except perturbed by a small shear")

def main(parser, known):

    #Set up 'parent_output' - directory where output from all processes/sims goes
    #Within this set up and output directory for each process, containing a log directory
    try:
        os.makedirs(known.parent_output)
    except OSError:
        if not os.path.isdir(known.parent_output):
            raise

    comm = None
    if not known.no_mpi:
        #from mpi4py import MPI
        #print 'importing mpi4py'
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        known.outdir = os.path.join(known.parent_output,known.tilename+"-rank.%d"%rank)
    else:
        known.outdir = os.path.join(known.parent_output,known.tilename)
    if not os.path.isdir(known.outdir):
        os.mkdir(known.outdir)
    known.logdir = os.path.join(known.outdir, 'balrog_log')
    CreateSubDir(known.logdir)
    known.logs = SetupLogger(known)

    #Run the main method
    try:
        Run(parser, known, comm=comm)
    except:
        print known.logs[0]
        RaiseException(known.logs[0], fulltraceback=known.fulltraceback)


if __name__ == "__main__":
   
    parser = GetNativeOptions()
    known = GetKnown(parser)
    main(parser,known)
