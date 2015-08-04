from balrog import * #DefaultArgs, CreateSubDir, SetupLogger, RaiseException
import pywcs
import os
import pyfits
import argparse
import sys
import shutil
import astropy.io.fits as astropyfits
import multi_epoch_tools.coords as coords
import multi_epoch_tools.file_tools as file_tools
import multi_epoch_tools.des_meds as des_meds
import pylab

MAGZP_REF=30.

#Read in tile name and balrog config file
#Generate x,y and use coadd wcs to convert to ra,dec 

class SEBalrogSetup(dict):
    def __init__(self,CoaddSetup,srclist_entry,funpack=True):
        #for now copy all coadd settings to get defaults correctly
        for key in CoaddSetup.keys():
            self[key]=CoaddSetup[key]
        #Now overwrite the relevant ones
        #For now save single epoch images in copy of red directory
        imgdir=os.path.dirname(srclist_entry['red_image'])+"_balrog_image"
        img_path=srclist_entry['red_image']
        if funpack:
            #funpack dat shit
            file_tools.funpack(srclist_entry['red_image'])
            img_path=img_path.strip(".fz")
            img_suffix=".fits"
            self['imageext']=0
            self['weigtext']=1
        else:
            img_suffix='.fits.fz'
            self['imageext']=1
            self['weightext']=2

        img_basename=os.path.basename(img_path)
        print 'img_basename',img_basename
        if not os.path.isdir(imgdir):
            os.mkdir(imgdir)
        self['logdir']=os.path.dirname(srclist_entry['red_image'])+"_balrog_log"
        self["logs"] = SetupLogger(known)
        self.runlogger = self["logs"][0]
        self.sexlogger = self["logs"][1]
        self.arglogger = self["logs"][2]

        if not os.path.isdir(self['logdir']):
            os.mkdir(self['logdir'])
        catdir=os.path.dirname(srclist_entry['red_image'])+"_balrog_cat"
        if not os.path.isdir(catdir):
            os.mkdir(catdir)
        sexconfigdir=os.path.dirname(srclist_entry['red_image'])+"_sexconfig"
        if not os.path.isdir(sexconfigdir):
            os.mkdir(sexconfigdir)
        self['detimageext']=self['imageext']
        self['imgdir']=imgdir
        self['sexlog']=os.path.join(self['logdir'],img_basename.replace(img_suffix,'.sexlog'))
        self['zeropoint']=srclist_entry['magzp']
        self['weight']=img_path
        self['image']=img_path
        self['runlog']=os.path.join(self['logdir'],img_basename.replace(img_suffix,'.runlog'))
        #self['weightout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.nosim.fits'))
        self['imageout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.sim.fits'))
        self['detimageout']=self['imageout']
        self['psfout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.psf'))
        self['detpsf']= img_path.replace(img_suffix,'_psfcat.psf')
        #self['balrog_output_dir']= test
        #self['noassoc']=False
        self['catalogmeasured']=os.path.join(catdir,img_basename.replace(img_suffix,'.measuredcat.sim.fits'))
        self['detimagein']=img_path
        self['catalogtruth']=os.path.join(catdir,img_basename.replace(img_suffix,'.truthcat.sim.fits'))
        self['assoc_simfile']=os.path.join(sexconfigdir,img_basename.replace(img_suffix,'.assoc.sim.txt'))
        self['detweight']=img_path
        self['assoc_nosimfile']=os.path.join(sexconfigdir,img_basename.replace(img_suffix,'.assoc.nosim.txt'))
        self['psf']= img_path.replace(img_suffix,'_psfcat.psf')
        #catalog None
        #catext 2
        #sersicindex SERSIC_INDEX
        self['gain']=srclist_entry['gaina']
        self['sexdir']= sexconfigdir
        #self['outdir']= /home/maccrann/code/Balrog/default_example/output
        #sexnnw /home/maccrann/code/Balrog/astro_config/sex.nnw
        #detweightext 1
        #self['detweightout']= os.path.join(imgdir,img_basename.replace(img_suffix,'.sim.fits'))
        #sexlogger <logging.Logger object at 0x2b8c050>
        self['scampheader']=srclist_entry['astro_refine']
        #logverbosity n
        self['detimage']= self['image']
        self['detpsfout']= self['psfout']
        self['nosim_detimageout']= self['image']
        #sexconv /home/maccrann/code/Balrog/astro_config/sex.conv
        #ext 1
        self['nosim_catalogmeasured']=os.path.join(catdir,img_basename.replace(img_suffix,'.measuredcat.nosim.fits'))
        #nosimsexparam /home/maccrann/code/Balrog/astro_config/sex.param
        #clean False
        #weightext 1
        #nodraw False
        self['nosim_imageout']=os.path.join(imgdir,img_basename.replace(img_suffix,'.nosim.fits'))
        self['arglog']=os.path.join(self['logdir'],img_basename.replace(img_suffix,'.arglog'))
        self['xmax']=pyfits.open(self['image'])[self['imageext']].header['NAXIS1']
        self['ymax']=pyfits.open(self['image'])[self['imageext']].header['NAXIS2']
        for key,val in self.iteritems():
            setattr(self,key,val)

        print self.sexconfig, self.sexdir
        shutil.copy(self.sexconfig, self.sexdir)
        shutil.copy(self.sexparam, self.sexdir)
        shutil.copy(self.sexnnw, self.sexdir)
        shutil.copy(self.sexconv, self.sexdir)
        shutil.copy(self.nosimsexparam, self.sexdir)
        if os.path.lexists(self.pyconfig):
            shutil.copy(self.pyconfig, self.logdir)

        if self['catfitstype']=='FITS_LDAC':
            self['catext'] = 2
        elif self['catfitstype']=='FITS_1.0':
            self['catext'] = 1

        #Set keys as attributes:
        for key,val in self.iteritems():
            setattr(self,key,val)
        print self.arglog

    def ReadImages(self):
        _,self.input_hdulist,_=galsim.fits.readFile(self.image)
        print self.input_hdulist
        image = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.imageext)#,compression='rice')
        weight = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.weightext)#,compression='rice') 
        badpix = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.weightext+1)
        if self.scampheader is not None:
            h=pyfits.Header()
            h.fromTxtFile(self.scampheader)
            wcs=galsim.wcs.readFromFitsHeader(h)[0]
        else:
            wcs=image.wcs
        print image.wcs
        image.wcs = wcs
        print image.wcs
        weight.wcs = image.wcs

        subBounds = galsim.BoundsI(self.xmin,self.xmax,self.ymin,self.ymax)
        image = image[subBounds]
        weight = weight[subBounds]
        badpix = badpix[subBounds]
        psfmodel = galsim.des.DES_PSFEx(self.psf, wcs=wcs)
        return image, weight, badpix, psfmodel, wcs

    def WriteImages(self, image, nosim=False):
        print 'nosim',nosim
        if nosim:
            imageout = self.nosim_imageout
        else:
            imageout = self.imageout

        #Replace image hdu data
        self.input_hdulist[self.imageext].data=np.zeros_like(image.array)
        print self.input_hdulist
        #Write to imageout
        print  imageout
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
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    known.outdir = known.tilename+"_output"
    known.logdir = os.path.join(known.outdir, 'balrog_log')
    CreateSubDir(known.logdir)
    known.logs = SetupLogger(known)
    return known

def TileArgs(parser):
    parser.add_argument("tilename",type=str)
    parser.add_argument("balrog_output_dir",type=str)
    parser.add_argument("--band",type=str,default='i')
    parser.add_argument("--sync_coadd",default=False)
    parser.add_argument("--make_meds",action='store_true',default=True)

def NativeParse(parser, known):
    cmdline_opts = parser.parse_args()
    known.logs[2].info('# Exact command call')
    known.logs[2].info(' '.join(sys.argv))
    LogCmdlineOpts(cmdline_opts, cmdline_opts, known.logs[2], '\n# Values received for each possible command line option, filling with defaults if necessary')
    
    #ParseDefaultArgs(cmdline_opts, known)
    SetupCoadd(cmdline_opts,known.logs[0], known)
    #Get derived args - output files mostly
    #dargs=DerivedArgs(cmdline_opts,known)
    BalrogSetup = ConfigureBalrog(cmdline_opts, known)
    return cmdline_opts, BalrogSetup

def SetupCoadd(args, log, known):
    args.pyconfig = known.pyconfig
    args.outdir = known.outdir
    thisdir = os.path.dirname( os.path.realpath(__file__) )

    #Setup coadd image
    args.image = known.coadd['image_url_funpacked']
    args.detimage = args.image
    args.imageext = 0
    args.detimageext = 0
    args.weight = args.image
    args.weightext = 1
    args.detweight = args.image
    args.detweightext = 1
    args.psf = known.coadd['psf_path']
    args.detpsf = args.psf

    #Setup logs here
    known.logdir = os.path.join(known.outdir, 'balrog_log')
    CreateSubDir(known.logdir)
    known.logs = SetupLogger(known)

    attrs, ss, eattrs, ess, codes = get_check_array_ext(args)
    for attr,s, eattr,es, code in zip(attrs,ss, eattrs,ess, codes):
        images_existence_check(attr, code, s)
        images_are_fits(attr, code+10, s)
        #if s!='psf':
        exts_exist(attr, s, eattr, es, code+20)

    SizesOK(args, log)
    
    #Same as normal Balrog
    configdir = os.path.join(thisdir, 'astro_config')
    ParseFloatKeyword(args, known.logs[0]) 
    ParseSex(args, known.logs[0], configdir)


def ParseDefaultArgs(args,known):
    
    ParseImages(args, known.logs[0], indir)
    ParseFloatKeyword(args, known.logs[0]) 
    ParseSex(args, known.logs[0], configdir)

    return args

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

class MedsInput(dict):
    def __init__(self,id=0):
        self['id']=id
        self['images']=[]
        self['weights']=[]
        self['badpixs']=[]
        self['segs']=[]
        self['wcss']=[]
        self['orig_rows']=[]
        self['orig_cols']=[]
        self['orig_start_rows']=[]
        self['orig_start_cols']=[]
        self['cutout_rows']=[]
        self['cutout_cols']=[]

    def add(self,image,weight,badpix,seg,wcs,orig_row,orig_col,orig_start_row,orig_start_col,cutout_row,cutout_col):
        self['images'].append(image)
        self['weights'].append(weight)
        self['badpixs'].append(badpix)
        self['segs'].append(seg)
        self['wcss'].append(wcs)
        self['orig_rows'].append(orig_row)
        self['orig_cols'].append(orig_col)
        self['orig_start_rows'].append(orig_start_row)
        self['orig_start_cols'].append(orig_start_col)
        self['cutout_rows'].append(cutout_row)
        self['cutout_cols'].append(cutout_col)

    def make_MEobj(self):
        return des_meds.MultiExposureObject(images=self['images'], weights=self['weights'], badpix=self['badpixs'], segs=self['segs'], 
                                            wcs=self['wcss'], id=self['id'], number=0, orig_rows=self['orig_rows'],
                                            orig_cols=self['orig_cols'], orig_start_rows=self['orig_start_rows'], 
                                            orig_start_cols=self['orig_start_cols'], cutout_rows=self['cutout_rows'], 
                                            cutout_cols=self['cutout_cols'])

def get_stamp(image,x_pos,y_pos,stamp_size,offset=(0,0),scale=1.):
    "Get postage stamp to run im3shape on"
    "Expects galsim image, SExtractor X_IMAGE,Y_IMAGE and stamp size as input"
    "Returns stamp, x0 and y0 - coordinates of object in stamp (starting from 0)"
    "and galsim bounds object"
    half=stamp_size/2
    x_min,x_max,y_min,y_max = int(x_pos+offset[0]-half),int(x_pos++offset[0]+half),int(y_pos+offset[1]-half),int(y_pos+offset[1]+half)
    #over_edge = (x_min<0) | (y_min<0) | (x_max>image.array.shape[1]) | (y_max>image.array.shape[0])
    subBounds = galsim.BoundsI(x_min,x_max-1,y_min,y_max-1)
    stamp = image[subBounds] * scale
    x0,y0 = x_pos - x_min, y_pos - y_min
    return stamp, x0, y0, subBounds
        
def RunBalrog(parser, known):

    boxsize=64
    #First step is to read in the coadd and generate the truth catalog based on this:

    # Find the user's config file
    config = GetConfig(known)

    # Add the user's command line options
    AddCustomOptions(parser, config, known.logs[0])

    #Setup coadd object to get coadd path etc. Use det image as default
    coadd=file_tools.setup_tile(known.tilename, band=known.band, sync=known.sync_coadd)
    coadd.funpack()
    known.coadd = coadd
    #print coadd

    # Parse the command line arguments and interpret the user's settings for the simulation
    cmdline_opts, CoaddSetup = NativeParse(parser, known)
    rules, extra_sex_config, cmdline_opts_copy, TruthCatExtra = CustomParse(cmdline_opts, CoaddSetup, config)

    # Take the the user's configurations and build the simulated truth catalog out of them.
    catalog, gspcatalog, extracatalog, TruthCatExtra = GetSimulatedGalaxies(CoaddSetup, rules, config, cmdline_opts_copy, TruthCatExtra)

    if known.make_meds:
        meds_inputs=[]
        for i in range(known.ngal):
            meds_inputs.append(MedsInput(id=i))

    # Insert galaxies into coadd. For now just do one band. Next step would be to insert into all detection bands
    print CoaddSetup
    bigImage, subWeight, psfmodel, wcs = ReadImages(CoaddSetup)

    if not CoaddSetup.imageonly:
        # If desired, run sextractor over the image prior to inserting any simulated galaxies.
        if not CoaddSetup.nonosim:
            NosimRunSextractor(CoaddSetup, bigImage, subWeight, extra_sex_config, catalog)
    # Insert simulated galaxies.
    if not CoaddSetup.nodraw:
        bigImage,wcs_list_coadd = InsertSimulatedGalaxies(bigImage, catalog, psfmodel, CoaddSetup, wcs, 
                                                    gspcatalog, return_wcs_list=True)
        #WriteImages(CoaddSetup, bigImage, subWeight)
        WriteCatalog(catalog, CoaddSetup, txt=None, fits=True, TruthCatExtra=TruthCatExtra, extracatalog=extracatalog)

    #If making meds on the fly, need to get image, weight and badpix cutouts
    for i,(x,y) in enumerate(zip(catalog.galaxy['x'],catalog.galaxy['y'])):
        print i
        if catalog.galaxy['not_drawn'][i]:
            continue
        image_stamp,x0,y0,bounds=get_stamp(bigImage,x,y,boxsize)
        print 'type(image)',type(image_stamp)
        weight_stamp,_,_,_=get_stamp(subWeight,x,y,boxsize)
        seg_stamp=badpix_stamp=galsim.ImageI(np.zeros_like(image_stamp.array))
        meds_inputs[i].add(image_stamp,weight_stamp,badpix_stamp,seg_stamp,wcs_list_coadd[i],y-1,x-1,bounds.ymin-1,bounds.xmin-1,y0,x0)
        pylab.imshow(image_stamp.array,interpolation='nearest',origin='lower')
        pylab.show()


    """
    else:
        WriteImages(CoaddSetup, bigImage, subWeight)    
        
    if not CoaddSetup.imageonly:

        # Run sextractor over the simulated image.
        RunSextractor(CoaddSetup, extra_sex_config, catalog)

    # If chosen, clean up image files you don't need anymore
    if CoaddSetup.clean:
        Cleanup(CoaddSetup)

    # Log some  extra stuff Balrog used along the way
    LogDerivedOpts(cmdline_opts, CoaddSetup, '\n#Psuedo-args. Other values derived from the command line arguments.')
    """
    #Print out the runs needed to get all the single exposures for this tile
    runs=[]
    for s in coadd.srclist:
        runs.append(s['run'])
    print np.unique(runs)

    # get galaxy ra,dec
    header=pyfits.getheader(coadd['image_url_funpacked'])
    coadd_x,coadd_y=catalog.galaxy['x'],catalog.galaxy['y']
    catalog.galaxy['ra'],catalog.galaxy['dec']=coords.get_wcoords(coadd_x,coadd_y,header)
    """
    print 'coadd_x, coadd_y:', coadd_x,coadd_y
    print 'ra, dec:',catalog.galaxy['ra'],catalog.galaxy['dec']
    """
    #Loop through single epoch images, inserting galaxies
    #If in meds mode, we also want to add image, weight etc. postage stamp to the appropriate lists, with which a
    #meds object will be created afterwards
    
    for i,se_info in enumerate(coadd.srclist[0:15]):
        print i,se_info
        BalrogSetup=SEBalrogSetup(CoaddSetup,se_info)
        print BalrogSetup['zeropoint']
        #print '**************************'
        #for key,val in BalrogSetup.iteritems():
        #    print key,val

        # Get the subsampled flux and weightmap images, along with the PSF model and WCS.
        bigImage, subWeight, badpix, psfmodel, wcs = BalrogSetup.ReadImages()
        #print bigImage.array.shape
        bigImage=galsim.ImageD(np.zeros_like(bigImage.array))
        # Get x and y for this exposure
        x,y=coords.get_pcoords(catalog.galaxy['ra'],catalog.galaxy['dec'],wcs)
        catalog.galaxy['x'],catalog.galaxy['y']=x,y
        print catalog.galaxy['x'],catalog.galaxy['y']
        """
        if not BalrogSetup.imageonly:
            # If desired, run sextractor over the image prior to inserting any simulated galaxies.
            if not BalrogSetup.nonosim:
                NosimRunSextractor(BalrogSetup, bigImage, subWeight, extra_sex_config, catalog)
        """
        # Insert simulated galaxies.
        if not BalrogSetup.nodraw:
            bigImage,wcs_list = InsertSimulatedGalaxies(bigImage, catalog, psfmodel, BalrogSetup, wcs, gspcatalog, return_wcs_list=True)
            BalrogSetup.WriteImages(bigImage)
            WriteCatalog(catalog, BalrogSetup, txt=None, fits=True, TruthCatExtra=TruthCatExtra, extracatalog=extracatalog)
        else:
            BalrogSetup.WriteImages(bigImage, subWeight)
        
        #If making meds on the fly, need to get image, weight and badpix cutouts
        #First calculate scaling factor to set to magzp_ref
        scale=10**(0.4*(MAGZP_REF-se_info['magzp']))
        for i,(x,y) in enumerate(zip(catalog.galaxy['x'],catalog.galaxy['y'])):
            print i
            if catalog.galaxy['not_drawn'][i]:
                continue
            image_stamp,x0,y0,bounds=get_stamp(bigImage,x,y,boxsize,offset=(5,0),scale=scale)
            print 'type(image)',type(image_stamp)
            weight_stamp,_,_,_=get_stamp(subWeight,x,y,boxsize)
            badpix_stamp,_,_,_=get_stamp(badpix,x,y,boxsize)
            seg_stamp=galsim.ImageI(np.zeros_like(image_stamp.array))
            meds_inputs[i].add(image_stamp,weight_stamp,badpix_stamp,seg_stamp,wcs_list[i],y-1,x-1,bounds.ymin-1,bounds.xmin-1,y0,x0)
            pylab.imshow(image_stamp.array,interpolation='nearest',origin='lower')
            pylab.show()


        # If chosen, clean up image files you don't need anymore
        if BalrogSetup.clean:
            Cleanup(BalrogSetup)

        # Log some  extra stuff Balrog used along the way
        #LogDerivedOpts(cmdline_opts, BalrogSetup, '\n#Psuedo-args. Other values derived from the command line arguments.')
    MEobjs=[]
    for i,meds_input in enumerate(meds_inputs):
        if len(meds_input['images'])<1: 
            print "no exposures for galaxy",i
            continue
        print "type(meds_input['images'][0])",type(meds_input['images'][0])
        MEobjs.append(meds_input.make_MEobj())
        des_meds.write_meds('meds_test', MEobjs, clobber=True)


    '''
    del bigImage
    del subWeight

    import gc
    gc.collect()
    '''

if __name__ == "__main__":
   
    # First get the needed info to setup up the logger, which allows everything to be logged even if things fail at very early stages.
    # Only a writable outdir is required to be able to get the output log file.

    #import resource
    parser = GetNativeOptions()
    known = GetKnown(parser)
    print parser
    print known
    try:
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000
        RunBalrog(parser, known)
        #print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000
    except:
        print known.logs[0]
        RaiseException(known.logs[0], fulltraceback=known.fulltraceback)

