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
import pylab
import psutil
import functools
import copy

MAGZP_REF=30
BALROG_DESDATA="/home/maccrann/balrog_desdata"

class SEBalrogSetup(dict):
    def __init__(self,CoaddSetup,srclist_entry,funpack=True):
        #for now copy all coadd settings to get defaults
        for key in CoaddSetup.keys():
            self[key]=CoaddSetup[key]

        #Now overwrite the relevant ones
        img_path=srclist_entry['red_image']
        sky_path=srclist_entry['red_bkg']
        if funpack:
            img_path=img_path.strip(".fz")
            if not os.path.isfile(img_path):
                file_tools.funpack(srclist_entry['red_image'])
                file_tools.funpack(srclist_entry['red_bkg'])
            self['sky_path']=sky_path.strip(".fz")            
            img_suffix=".fits"
            self['skyext']=0
            self['imageext']=0
            self['weightext']=2
        else:
            img_suffix='.fits.fz'
            self['imageext']=1
            self['weightext']=3

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
        self['xmax']=pyfits.open(self['image'])[self['imageext']].header['NAXIS1']
        self['ymax']=pyfits.open(self['image'])[self['imageext']].header['NAXIS2']
        for key,val in self.iteritems():
            setattr(self,key,val)
        if self['catfitstype']=='FITS_LDAC':
            self['catext'] = 2
        elif self['catfitstype']=='FITS_1.0':
            self['catext'] = 1

        #Set keys as attributes:
        for key,val in self.iteritems():
            setattr(self,key,val)

    def ReadImages(self):
        _,self.input_hdulist,_=galsim.fits.readFile(self.image)
        image = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.imageext)#,compression='rice')
        weight = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.weightext)#,compression='rice') 
        badpix = galsim.fits.read(hdu_list=self.input_hdulist,hdu=self.imageext+1)
        sky = galsim.fits.read(file_name=self.sky_path, hdu=self.skyext)
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
    SetupCoadd(cmdline_opts,known.logs[0], known)
    BalrogSetup = ConfigureBalrog(cmdline_opts, known)
    return cmdline_opts, BalrogSetup

def SetupCoadd(args, log, known):
    args.pyconfig = known.pyconfig
    args.outdir = known.outdir
    thisdir = os.path.dirname( os.path.realpath(__file__) )

    #Setup coadd image
    try:
        args.image = known.coadd['image_url_funpacked']
    except KeyError:
        known.coadd['image_url_funpacked']=(known.coadd['image_url']).rsplit('.',1)[0]
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

    #Also set path to segmentation map
    imgdir = os.path.dirname(known.coadd['image_url'])
    imgbase = os.path.basename(known.coadd['image_url'])
    args.seg = os.path.join(os.path.dirname(imgdir),'QA','segmap',imgbase.replace('.fits.fz','_seg.fits.fz'))

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
        exts_exist(attr, s, eattr, es, code+20)

    SizesOK(args, log)

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

def catalog_to_table(sample, BalrogSetup, TruthCatExtra=None, extracatalog=None):
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
        elif key=='ra':
            unit = 'deg'
        elif key=='dec':
            unit = 'deg'
        else:
            unit = 'dimensionless'
        col = astropy.table.Column(data=arr, name=name,format='E', unit=unit, dtype=float)
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
                
            col = astropy.table.Column(data=extracatalog.galaxy[name], name=name, format=fmt, unit=unit)
            columns.append(col)


    for i in range(len(sample.component)):
        for key in sample.component[i].keys():
            name = '%s_%i' %(key,i)
            if key.find('halflightradius')!=-1:
                col = astropy.table.Column(name=name, data=sample.component[i][key]/np.sqrt(sample.component[i]['axisratio']), format='E', unit='arcsec')
            else:
                if key.find('sersicindex')!=-1:
                    unit = 'dimensionless'
                if key.find('flux')!=-1:
                    unit = 'ADU'
                if key.find('beta')!=-1:
                    unit = 'deg'
                if key.find('axisratio')!=-1:
                    unit = 'dimensionless'
                col = astropy.table.Column(name=name, data=sample.component[i][key],format='E', unit=unit)
            columns.append(col)

    return astropy.table.Table(data=columns)

def CoaddRun(CoaddSetup):

    # Take the the user's configurations and build the simulated truth catalog out of them.
    catalog, gspcatalog, extracatalog, TruthCatExtra = GetSimulatedGalaxies(CoaddSetup, CoaddSetup.rules, 
                                                                            CoaddSetup.config, CoaddSetup.cmdline_opts_copy, 
                                                                            CoaddSetup.TruthCatExtra)
    coadd_x,coadd_y=catalog.galaxy['x'],catalog.galaxy['y']
    catalog.galaxy['ra'],catalog.galaxy['dec']=coords.get_wcoords(coadd_x,coadd_y,CoaddSetup.wcs)

    meds_inputs=[]
    for i in range(CoaddSetup.ngal):
        meds_inputs.append(des_meds.MedsInput(id=CoaddSetup.indexstart+i))


    # If desired, run sextractor over the image prior to inserting any simulated galaxies.
    NosimRunSextractor(CoaddSetup, CoaddSetup.bigImage, CoaddSetup.subWeight, 
                       CoaddSetup.extra_sex_config, catalog)

    # Insert simulated galaxies.
    bigImage,wcs_list_coadd = InsertSimulatedGalaxies(CoaddSetup.bigImage, catalog, CoaddSetup.psfmodel, 
                                                      CoaddSetup, CoaddSetup.wcs, 
                                                      gspcatalog, return_wcs_list=True)
    
    CoaddSetup.runlogger.info("Inserted Simulated galaxies")
    WriteImages(CoaddSetup, CoaddSetup.bigImage, CoaddSetup.subWeight)
    WriteCatalog(catalog, CoaddSetup, txt=None, fits=True, TruthCatExtra=CoaddSetup.TruthCatExtra, extracatalog=extracatalog)
    truth_data=catalog_to_table(catalog, CoaddSetup, TruthCatExtra=CoaddSetup.TruthCatExtra, extracatalog=extracatalog)
    RunSextractor(CoaddSetup, CoaddSetup.extra_sex_config, catalog)            

    #Get new combined seg
    seg, blended_with_orig=combined_seg(CoaddSetup)

    sex_data=astropy.table.Table(pyfits.getdata(CoaddSetup.catalogmeasured, CoaddSetup.catext))
    #replace vector_assoc with balrog_index
    balrog_index_det=sex_data['VECTOR_ASSOC'][:,CoaddSetup.assocnames.index('balrog_index')].astype(int)
    col=astropy.table.Column(balrog_index_det, name='balrog_index')
    sex_data.remove_column('VECTOR_ASSOC')
    sex_data.add_column(col)
    CoaddSetup.runlogger.info("Saved images and ran SExtractor")

    #Make new column specifying whether object is blended with original object (using seg map)
    blended_with_orig=list(blended_with_orig)
    new_sex_ids=np.array(sex_data['NUMBER'])
    col_data=np.zeros_like(new_sex_ids)
    for i,new_id in enumerate(new_sex_ids):
        if new_id in blended_with_orig:
            col_data[i]=1        
    blended_with_orig_col=astropy.table.Column(col_data, name='blended_with_orig')
    sex_data.add_column(blended_with_orig_col)

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
    if not CoaddSetup.no_meds:
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
                image_stamp,x0,y0,bounds=get_stamp(CoaddSetup.bigImage,x,y,box_size[i])
            except OutOfBoundsError:
                continue
            weight_stamp,_,_,_=get_stamp(CoaddSetup.subWeight,x,y,box_size[i])
            seg_stamp,_,_,_=get_stamp(seg,x,y,box_size[i])
            badpix_stamp=galsim.ImageI(np.zeros_like(image_stamp.array))
            meds_inputs[i].add(image_stamp,weight_stamp,badpix_stamp,seg_stamp,wcs_list_coadd[i],y-1,x-1,bounds.ymin-1,bounds.xmin-1,y0,x0,0)

    coadd_not_drawn=np.copy(catalog.galaxy['not_drawn'])
    CoaddSetup.runlogger.info("%d objects not drawn"%(coadd_not_drawn).sum())
    Cleanup(CoaddSetup)
    
    try:
        shutil.copytree(CoaddSetup.outdir,CoaddSetup.sim_outdir)
    except OSError:
        shutil.rmtree(CoaddSetup.sim_outdir)
        shutil.copytree(CoaddSetup.outdir,CoaddSetup.sim_outdir)
    CoaddSetup.runlogger.info("Copied output files to sim dir")
    
    return catalog, sex_data, truth_data, meds_inputs, gspcatalog

def combined_seg(BalrogSetup):
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
    seg=galsim.fits.read(BalrogSetup.seg,hdu=1)
    #Read new seg map
    seg_new=galsim.fits.read(BalrogSetup.extra_sex_config['CHECKIMAGE_NAME'])
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


def SERun(catalog, se_info, BalrogSetup):

    x,y=coords.get_pcoords(catalog.galaxy['ra'],catalog.galaxy['dec'],BalrogSetup.wcs)
    catalog.galaxy['x'],catalog.galaxy['y']=x,y
    """
    if not BalrogSetup.imageonly:
        # If desired, run sextractor over the image prior to inserting any simulated galaxies.
        if not BalrogSetup.nonosim:
            NosimRunSextractor(BalrogSetup, bigImage, subWeight, extra_sex_config, catalog)
    """
    # Insert simulated galaxies.
    bigImage, wcs_list = InsertSimulatedGalaxies(BalrogSetup.bigImage, 
                                                 catalog, BalrogSetup.psfmodel, 
                                                 BalrogSetup, BalrogSetup.wcs, BalrogSetup.gspcatalog, 
                                                 return_wcs_list=True, out_of_bounds_threshold=48)
    """
    Would save SE image, run SExtractor, and save catalogs here
    BalrogSetup.WriteImages(bigImage)
    WriteCatalog(catalog, BalrogSetup, txt=None, fits=True, TruthCatExtra=TruthCatExtra, extracatalog=extracatalog)
    """

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
        badpix_stamp,_,_,_=get_stamp(BalrogSetup.badpix,x,y,catalog.galaxy['box_size'][i_obj])
        sky_stamp,_,_,_=get_stamp(BalrogSetup.sky,x,y,catalog.galaxy['box_size'][i_obj])
        #sky subtract, rescale image stamp and set badpix to zero
        if not known.no_noise:
            image_stamp-=sky_stamp
        image_stamp*=scale
        weight_stamp/=(scale*scale)
        image_stamp_array=image_stamp.array
        badpix_array=badpix_stamp.array
        image_stamp_array[badpix_array!=0]=0.
        image_stamp=galsim.ImageD(image_stamp_array)
        seg_stamp=None
        se_object=des_meds.SEObject(image_stamp,weight_stamp,badpix_stamp,seg_stamp,wcs_list[i_obj],y-1,x-1,
                               bounds.ymin-1,bounds.xmin-1,y0,x0,se_info['file_id'])
        se_objects.append((i_obj, se_object))

    return se_objects


def Run(parser,known, comm=None):

    #Get number of processes. If no mpi, set rank to 0 and everything else *should* work
    if comm is not None:
        rank = comm.Get_rank()
        size = comm.Get_size()
    else:
        rank, size = 0,1

    #First step is to read in the coadd and generate the truth catalog based on this:
    # Find the user's config file
    config = GetConfig(known)

    # Add the user's command line options
    AddCustomOptions(parser, config, known.logs[0])

    #Setup coadd object to get coadd path etc.
    coadd=file_tools.setup_tile(known.tilename, band=known.band, sync=False)
    known.coadd = coadd
    known.coadd.funpack()

    if rank==0:
        print 'coadd:',coadd

    # Parse the command line arguments and interpret the user's settings for the simulation
    cmdline_opts, CoaddSetup = NativeParse(parser, known)
    CoaddSetup.config = config
    CoaddSetup.rules, CoaddSetup.extra_sex_config, CoaddSetup.cmdline_opts_copy, CoaddSetup.TruthCatExtra = CustomParse(cmdline_opts, CoaddSetup, config)

    # Add sex options to save seg mask
    CoaddSetup.segout=CoaddSetup.imageout.replace('.fits','_seg.fits')
    CoaddSetup.extra_sex_config['CHECKIMAGE_TYPE']='SEGMENTATION'
    CoaddSetup.extra_sex_config['CHECKIMAGE_NAME']=CoaddSetup.segout

    #Calculate how many coadd sims each process has to do. If none, set idle to True
    idle=False
    if known.coadd_nproc is None or known.coadd_nproc>size:
        coadd_size = size
    else:
        coadd_size = known.coadd_nproc
    sim_inds=np.arange(known.n_sims)
    local_sim_inds=np.array_split(sim_inds,coadd_size)[rank]
    CoaddSetup.runlogger.info('process %d doing %d sims'%(rank,len(local_sim_inds)))
    if len(local_sim_inds) < 1:
        print 'process %d has no sims assigned, nothing to do'%rank
        comm.send(None)
        idle=True

    catalogs_allsims, gspcatalog = None, None

    #If not idle, read in images and balrog (using the 'CoaddRun' method) the given number of times
    if not idle:
        local_catalogs=[]
        local_meds_inputs=[]
        local_truthcats=[]
        local_sexcats=[]
        #Read in images
        CoaddSetup.runlogger.info("Reading orignal image/psf etc.")
        CoaddSetup.bigImage, CoaddSetup.subWeight, CoaddSetup.psfmodel, CoaddSetup.wcs = ReadImages(CoaddSetup)
        if known.no_noise:
            CoaddSetup.bigImage *= 0
            CoaddSetup.runlogger.info("Set image to zero (no_noise mode)")
        CoaddSetup.runlogger.info("Done reading image etc.")

        for i_sim in local_sim_inds:
            CoaddSetup.sim_outdir=os.path.join(CoaddSetup.parent_output,CoaddSetup.tilename+".%d"%i_sim)
            CoaddSetup.runlogger.info("Doing sim %d"%i_sim)
            #Want each object in each sim to have a unique balrog index (for matching later)
            #So start at i_sim * CoaddSetup.ngal
            CoaddSetup.indexstart = i_sim * known.ngal
            catalog, sex_data, truth_data, meds_inputs, gspcatalog = CoaddRun(CoaddSetup)
            local_catalogs.append(copy.copy(catalog))
            local_meds_inputs.append(meds_inputs)
            local_truthcats.append(truth_data)
            local_sexcats.append(sex_data)

        CoaddSetup.bigImage,CoaddSetup.subWeight=None,None
        #rank 0 collects output from all sims
        if rank==0:
            catalogs_allsims, meds_inputs_allsims, truthcats_allsims, sexcats_allsims = local_catalogs, local_meds_inputs, local_truthcats, local_sexcats
            for i in range(1,size):
                rec_data = comm.recv(source=i)
                if rec_data is not None:
                    rec_catalogs, rec_meds_inputs, rec_truthcats, rec_sexcats = rec_data
                    catalogs_allsims+=rec_catalogs
                    meds_inputs_allsims+=rec_meds_inputs
                    truthcats_allsims+=rec_truthcats
                    sexcats_allsims+=rec_sexcats
            assert len(catalogs_allsims)==known.n_sims

            #Save stacked truth/measured catalogs from all sims for convenience
            truthcat_all=astropy.table.vstack(truthcats_allsims)
            sexcat_all=astropy.table.vstack(sexcats_allsims)
            truthcat_all.write(os.path.join(known.parent_output,"truthcat_full.fits"),format='fits',overwrite=True)
            sexcat_all.write(os.path.join(known.parent_output,"measuredcat_full.fits"),format='fits',overwrite=True)
            CoaddSetup.runlogger.info("Concatenated and saved truth/measured catalogs from all sims")
        else:
            send_data=(local_catalogs, local_meds_inputs, local_truthcats, local_sexcats)
            comm.send(send_data)

    #So...now we've done the coadd bit, we need to allocate the SE images across the processes
    #Send sim catalogs to all processes
    catalogs_allsims = comm.bcast(catalogs_allsims, root=0)
    #Get 'not_drawn' array for each sim (need to copy since this will be overwritten)
    coadd_not_drawn_allsims = [np.copy(c.galaxy['not_drawn']) for c in catalogs_allsims]
    #Also send gspcatalog (galsim parameters)
    gspcatalog = comm.bcast(gspcatalog, root=0)
    CoaddSetup['gspcatalog']=gspcatalog

    #Now split se images across available processes - the idea is that each process gets a 
    #certain number of single-epoch images, and then for each image, loops through all the 
    #simulations, thus minimising the number of times the images are read.
    se_image_inds=np.arange(0,len(coadd.srclist)) #Use this to index source list
    local_se_image_inds=np.array_split(se_image_inds,size)[rank]
    print 'rank %d doing images:'%rank, local_se_image_inds
    if len(local_se_image_inds)==0:
        comm.send(None)
        return
    coadd.srclist.sort(key=lambda x: x['red_cat']) #Sort to make sure each process has source list in same order
    local_srclist=[coadd.srclist[i] for i in local_se_image_inds]
    local_file_ids=1+local_se_image_inds #'file_id' used in meds creation to find correct image/psf for each stamp. Add 1 because we add coadd to beginning of source list later.

    #set up output - just collect single-epoch objects in a list
    local_se_objects=[]
    CoaddSetup.runlogger.info('rank %d doing %d SE images'%(rank,len(local_se_image_inds)))
    for i_se,se_info in enumerate(local_srclist):
        CoaddSetup.runlogger.info('processing image %s'%se_info['red_image'])
        file_id=local_file_ids[i_se]
        se_info['file_id']=file_id
        BalrogSetup=SEBalrogSetup(CoaddSetup,se_info)
        # Get the subsampled flux and weightmap images, along with the PSF model and WCS.
        BalrogSetup.bigImage, BalrogSetup.weight, BalrogSetup.badpix, BalrogSetup.sky, BalrogSetup.psfmodel, BalrogSetup.wcs = BalrogSetup.ReadImages()
        if known.no_noise:
            BalrogSetup.bigImage*=0.

        #Loop through sim catalogs
        for i_sim,catalog in enumerate(catalogs_allsims):
            BalrogSetup.runlogger.info('i_sim = %d'%i_sim)
            sim_se_objects = SERun(catalog, se_info, BalrogSetup)
            BalrogSetup.runlogger.info('%d objects in image %s'%(len(sim_se_objects),se_info['red_image']))
            for s in sim_se_objects:
                local_se_objects.append((i_sim, s[0], s[1]))
        Cleanup(BalrogSetup)

    #Again, rank 0 collects outputs - in this case single-epoch objects to be added to multi-epoch objects, which
    #are then added to meds file.
    if rank==0:
        BalrogSetup.runlogger.info('rank 0 collecting se_objects')
        se_objects=local_se_objects
        BalrogSetup.runlogger.info('rank, # of objects')
        BalrogSetup.runlogger.info('0, %d'%len(se_objects))
        for i in range(1,size):
            se_objects += comm.recv(source=i)
            BalrogSetup.runlogger.info('%d, %d'%(i,len(se_objects)))
        for se_obj in se_objects:
            i_sim,i_obj=se_obj[0],se_obj[1]
            meds_inputs_allsims[i_sim][i_obj].add_SEObject(se_obj[2])

        srclist=[{'id':coadd['image_id'], 'red_image':coadd['image_url']}]+coadd.srclist

        for i_sim,meds_inputs in enumerate(meds_inputs_allsims):
            MEobjs=[]
            not_drawn=coadd_not_drawn_allsims[i_sim]
            for i,meds_input in enumerate(meds_inputs):
                if len(meds_input['images'])<1: 
                    BalrogSetup.runlogger.info("no exposures for galaxy",i)
                    continue
                if not_drawn[i]:
                    print 'object %d not drawn in coadd, skipping'%i
                    continue
                MEobjs.append(meds_input.make_MEobj(dummy_segs=False))
            #append coadd info to srclist
            meds_name=os.path.join(known.parent_output,"%s-%s-meds.%d.fits")%(known.tilename,known.band,i_sim)
            des_meds.write_meds(meds_name, MEobjs, srclist=srclist, clobber=True)
            BalrogSetup.runlogger.info("wrote %d objects to meds file %s"%(len(meds_inputs),meds_name))
    else:
        comm.send(local_se_objects)
            
def TileArgs(parser):
    parser.add_argument("tilename",type=str)
    parser.add_argument("parent_output",type=str,help="the output directory (each simulation will be saved to a different directory within this one")
    parser.add_argument("--no_mpi", default=False, action='store_true', help="run on single core")
    parser.add_argument("--band",type=str,default='i')
    parser.add_argument("--coadd_nproc",type=int,default=None,
                        help='limit the coadd part (more memory intensive??) to this many processes')
    parser.add_argument("--no_meds",action='store_true',default=False)
    parser.add_argument("--no_noise",action='store_true',default=False,help="make 'noiseless' images i.e. set background image to zero everywhere")
    parser.add_argument("--n_sims",type=int,default=1, help="number of simulations (of same tile)")

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
        from mpi4py import MPI
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
