import numpy
import galsim
import os

# these image stamp sizes are available in MEDS format
BOX_SIZES = [32,48,64,96,128,192,256]
# while creating the meds file, all the data is stored in memory, and then written to disc once all
# the necessary images have been created.
# You can control the amound of memory used to prevent jamming your system.
MAX_MEMORY = 1e9
# Maximum number of exposures allowed per galaxy (incl. coadd)
MAX_NCUTOUTS = 10
# flags for unavailable data
EMPTY_START_INDEX = 9999
EMPTY_JAC_diag    = 1
EMPTY_JAC_offdiag = 0
EMPTY_SHIFT = 0

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
        self['file_ids']=[]

    def add(self,image,weight,badpix,seg,wcs,orig_row,orig_col,orig_start_row,orig_start_col,cutout_row,cutout_col,file_id):
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
        self['file_ids'].append(file_id)

    def add_SEObject(self, se):
        self['images'].append(se.image)
        self['weights'].append(se.weight)
        self['badpixs'].append(se.badpix)
        self['segs'].append(se.seg)
        self['wcss'].append(se.wcs)
        self['orig_rows'].append(se.orig_row)
        self['orig_cols'].append(se.orig_col)
        self['orig_start_rows'].append(se.orig_start_row)
        self['orig_start_cols'].append(se.orig_start_col)
        self['cutout_rows'].append(se.cutout_row)
        self['cutout_cols'].append(se.cutout_col)
        self['file_ids'].append(se.file_id)

    def make_MEobj(self,dummy_segs=False):
        if dummy_segs:
            self['segs']=None
        return MultiExposureObject(images=self['images'], weights=self['weights'], badpix=self['badpixs'], segs=self['segs'], 
                                            wcs=self['wcss'], id=self['id'], orig_rows=self['orig_rows'],
                                            orig_cols=self['orig_cols'], orig_start_rows=self['orig_start_rows'], 
                                            orig_start_cols=self['orig_start_cols'], cutout_rows=self['cutout_rows'], 
                                            cutout_cols=self['cutout_cols'],file_ids=self['file_ids'])

class SEObject(object):
    def __init__(self, image,weight,badpix,seg,wcs,orig_row,orig_col,orig_start_row,
                 orig_start_col,cutout_row,cutout_col,file_id):
        self.image=image
        self.weight=weight
        self.badpix=badpix
        self.seg=seg
        self.wcs=wcs
        self.orig_row=orig_row
        self.orig_col=orig_col
        self.orig_start_row=orig_start_row
        self.orig_start_col=orig_start_col
        self.cutout_row=cutout_row
        self.cutout_col=cutout_col
        self.file_id=file_id

class MultiExposureObject(object):
    """
    A class containing exposures for single object, along with other information.

    Initialization
    --------------

    @param images       List of images of the object (GalSim Images).
    @param weights      List of weight maps (GalSim Images). [default: None]
    @param badpix       List of bad pixel masks (GalSim Images). [default: None]
    @param segs         List of segmentation maps (GalSim Images). [default: None]
    @param wcs          List of WCS transformations (GalSim AffineTransforms). [default: None]    
    @param ids           Galaxy id. [default: 0]
    @param number       Sextractor NUMBER [default: 0]
    #@param file_ids    #Not sure what this is so not doing it 
    @param orig_rows    List of orig_row (see MEDS documentation for defn of orig_row and following params)
    @param orig_cols
    @param orig_start_rows
    @param orig_start_cols
    @param cutout_rows
    @param cutout_cols
    @param image_ids

    Attributes
    ----------

    self.images         List of images of the object (GalSim Images).
    self.weights        List of weight maps (GalSim Images).
    self.segs           List of segmentation masks (GalSim Images).
    self.wcs            List of WCS transformations (GalSim AffineTransforms).
    self.n_cutouts      Number of exposures.
    self.box_size       Size of each exposure image.
    self.ids
    self.number
    self.orig_rows
    self.orig_cols
    self.orig_start_rows
    self.orig_start_cols
    self.cutout_rows
    self.cutout_cols
    self.image_ids

    Module level variables
    ----------------------

    Images, weights and segs have to be square numpy arrays with size in
    BOX_SIZES = [32,48,64,96,128,196,256].
    Number of exposures for all lists (images,weights,segs,wcs) have to be the same and smaller 
    than MAX_NCUTOUTS (default 11).
    """

    def __init__(self, images, weights=None, badpix=None, segs=None, wcs=None, id=0, orig_rows=None,
                 orig_cols=None, orig_start_rows=None, orig_start_cols=None, cutout_rows=None, cutout_cols=None, 
                 file_ids=None, number=None):
        
        # assign the ID
        self.id = id
        if number==None:
            self.number = id+1
        else:
            self.number = number

        # check if images is a list
        if not isinstance(images,list):
            raise TypeError('images should be a list')

        # get number of cutouts from image list
        self.images = images
        # get box size from the first image
        self.box_size = self.images[0].array.shape[0]
        self.n_cutouts = len(self.images)

        for entry in ["orig_rows","orig_cols","orig_start_rows","orig_start_cols","cutout_rows","cutout_cols","file_ids"]:
            if eval(entry) is None:
                setattr(self,entry,[-9999]*self.n_cutouts)
            else:
                assert len(eval(entry))==self.n_cutouts
                setattr(self,entry,eval(entry))

        # see if there are cutouts
        if self.n_cutouts < 1:
            raise ValueError('no cutouts in this object')

        # check if the box size is correct
        if self.box_size not in BOX_SIZES:
            # raise ValueError('box size should be in  [32,48,64,96,128,196,256], is %d' % box_size)
            raise ValueError( 'box size should be in '+str(BOX_SIZES)+', is '+str(self.box_size) )

        # check if weights, segs and wcs were supplied. If not, create sensible values.
        if weights != None:
            self.weights = weights
        else:
            self.weights = [galsim.Image(self.box_size, self.box_size, init_value=1)]*self.n_cutouts

        # check segmaps
        if segs != None:
            self.segs = segs
        else:
            self.segs = [galsim.ImageI(self.box_size, self.box_size, init_value=self.number)]*self.n_cutouts

        # check wcs
        if wcs != None:
            self.wcs = wcs
        else:
            # Get the wcs from the images.  Probably just the pixel scale.
            self.wcs = [ im.wcs.jacobian().setOrigin(im.trueCenter()) for im in self.images ]

         # check if weights,segs,jacks are lists
        if not isinstance(self.weights,list):
            raise TypeError('weights should be a list')
        if not isinstance(self.segs,list):
            raise TypeError('segs should be a list')
        if not isinstance(self.wcs,list):
            raise TypeError('wcs should be a list')


        # loop through the images and check if they are of the same size
        for extname in ('images','weights','segs'):

            # get the class field
            ext = eval('self.' + extname )

            # loop through exposures
            for icutout,cutout in enumerate(ext):

                # get the sizes of array
                nx=cutout.array.shape[0]
                ny=cutout.array.shape[1]

                # x and y size should be the same
                if nx != ny:
                    raise ValueError('%s should be square and is %d x %d' % (extname,nx,ny))

                # check if box size is correct

                    raise ValueError('%s object %d has size %d and should be %d' % 
                            ( extname,icutout,nx,self.box_size ) )

        # see if the number of Jacobians is right
        if len(self.wcs) != self.n_cutouts:
            raise ValueError('number of Jacobians is %d is not equal to number of cutouts %d'%
                    ( len(self.wcs),self.n_cutouts ) )

        # check each Jacobian
        for jac in self.wcs:
            # should ba an AffineTransform instance
            if not (isinstance(jac, galsim.AffineTransform) or isinstance(jac, galsim.wcs.LocalWCS)):
                raise TypeError('wcs list should contain AffineTransform or LocalWCS objects')
            

def write_meds(file_name, obj_list, srclist=None, clobber=True):
    """
    @brief Writes the galaxy, weights, segmaps images to a MEDS file.

    Arguments:
    ----------
    @param file_name:    Name of meds file to be written
    @param obj_list:     List of MultiExposureObjects
    @param clobber       Setting `clobber=True` when `file_name` is given will silently overwrite 
                         existing files. (Default `clobber = True`.)
    """

    import numpy
    import sys
    from galsim import pyfits

    # initialise the catalog
    cat = {}
    cat['ncutout'] = []
    cat['box_size'] = []
    cat['start_row'] = []
    cat['orig_row']=[]
    cat['orig_col']=[]
    cat['orig_start_row']=[]
    cat['orig_start_col']=[]
    cat['cutout_row']=[]
    cat['cutout_col']=[]
    cat['id'] = []
    cat['number'] = []
    cat['dudrow'] = []
    cat['dudcol'] = []
    cat['dvdrow'] = []
    cat['dvdcol'] = []
    cat['file_id'] = []

    # initialise the image vectors
    vec = {}
    vec['image'] = []
    vec['seg'] = []
    vec['weight'] = []

    # initialise the image vector index
    n_vec = 0
    
    # get number of objects
    n_obj = len(obj_list)

    # loop over objects
    for obj in obj_list:

        # initialise the start indices for each image
        start_row = numpy.ones(MAX_NCUTOUTS)*EMPTY_START_INDEX
        dudrow = numpy.ones(MAX_NCUTOUTS)*EMPTY_JAC_diag 
        dudcol = numpy.ones(MAX_NCUTOUTS)*EMPTY_JAC_offdiag
        dvdrow = numpy.ones(MAX_NCUTOUTS)*EMPTY_JAC_offdiag
        dvdcol = numpy.ones(MAX_NCUTOUTS)*EMPTY_JAC_diag
        orig_row = numpy.ones(MAX_NCUTOUTS)*EMPTY_SHIFT
        orig_col = numpy.ones(MAX_NCUTOUTS)*EMPTY_SHIFT
        orig_start_row = numpy.ones(MAX_NCUTOUTS)*EMPTY_SHIFT
        orig_start_col = numpy.ones(MAX_NCUTOUTS)*EMPTY_SHIFT
        cutout_row = numpy.ones(MAX_NCUTOUTS)*EMPTY_SHIFT
        cutout_col = numpy.ones(MAX_NCUTOUTS)*EMPTY_SHIFT
        file_id = numpy.ones(MAX_NCUTOUTS)*-9999

        # get the number of cutouts (exposures)
        n_cutout = obj.n_cutouts
        
        # append the catalog for this object
        cat['ncutout'].append(n_cutout)
        cat['box_size'].append(obj.box_size)
        cat['id'].append(obj.id)
        cat['number'].append(obj.number)

        # loop over cutouts
        for i in range(n_cutout):
            
            # assign the start row to the end of image vector
            start_row[i] = n_vec
            # update n_vec to point to the end of image vector
            n_vec += len(obj.images[i].array.flatten()) 


            # append the image vectors
            vec['image'].append(obj.images[i].array.flatten())
            vec['seg'].append(obj.segs[i].array.flatten())
            vec['weight'].append(obj.weights[i].array.flatten())


            # append the Jacobian
            dudrow[i] = obj.wcs[i].dudy
            dudcol[i] = obj.wcs[i].dudx
            dvdrow[i] = obj.wcs[i].dvdy
            dvdcol[i] = obj.wcs[i].dvdx
            orig_row[i] = obj.orig_rows[i]
            orig_col[i]= obj.orig_cols[i]
            orig_start_row[i]= obj.orig_start_rows[i]
            orig_start_col[i]= obj.orig_start_cols[i]
            cutout_row[i]= obj.cutout_rows[i]
            cutout_col[i]= obj.cutout_cols[i]
            file_id[i]= obj.file_ids[i]

            # check if we are running out of memory
            if sys.getsizeof(vec) > MAX_MEMORY:
                raise MemoryError(
                    'Running out of memory > %1.0fGB '%MAX_MEMORY/1.e9 +
                    '- you can increase the limit by changing MAX_MEMORY')

        # update the start rows fields in the catalog
        cat['start_row'].append(start_row)

        # add lists of Jacobians
        cat['dudrow'].append(dudrow)
        cat['dudcol'].append(dudcol)
        cat['dvdrow'].append(dvdrow)
        cat['dvdcol'].append(dvdcol)
        cat['orig_row'].append(orig_row)
        cat['orig_col'].append(orig_col)
        cat['orig_start_row'].append(orig_start_row)
        cat['orig_start_col'].append(orig_start_col)
        cat['cutout_row'].append(cutout_row)
        cat['cutout_col'].append(cutout_col)
        cat['file_id'].append(file_id)

    # concatenate list to one big vector
    vec['image'] = numpy.concatenate(vec['image'])
    vec['seg'] = numpy.concatenate(vec['seg'])
    vec['weight'] = numpy.concatenate(vec['weight'])

    # get the primary HDU
    primary = pyfits.PrimaryHDU()

    # second hdu is the object_data
    cols = []
    cols.append( pyfits.Column(name='ncutout', format='I', array=cat['ncutout'] ) )
    cols.append( pyfits.Column(name='id', format='K', array=cat['id'] ) )
    cols.append( pyfits.Column(name='number', format='K', array=cat['number'] ) )
    cols.append( pyfits.Column(name='box_size', format='J', array=cat['box_size'] ) )
    cols.append( pyfits.Column(name='file_id', format='%dI' % MAX_NCUTOUTS, array=cat['file_id']) )
    cols.append( pyfits.Column(name='start_row', format='%dJ' % MAX_NCUTOUTS,
                               array=numpy.array(cat['start_row'])) )
    cols.append( pyfits.Column(name='orig_row', format='%dD' % MAX_NCUTOUTS, array=numpy.array(cat['orig_row'])) )
    cols.append( pyfits.Column(name='orig_col', format='%dD' % MAX_NCUTOUTS, array=numpy.array(cat['orig_col'])) )
    cols.append( pyfits.Column(name='orig_start_row', format='%dI' % MAX_NCUTOUTS, array=numpy.array(cat['orig_start_row'])) )
    cols.append( pyfits.Column(name='orig_start_col', format='%dI' % MAX_NCUTOUTS, array=numpy.array(cat['orig_start_col'])) )
    cols.append( pyfits.Column(name='dudrow', format='%dD'% MAX_NCUTOUTS,
                               array=numpy.array(cat['dudrow']) ) )
    cols.append( pyfits.Column(name='dudcol', format='%dD'% MAX_NCUTOUTS,
                               array=numpy.array(cat['dudcol']) ) )
    cols.append( pyfits.Column(name='dvdrow', format='%dD'% MAX_NCUTOUTS,
                               array=numpy.array(cat['dvdrow']) ) )
    cols.append( pyfits.Column(name='dvdcol', format='%dD'% MAX_NCUTOUTS,
                               array=numpy.array(cat['dvdcol']) ) )
    cols.append( pyfits.Column(name='cutout_row', format='%dD'% MAX_NCUTOUTS,
                               array=numpy.array(cat['cutout_row']) ) )
    cols.append( pyfits.Column(name='cutout_col', format='%dD'% MAX_NCUTOUTS,
                               array=numpy.array(cat['cutout_col']) ) )


    object_data = pyfits.new_table(pyfits.ColDefs(cols))
    object_data.update_ext_name('object_data')

    # third hdu is image_info
    cols = []
    if srclist is not None:
        image_ids=[entry['id'] for entry in srclist]
        image_paths=[entry['red_image'] for entry in srclist]
        cols.append( pyfits.Column(name='image_id', format='J', array=numpy.array(image_ids)) )
        cols.append( pyfits.Column(name='image_path', format='A256',   array=numpy.array(image_paths) ) )
    else:
        cols.append( pyfits.Column(name='image_id', format='I', array=numpy.array([0])))
        cols.append( pyfits.Column(name='image_path', format='A256',   array= ['not_set']) )
    cols.append( pyfits.Column(name='sky_path',   format='A256',   array=['generated_by_galsim'] ) )
    cols.append( pyfits.Column(name='seg_path',   format='A256',   array=['generated_by_galsim'] ) )
    image_info = pyfits.new_table(pyfits.ColDefs(cols))
    image_info.update_ext_name('image_info')

    # fourth hdu is metadata
    cols = []
    cols.append( pyfits.Column(name='DESDATA',       format='A256', array=[os.environ['DESDATA']] ))
    cols.append( pyfits.Column(name='cat_file',      format='A256', array=['generated_by_galsim'] ))
    cols.append( pyfits.Column(name='coadd_file',    format='A256', array=['generated_by_galsim'] ))
    cols.append( pyfits.Column(name='coadd_hdu',     format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='coadd_seg_hdu', format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='coadd_srclist', format='A256', array=['generated_by_galsim'] ))
    cols.append( pyfits.Column(name='coadd_wt_hdu',  format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='coaddcat_file', format='A256', array=['generated_by_galsim'] ))
    cols.append( pyfits.Column(name='coaddseg_file', format='A256', array=['generated_by_galsim'] ))
    cols.append( pyfits.Column(name='cutout_file',   format='A256', array=['generated_by_galsim'] ))
    cols.append( pyfits.Column(name='max_boxsize',   format='A3',   array=['x']                   ))
    cols.append( pyfits.Column(name='medsconf',      format='A3',   array=['x']                   ))
    cols.append( pyfits.Column(name='min_boxsize',   format='A2',   array=['x']                   ))
    cols.append( pyfits.Column(name='se_badpix_hdu', format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='se_hdu',        format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='se_wt_hdu',     format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='seg_hdu',       format='A1',   array=['x']                   ))
    cols.append( pyfits.Column(name='sky_hdu',       format='A1',   array=['x']                   ))
    metadata = pyfits.new_table(pyfits.ColDefs(cols))
    metadata.update_ext_name('metadata')

    # rest of HDUs are image vectors
    image_cutouts   = pyfits.ImageHDU( vec['image'] , name='image_cutouts'  )
    weight_cutouts  = pyfits.ImageHDU( vec['weight'], name='weight_cutouts' )
    seg_cutouts     = pyfits.ImageHDU( vec['seg']   , name='seg_cutouts'    )

    # write all
    hdu_list = pyfits.HDUList([
        primary,
        object_data,
        image_info,
        metadata,
        image_cutouts, 
        weight_cutouts,
        seg_cutouts
    ])
    hdu_list.writeto(file_name,clobber=clobber)


# Now add this to the config framework.
import galsim.config

# Make this a valid output type:
galsim.config.process.valid_output_types['des_meds'] = (

    'galsim.des.BuildMEDS',      # Function that builds the objects using config
    'galsim.des.GetNObjForMEDS', # Function that calculates the number of objects
    True,   # Takes nproc argument
    False,  # Takes *_file_name arguments for psf, weight, badpix
    False)  # Takes *_hdu arguments for psf, weight, badpix

def BuildMEDS(file_name, config, nproc=1, logger=None, file_num=0, image_num=0, obj_num=0):
    """
    Build a meds file as specified in config.

    @param file_name         The name of the output file.
    @param config            A configuration dict.
    @param nproc             How many processes to use. [default: 1]
    @param logger            If given, a logger object to log progress. [default: None]
    @param file_num          If given, the current file_num. [default: 0]
    @param image_num         If given, the current image_num. [default: 0]
    @param obj_num           If given, the current obj_num. [default: 0]

    @returns the time taken to build file
    """
    import time
    t1 = time.time()

    config['seq_index'] = file_num
    config['file_num'] = file_num

    ignore = [ 'file_name', 'dir', 'nfiles', 'psf', 'weight', 'badpix', 'nproc' ]
    req = { 'nobjects' : int , 'nstamps_per_object' : int }
    params = galsim.config.GetAllParams(config['output'],'output',config,ignore=ignore,req=req)[0]

    nobjects = params['nobjects']
    nstamps_per_object = params['nstamps_per_object']
    ntot = nobjects * nstamps_per_object

    all_images = galsim.config.BuildImages(
        ntot, config=config, nproc=nproc, logger=logger, obj_num=obj_num,
        make_psf_image=False, make_weight_image=True, make_badpix_image=True)

    main_images = all_images[0]
    weight_images = all_images[2]
    badpix_images = all_images[3]

    obj_list = []
    for i in range(nobjects):
        k1 = i*nstamps_per_object
        k2 = (i+1)*nstamps_per_object
        obj = MultiExposureObject(images = main_images[k1:k2], 
                                  weights = weight_images[k1:k2],
                                  badpix = badpix_images[k1:k2])
        obj_list.append(obj)

    write_meds(file_name, obj_list)

    t2 = time.time()
    return t2-t1

def GetNObjForMEDS(config, file_num, image_num):
    ignore = [ 'file_name', 'dir', 'nfiles', 'psf', 'weight', 'badpix', 'nproc' ]
    req = { 'nobjects' : int , 'nstamps_per_object' : int }
    params = galsim.config.GetAllParams(config['output'],'output',config,ignore=ignore,req=req)[0]
    config['seq_index'] = file_num

    if 'image' in config and 'type' in config['image']:
        image_type = config['image']['type']
        if image_type != 'Single':
            raise AttibuteError("MEDS files are not compatible with image type %s."%image_type)

    nobjects = params['nobjects']
    nstamps_per_object = params['nstamps_per_object']
    ntot = nobjects * nstamps_per_object

    # nobj is a list of nobj per image.
    # The MEDS file is considered to only have a single image, so the list has only 1 element.
    nobj = [ ntot ]
    return nobj
