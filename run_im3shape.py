#Single epoch driver for im3shape (a model-fitting shape measurement code) for use on e.g. a coadd
#You need to install im3shape (https://bitbucket.org/joezuntz/im3shape) to use it
#The SExtractor segmentation map can (and usually should!) be used to mask neighbouring objects when fitting
#This process is a bit more complicated in Balrog, since within Balrog SExtractor will be run in assoc mode,
#which won't produce the necessary segmentation mask (you need to run SExtractor in normal mode, after inserting
#galaxies). Having a segmentation map both from the full run, and the assoc run allows both masking, and the matching
#of objects identified in the full run, to the assoc run. This allows the use of uberseg. Hence the option '--seg_filnames' 
#can be used with two segmentation maps.

import py3shape
import numpy as np
import galsim
import galsim.des
import os
import time
import pyfits
import argparse
import logging
from scipy.stats import mode
from py3shape.analyse import get_psf_params

#im3shape defaults
class DEFAULT:
    BALROG_DIR = '/'.join( os.path.realpath(__file__).split('/')[:-1] )
    ASTRO_DIR = os.path.join( BALROG_DIR, 'astro_config' )
    IM3_OPTIONS = os.path.join( ASTRO_DIR, 'im3_options.ini' )

#create galaxy catalog dictionary from fits catalog
def get_gal_cat(args):
    gal_cat={}
    ext=2
    if args.cat_filename[-1]==']' and args.cat_filename[-3]=='[':
        try:
            ext=int(args.cat_filename[-2])
            args.cat_filename = args.cat_filename[:-3]
        except TypeError:
            print 'specify catalog as fits filename (optionally) with extension appended in square brackets'
            print 'e.g. bullet_cat.fits[2], 2 is the default extension'
    gal_data = pyfits.getdata(args.cat_filename,ext)
    gal_cat['X_IMAGE'] = gal_data[args.x_image_col]
    gal_cat['Y_IMAGE'] = gal_data[args.y_image_col]
    gal_cat['ID'] = gal_data[args.id_col]
    gal_cat['RA'] = gal_data[args.ra_col]
    gal_cat['DEC'] = gal_data[args.dec_col]
    return gal_cat

#Read Balrog output catalog to get positions of galaxies
def read_Balrog_output(catalog):
	return

def getPSFExarray(psfex, pos, nsidex, nsidey, upsampling=1, offset=None):
    """Return an image of the PSFEx model of the PSF as a NumPy array.

    Arguments
    ---------
    psfex       A galsim.des.PSFEx instance opened using, for example,
                `psfex = galsim.des.DES_PSFEx(psfex_file_name)`.
    pos         galsim PositionD

    nsidex      Size of PSF image along x [pixels]
    nsidey      Size of PSF image along y [pixels]
    upsampling  Upsampling (see Zuntz et al 2013)

    Returns a NumPy array with shape (nsidey, nsidex) - note the reversal of y and x to match the
    NumPy internal [y, x] style array ordering.  This is to ensure that `pyfits.writeto()` using the
    ouput array creates FITS-compliant output.
    """
    image = galsim.ImageD(nsidex, nsidey)    
    psf = psfex.getPSF(pos)
    psf.draw(image, scale=1./upsampling, offset=offset)
    return image.array

def get_stamp(image,x_pos,y_pos,stamp_size):
    "Get postage stamp to run im3shape on"
    "Expects galsim image, SExtractor X_IMAGE,Y_IMAGE and stamp size as input"
    "Returns stamp, x0 and y0 - coordinates of object in stamp (starting from 0)"
    half=stamp_size/2.
    x_min,x_max,y_min,y_max = int(x_pos-half),int(x_pos+half),int(y_pos-half),int(y_pos+half)
    over_edge = (x_min<0) | (y_min<0) | (x_max>image.array.shape[1]) | (y_max>image.array.shape[0])
    if over_edge:
        logging.warning('galaxy stamp overlaps edge of image, skipping')
        return 1
    subBounds = galsim.BoundsI(x_min,x_max-1,y_min,y_max-1)
    stamp = image[subBounds]
    #stamp = img_data[ypos-half:ypos+half, xpos-half:xpos+half]
    x0,y0 = x_pos - x_min, y_pos - y_min
    return stamp.array, x0, y0



def convert_g_image2sky(local, g1image, g2image):
    """Return the ellipticity (g1sky, g2sky) in sky coordinates corresponding to the input
    (g1image, g2image).

    Uses the ellipticity convention |g| = (a-b)/(a+b).
    
    Currently only works for scalar input g1image, g2image, but can be called in list
    comprehensions.
    """
    Aimage2sky = np.array([[local.dudx,local.dudy],[local.dvdx,local.dvdy]])
    e1image, e2image = py3shape.utils.convert_e_linear_to_quadratic(g1image, g2image)
    # Build the ellipticity matrix
    Ei = np.array(((1. + e1image, e2image), (e2image, 1. - e1image)))
    # Perform the transformation
    Es = np.dot(np.dot(Aimage2sky, Ei), Aimage2sky.T)
    # Extract the ellipticities and convert back to the linear |g|=(a-b)/(a+b) ellips
    Estrace = Es.trace()
    e1sky, e2sky = (Es[0, 0] - Es[1, 1]) / Estrace, 2. * Es[0, 1] / Estrace
    g1sky, g2sky = py3shape.utils.convert_e_quadratic_to_linear(e1sky, e2sky)
    return g1sky, g2sky

def seg_to_mask_basic(identifier,seg_stamp):
    if identifier not in seg_stamp:
        print 'target object not in segmentation mask...'
    mask=np.ones(seg_stamp.shape)
    mask[((seg_stamp!=0) & (seg_stamp!=identifier))] = 0
    return mask

def seg_to_mask_uber(identifier,seg_stamp):
    #Object id in seg map should be 
    #First check that expected object is in seg map
    if identifier not in seg_stamp:
        print 'ID not in seg...'
        raise ValueError
    #First get all indices of all seg map pixels which contain an object i.e. are not equal to zero
    obj_inds = np.where(seg_stamp!=0)
    mask=np.ones(seg_stamp.shape)
    #Then loop through pixels in seg map, check which obj ind it is closest to.
    #If the closest obj ind does not correspond to the target, set this pixel in the weight map to zero.
    for i,row in enumerate(seg_stamp):
        for j, element in enumerate(row):
            obj_dists = (i-obj_inds[0])**2 + (j-obj_inds[1])**2
            ind_min=np.argmin(obj_dists)
            if seg_stamp[obj_inds[0][ind_min],obj_inds[1][ind_min]] != identifier:
                mask[i,j] = 0.
    return mask

def get_fits_extension(input_filename):
    ext=0
    output_filename=input_filename
    if input_filename[-1]==']' and input_filename[-3]=='[':
        try:
            ext=int(input_filename[-2])
            output_filename = input_filename[:-3]
        except TypeError:
            print 'specify image/weight/seg as fits filename with extension appended in square brackets'
            print 'e.g. bullet.fits[2], default extension is 0'
            raise
    return output_filename,ext


def main(args):

    # load the options file
    options=py3shape.Options(args.ini_filename)
    options.validate()

    #Get extra command line ini file options
    if args.extra_options is not None:
        for opt in args.extra_options:
            key,val = opt.split('=')
            setattr(options, key, val)

    # Read in the FITS data
    img_filename,img_ext = get_fits_extension(args.img_filename)
    if args.weight_filename:
        weight_filename,weight_ext = get_fits_extension(args.weight_filename)
        weight_gs = galsim.fits.read(weight_filename,hdu=weight_ext)

    if args.seg_filenames:
        seg_imgs=[]
        for seg_file in args.seg_filenames:
            seg_filename, seg_ext = get_fits_extension(seg_file)
            seg_imgs.append(galsim.fits.read(seg_filename,hdu=seg_ext))

    #Read in image
    image_gs=galsim.fits.read(img_filename,hdu=img_ext)
    img_data=image_gs.array
    #Get wcs info:
    wcs = image_gs.wcs
    #Read psfex file
    try:
        psfexer = galsim.des.DES_PSFEx(args.psf_filename,args.img_filename)
    except Exception:
        logging.error('failed to read psf file %s',(args.psf_filename))
        raise


    #Read in catalog data
    gal_cat = get_gal_cat(args)
        
    # Create i3_image of certain stamp size
    stamp_size = options.stamp_size

    #overwrite the output filename
    options.output_filename = args.out_filename
    options.save_output = False
    extra_cols=['e1_sky','e2_sky']
    if args.psf_props:
        extra_cols+=['psf_fwhm','psf_e1','psf_e2']
    if args.masking_type=='all':
        extra_cols+=['e1_nomask','e2_nomask','e1_uber','e2_uber']

    #exclude following columns...this is messy, probably better to define new single epoch output object in im3shape...
    excluded_cols = ['exposure_x','exposure_y','exposure_e1','exposure_e2','exposure_chi2',
                     'mean_flux','exposure_residual_stdev']
    output = py3shape.output.Output(args.out_filename, options, excluded_cols=excluded_cols)

    extra_lines=['driver: Balrog/run_im3shape.py', 
                'ini file: %s' % (args.ini_filename,),
                'catalog file: %s' % (args.cat_filename,),
                'image file: %s' % (args.img_filename,),
                'psfex file: %s' % (args.psf_filename,),
                'first object: %s' % (args.first,),
                'last object: %s' % (args.last,)]

    output.write_header(include_radec=True, extra_cols=extra_cols,extra_lines=extra_lines)

    #Write extra lines to info() log
    for extra_line in extra_lines:
        logging.info(extra_line)

    #main galaxy loop
    ID=gal_cat['ID']
    X_IMAGE=gal_cat['X_IMAGE']
    Y_IMAGE=gal_cat['Y_IMAGE']
    RA,DEC=gal_cat['RA'],gal_cat['DEC']
    if args.last==None or args.last>len(ID)-1:
        args.last = len(ID)-1

    # Loop over objects in catalog
    logging.info('Analyzing %d/%d galaxies' % (min(args.last-args.first+1,len(gal_cat['ID'])), len(gal_cat['ID'])))
    start_time = time.clock()

    for i in range(args.first,args.last+1):
        # Read galaxy catalog entry
        identifier = ID[i]
        print identifier
        xpos = X_IMAGE[i]
        ypos = Y_IMAGE[i]
        #Get image stamp
        stamp_stuff = get_stamp(image_gs,X_IMAGE[i],Y_IMAGE[i],options.stamp_size)
        try:
            stamp,x0,y0 = stamp_stuff
        except TypeError:
            #if stamp_stuff==1:
            #    logging.warning('galaxy %d stamp overlaps edge of image, skipping',identifier)
            #else:
            #    logging.warning('failed to get stamp for galaxy %d, skipping',identifier)
            continue

        #Get weight stamp and mask stamp:
        if args.weight_filename:
            weight_stamp,_,_ = get_stamp(weight_gs,xpos,ypos,options.stamp_size)
            #Set negative values to zero
            weight_stamp[weight_stamp<0]=0
        else:
            weight_stamp=None
        if not args.masking_type=='none':
            #Get seg stamp data...decide what to do depending on how many seg files provided.
            #If only one provided, this should be the seg map from the original image. Set all non-zero pixels in this image
            #to -1, then create mask, that way they won't match ID of any detected simulated objects.
            if len(seg_imgs)==1:
                seg_stamp,_,_ = get_stamp(seg_imgs[0],xpos,ypos,options.stamp_size)
                #np.set_printoptions(threshold=np.nan)
                #print 'seg_stamp',seg_stamp
                seg_stamp[(seg_stamp!=0)] = -1
                mask=seg_to_mask_basic(identifier,seg_stamp)
                mask_stamp = py3shape.Image(mask)
            #If two provided, combine them in the following way:
            if len(seg_imgs)==2:
                #np.set_printoptions(threshold=np.nan)
                noassoc_seg_stamp,_,_ = get_stamp(seg_imgs[0],xpos,ypos,options.stamp_size)
                assoc_seg_stamp,_,_ = get_stamp(seg_imgs[1],xpos,ypos,options.stamp_size)
                #Find target object pixels in assoc stamp, and find which pixel value they overlap
                #most with in noassoc stamp. What if this is zero? Can't see why this would happen if same SExtractor settings
                #used for both seg maps, so for no just log a warning and skip object if this happens...
                assoc_obj_inds=np.where(assoc_seg_stamp==identifier)
                noassoc_seg_vals=noassoc_seg_stamp[assoc_obj_inds]
                try:
                    noassoc_identifier=mode(noassoc_seg_vals)[0]
                except UnboundLocalError:
                    logging.warning('No object found in seg map....skipping')
                    continue
                #Set pixels in noassoc seg stamp with pixel value noassoc_identifier to identifier, and others to -1
                noassoc_obj_inds=np.where(noassoc_seg_stamp==noassoc_identifier)
                noassoc_seg_stamp[(noassoc_seg_stamp!=0)] = -1
                noassoc_seg_stamp[noassoc_obj_inds] = identifier
                if args.masking_type=='seg':
                    mask=seg_to_mask_basic(identifier,noassoc_seg_stamp)
                    mask_stamp=py3shape.Image(mask)
                if args.masking_type=='uberseg':
                    mask = seg_to_mask_uber(identifier,noassoc_seg_stamp)
                    mask_stamp=py3shape.Image(mask)
                if args.masking_type=='all':
                    mask=seg_to_mask_basic(identifier,noassoc_seg_stamp)
                    uberseg_mask=seg_to_mask_uber(identifier,noassoc_seg_stamp)
                    mask_stamp = py3shape.Image(mask)
                    extra_mask_stamps = [py3shape.Image(uberseg_mask),None]       
        else:
            mask_stamp=None
        #np.set_printoptions(threshold=np.nan)
        #print 'mask_stamp',mask
        options.sersics_x0_start = x0
        options.sersics_y0_start = y0
        #print 'starting coords:',x0,y0
        #Get position in stamp for starting position:
        galaxy = py3shape.Image(stamp)

        #Get psf image
        print xpos,ypos
        pos = galsim.PositionD(float(xpos), float(ypos))  #Not quite sure why I need float() here...but get galism error if not
        psf_size=(options.stamp_size+options.padding)*options.upsampling
        try:
            psf = getPSFExarray(psfexer, pos, psf_size, psf_size, upsampling=options.upsampling)
            psf_Image=py3shape.Image(psf)
        except Exception:
            logging.warning('failed to get psf for galaxy %d, skipping',identifier)
            continue

        try:
            result,model = py3shape.analyze(galaxy, py3shape.Image(psf), options, weight=weight_stamp, mask=mask_stamp, ID=identifier)
        except Exception,emsg:
            print emsg
            logging.error('im3shape failed for galaxy %d, with following error message:',identifier)
            logging.error(emsg)
            continue

        if args.masking_type=='all':
            result_nomask,_=py3shape.analyze(galaxy, psf_Image, options, weight=weight_stamp, mask=extra_mask_stamps[1], ID=identifier)
            e1_nomask,e2_nomask = result_nomask.get_params().e1,result_nomask.get_params().e2
            result_uber,_=py3shape.analyze(galaxy, psf_Image, options, weight=weight_stamp, mask=extra_mask_stamps[0], ID=identifier)
            e1_uber,e2_uber = result_uber.get_params().e1,result_uber.get_params().e2

        #Convert e's to sky coordinates...not totally sure about this function...
        local_wcs = wcs.local(image_pos=pos)
        e1_sky,e2_sky = convert_g_image2sky(local_wcs,result.get_params().e1,result.get_params().e2)

        #Measure psf properties
        psf_fwhm,psf_e1,psf_e2=(get_psf_params([psf_Image], options, radialProfile=True, hsm=False))[0]

        if args.plot:   
            pylab.subplot(231)
            pylab.imshow(stamp,origin='lower',interpolation='nearest')
            pylab.subplot(232)
            pylab.imshow(model,origin='lower',interpolation='nearest')
            pylab.subplot(233)
            pylab.imshow(stamp-(model)*stamp.sum(),origin='lower',interpolation='nearest')
            pylab.subplot(234)
            pylab.imshow(mask,origin='lower',interpolation='nearest')
            if args.masking_type=='all':
                pylab.subplot(235)
                pylab.imshow(uberseg_mask,origin='lower',interpolation='nearest')            

            pylab.show() 
        
        extra_output=[e1_sky,e2_sky]
        if args.psf_props:
            extra_output+=[psf_fwhm,psf_e1,psf_e2]
        if args.masking_type=='all':
            extra_output+=[e1_nomask,e2_nomask,e1_uber,e2_uber]
        output.write_row(result, options, psf, extra_output=extra_output,ra=RA[i],dec=DEC[i])
        
    total_time = time.clock() - start_time
    ngal = args.last+1 - args.first
    logging.info("Total time for processing in Python:" , total_time)
    logging.info("Time per galaxy:" , total_time / ngal)

# Set up and parse the command line arguments
description = 'Im3shape measures the shapes of galaxies in astronomical survey images,\n \
taking into account that they have been distorted by a point-spread function.\n \
For more info visit https://bitbucket.org/joezuntz/im3shape\n\n \
Standard usage:\n \
image_file object_catalogue_file psf_file output_filename\n \
[first_object_to_process] [last_object_to_process] [additional options]'

parser = argparse.ArgumentParser(description=description, add_help=True,
                 formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument('img_filename', type=str, help='image file, with extension as filename[<ext>], otherwise defaults to 0')
parser.add_argument('cat_filename', type=str, help='object catalogue file, with extension as filename[<ext>], otherwise defaults to 2')
parser.add_argument('psf_filename', type=str, help='psfex file')
parser.add_argument('out_filename', type=str, help='output filename')
parser.add_argument('--first','-f', type=int, default=0, help='[first_image_to_process]')
parser.add_argument('--last', '-l', type=int, default=None, help='[last_image_to_process]')
parser.add_argument('--ini_filename', type=str, default=DEFAULT.IM3_OPTIONS, help='options filename')
parser.add_argument('--weight_filename', type=str, help='weight file, with extension as filename[<ext>], otherwise defaults to 0')
parser.add_argument('--seg_filenames', type=str, nargs='*', help="1 or 2 segmentation mask files (from SExtractor), \
with extension as filename[<ext>], otherwise defaults to 0. If only one, assume it is that from the original image, \
in which case all nonzero pixels are ignored in fit.If two seg files, first should be from simulated image full run, second \
should be from simulated image in assoc mode (allowing object matching/uberseg).")
parser.add_argument('--x_image_col', type=str, default='X_IMAGE', help='name of x column in catalog file, defaults to X_IMAGE')
parser.add_argument('--y_image_col', type=str, default='Y_IMAGE', help='name of y column in catalog file, defaults to X_IMAGE')
parser.add_argument('--id_col', type=str, default='NUMBER', help='name of id column in catalog file, defaults to NUMBER')
parser.add_argument('--ra_col', type=str, default='ALPHAPEAK_J2000', help='name of ra column in catalog file, defaults to ALPHAPEAK_J2000')
parser.add_argument('--dec_col', type=str, default='DELTAPEAK_J2000', help='name of dec column in catalog file, defaults to DELTAPEAK_J2000')
parser.add_argument('--log_file', type=str, default=None, help='name of log file, if not specified, no log file written')
parser.add_argument('--loglevel', type=str, default='INFO', help='python logging level (DEBUG,INFO,WARNING,ERROR...see https://docs.python.org/2/howto/logging.html#)')
parser.add_argument('--masking_type', type=str, default='seg',help="""One of 'seg' (neighbouring object segmentaiton pixels weighted to zero), 
                    'uberseg','none' or 'all'. Need two seg maps to use uberseg""")
parser.add_argument('--psf_props','-pp', action='store_true',help="""Save psf fwhm and ellipticities""")
parser.add_argument('--plot', action='store_true', default=False, help='display image, model and residuals for each galaxy (using pylab.imshow())')
parser.add_argument('-p', '--option', dest='extra_options', action='append',
    help='Additional options to be set. Can specify more than once in form -p option=value. Overrides ini file')

    
if __name__ == "__main__":
    args = parser.parse_args()
    if args.log_file:
        logging.basicConfig(filename=args.log_file,level=getattr(logging, args.loglevel.upper()))
    if args.plot:
        import pylab
    if args.masking_type not in ['seg','uberseg','none','all']:
        print "invalid masking type, should be one of 'seg','uberseg','none','all', exiting"
        logging.error("invalid masking type, should be one of 'seg','uberseg','none','all', exiting")
        exit(1)
    if args.masking_type==('uberseg' or 'all') and len(args.seg_filenames) < 2:
        print "Need two seg maps to use uberseg, reverting to default 'seg' masking"
        logging.warning("Need two seg maps to use uberseg, reverting to default 'seg' masking")
        args.masking_type='seg'
    if args.masking_type=='seg' and not args.seg_filenames:
        print "Need a seg map to do 'seg' masking, but none provided, proceeding without masking"
        logging.warning("Need a seg map to do 'seg' masking, but none provided, proceeding without masking")
        args.masking_type='none'    
    main(args)