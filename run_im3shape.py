import py3shape
import numpy as np
import galsim
import galsim.des
import os
import time
import pyfits
import argparse
import logging

#im3shape defaults
class DEFAULT:
    BALROG_DIR = '/'.join( os.path.realpath(__file__).split('/')[:-1] )
    ASTRO_DIR = os.path.join( BALROG_DIR, 'astro_config' )
    IM3_OPTIONS = os.path.join( ASTRO_DIR, 'im3_options.ini' )

#creat galaxy catalog dictionary from fits catalog
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
        print 'stamp over edge of image, skipping this object'
        return 1
    subBounds = galsim.BoundsI(x_min,x_max-1,y_min,y_max-1)
    stamp = image[subBounds]
    #stamp = img_data[ypos-half:ypos+half, xpos-half:xpos+half]
    x0,y0 = x_pos - x_min, y_pos - y_min
    return stamp, x0, y0



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
    ext=0
    if args.img_filename[-1]==']' and args.img_filename[-3]=='[':
        try:
            ext=int(args.img_filename[-2])
            args.img_filename = args.img_filename[:-3]
        except TypeError:
            print 'specify image as fits filename with extension appended in square brackets'
            print 'e.g. bullet.fits[2]'
            raise
    
    #Read in image
    image_gs=galsim.fits.read(args.img_filename,hdu=ext)
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

    output.write_header(include_radec=False, extra_cols=extra_cols,extra_lines=extra_lines)

    #Write extra lines to info() log
    for extra_line in extra_lines:
        logging.info(extra_line)

    # Loop over objects in catalog
    logging.info('Analyzing %d/%d galaxies' % (min(args.last-args.first+1,len(gal_cat['ID'])), len(gal_cat['ID'])))
    start_time = time.clock()

    #main galaxy loop
    ID=gal_cat['ID']
    X_IMAGE=gal_cat['X_IMAGE']
    Y_IMAGE=gal_cat['Y_IMAGE']

    for i,identifier in enumerate(ID):
        if i>args.last:
            logging.info('Done')
            break
        # Read galaxy catalog entry
        identifier = ID[i]
        xpos = X_IMAGE[i]
        ypos = Y_IMAGE[i]
        stamp_stuff = get_stamp(image_gs,X_IMAGE[i],Y_IMAGE[i],options.stamp_size)
        try:
            stamp,x0,y0 = stamp_stuff
        except TypeError:
            if stamp_stuff==1:
                logging.warning('galaxy %d stamp overlaps edge of image, skipping',identifier)
            else:
                logging.warning('failed to get stamp for galaxy %d, skipping',identifier)
            continue

            
        stamp_array = stamp.array
        options.sersics_x0_start = x0
        options.sersics_y0_start = y0
        #print 'starting coords:',x0,y0
        #Get position in stamp for starting position:
        galaxy = py3shape.Image(stamp_array)

        #Get psf image
        print xpos,ypos
        pos = galsim.PositionD(float(xpos), float(ypos))  #Not quite sure why I need float() here...but get galism error if not
        psf_size=(options.stamp_size+options.padding)*options.upsampling
        try:
            psf = getPSFExarray(psfexer, pos, psf_size, psf_size, upsampling=options.upsampling)
        except Exception:
            logging.warning('failed to get psf for galaxy %d, skipping',identifier)
            continue

        try:
            result,model = py3shape.analyze(galaxy, py3shape.Image(psf), options, ID=identifier)
        except Exception,emsg:
            logging.error('im3shape failed for galaxy %d, with following error message:',identifier)
            logging.error(emsg)
            continue

        #Convert e's to sky coordinates...not totally sure about this function...
        local_wcs = wcs.local(image_pos=pos)
        e1_sky,e2_sky = convert_g_image2sky(local_wcs,result.get_params().e1,result.get_params().e2)

        if args.plot:   
            pylab.subplot(131)
            pylab.imshow(stamp_array,origin='lower',interpolation='nearest')
            pylab.subplot(132)
            pylab.imshow(model,origin='lower',interpolation='nearest')
            pylab.subplot(133)
            pylab.imshow(stamp_array-(model)*stamp_array.sum(),origin='lower',interpolation='nearest')
            pylab.show()
        
        extra_output=[e1_sky,e2_sky]
        output.write_row(result, options, psf, extra_output=extra_output)
        # Compute rgpp/rp
        # rgpp_rp = compute_rgpp_rp(result, psf, options)
        
        # Write out results within python rather than from the C side
        
        
    total_time = time.clock() - start_time
    ngal = args.last+1 - args.first
    print "Total time for processing in Python:" , total_time
    print "Time per galaxy:" , total_time / ngal

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
parser.add_argument('--first', type=int, default=0, help='[first_image_to_process]')
parser.add_argument('--last', type=int, default=5, help='[last_image_to_process]')
parser.add_argument('--ini_filename', type=str, default=DEFAULT.IM3_OPTIONS, help='options filename')
parser.add_argument('--weight_filename', type=str, help='weight file, with extension as filename[<ext>], otherwise defaults to 1')
parser.add_argument('--x_image_col', type=str, default='X_IMAGE', help='name of x column in catalog file, defaults to X_IMAGE')
parser.add_argument('--y_image_col', type=str, default='Y_IMAGE', help='name of y column in catalog file, defaults to X_IMAGE')
parser.add_argument('--id_col', type=str, default='NUMBER', help='name of id column in catalog file, defaults to NUMBER')
parser.add_argument('--log_file', type=str, default=None, help='name of log file, if not specified, no log file written')
parser.add_argument('--loglevel', type=str, default='INFO', help='python logging level (DEBUG,INFO,WARNING,ERROR...see https://docs.python.org/2/howto/logging.html#)')
parser.add_argument('--plot', action='store_true', default=False, help='display image, model and residuals for each galaxy (using pylab.imshow())')
parser.add_argument('-p', '--option', dest='extra_options', action='append',
    help='Additional options to be set. Can specify more than once in form -p option=value. Overrides ini file')

        

if __name__ == "__main__":
    args = parser.parse_args()
    if args.log_file:
        logging.basicConfig(filename=args.log_file,level=getattr(logging, args.loglevel.upper()))
    if args.plot:
        import pylab
    main(args)