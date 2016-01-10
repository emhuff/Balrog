#!/usr/bin/env python

import os
import numpy as np
#import pywcs
import astropy.wcs as pywcs
#import pyfits
from model_class import *


### In this function you can define your own command line arguments.
### Google python argparse for help if the syntax is unclear.
def CustomArgs(parser):
    # Catlog to sample simulation parameters from
    parser.add_argument( "-cs", "--catalog", help="Catalog used to sample simulated galaxy parameter distriubtions from", type=str, default=None)
    parser.add_argument( "-ext", "--ext", help="Index of the data extension for sampling catalog", type=int, default=1)
    parser.add_argument( "-reff", "--reff", help="Column name when drawing half light radius from catalog", type=str, default="HALF_LIGHT_RADIUS")
    parser.add_argument( "-nsersic", "--sersicindex", help="Column name when drawing sersic index catalog", type=str, default="SERSIC_INDEX")
    parser.add_argument( "-mag", "--mag", help="Column name when drawing magnitude from catalog", type=str, default=None)
    parser.add_argument( "-b", "--band", help="Which filter band to choose from COSMOS catalog. Only relevant if --mag is not given and using COSMOS catlalog.", type=str, default='i', choices=['g','r','i','z'])


### Throughout the remainder of the file, you have access to your custom command line arguments in the attributes of args.
### Only changes to args made in CustomParseArgs persist throughout this file.

### You also have local access to the native Balrog command line arguments in args.
### However, to avoid accidentally breaking Balrog, any changes you make to these native args do not propagate outside this file.
###         e.g. You could change args.xmin for convenience if you wanted, but this would have no effect on the minimum x for your simulation area.
### Any native arguments you didn't specficy from the command line have already assumed their default values when these functions are called.

### How you want to parse your command line arguments
def CustomParseArgs(args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    if args.catalog==None:
        args.catalog = os.path.join(thisdir, 'cosmos.fits')
    if args.mag==None:
        args.mag = 'MAPP_%s_SUBARU' %(args.band.upper())
    

### Each galaxy has a single {x, y, g1, g2, magnification}.
### Each galaxy can have as many sersic components as you want,
###    i.e {beta, axisratio, sersicindex, halflightradius, magnitude} 
###    can be arrays. This example shows the simplest case, using the
###    default of a single component.

### You can simulate parameters in 4 (5 depending how you count) different ways:
###    1)   The same constant for each galaxy --> rules.g1 = 0
###    2)   An array with an element for each galaxy --> rules.g1 = np.zeros(args.ngal)
###    2.5) The same as another parameter --> rules.g2 = sampled.g1
###    3)   Sample drawing from a FITS catalog --> rules.magnitude = Catalog(file='somefitscatalog.fits', ext=1, col='IMAG')
###    4)   A function you defined --> rules.x = Function(function=rand, args=[args.xmin, args.xmax, args.ngal])

### Functions used with sampling type Function MUST return an array of length ngal
### Function accepts *args or **kwargs along the lines of usual python syntax,
###    where args is required arguments as a list and kwargs is a dictionary of keyword arguments.

### rules attributes are to be reassigned to define the simulation.
### sampled gives you access to what the parameters ultimately sample to. Sampled attributes cannot be reassigned.
### sampled attributes can be used as arguments to Function.
### Catalog arguments may also be used as Function arguments.
### Nesting a Function argument within another Function agument is also permissible.

### How you want to simulate your galaxies
def SimulationRules(args, rules, sampled, TruthCat):

    rules.x = Function(function=rand, args=[args.xmin, args.xmax, args.ngal])
    rules.y = Function(function=rand, args=[args.ymin, args.ymax, args.ngal])
    rules.g1 = 0
    rules.g2 = 0
    rules.magnification = np.ones(args.ngal)
    
    # Being precise, halflightradius is along the major axis
    tab = Table(file=args.catalog, ext=args.ext)
    rules.halflightradius = tab.Column(args.reff)
    rules.magnitude = Column(args.catalog, args.ext, args.mag)
    rules.sersicindex = Catalog(file=args.catalog, ext=args.ext, col=args.sersicindex)
    rules.beta = 0
    rules.axisratio = 1
    '''
    rules.halflightradius = tab.Column('halflightradius')
    rules.magnitude = Column(args.catalog, args.ext, 'Mapp_i_subaru')
    rules.sersicindex = Catalog(file=args.catalog, ext=args.ext, col='sersicindex')
    rules.beta = tab.Column('beta')
    rules.axisratio = tab.Column('axisratio')
    '''

    demo = Function(function=SampleFunction, args=[sampled.x, sampled.y])
    TruthCat.AddColumn(demo, name='demo', fmt='E')
    TruthCat.AddColumn(tab.Column('ID'))
    TruthCat.AddColumn(Column(args.catalog, args.ext, 'Z'))
    TruthCat.AddColumn(tab.Column('TYPE'), name='OBJTYPE')

    '''
    file = 'fiducialwcs.fits'
    fext = 0
    test = Function(function=Test, args=[sampled.x, sampled.y, file, fext])
    TruthCat.AddColumn(test, name='TEST')
    '''

    ### extra args/kwargs examples
    #rules.axisratio = Function(function=SampleFunction, args=[sampled.x, sampled.y, args.xmax, args.ymax])
    #rules.axisratio = Function(function=SampleFunction, args=[sampled.x, sampled.y], kwargs={'xmax':args.xmax, 'ymax':args.ymax})


def Test(x, y, file, ext):
    import pyfits
    hdu = pyfits.open(file)[ext]
    header = hdu.header
    wcs = pywcs.WCS(header)
    pcoords = np.dstack((x,y))[0]
    wcoords = wcs.wcs_pix2sky(pcoords, 1)
    r = wcoords[:,0]
    return x  


### Adjust the galsim GSParams
def GalsimParams(args, gsparams, galaxies):
    gsparams.alias_threshold = 1e-3
    gsparams.maximum_fft_size = 6144
    #gsparams.alias_threshold = Function( function=StupidSize, args=[galaxies.halflightradius] )


# These are extra configurations to give to sextractor which will override the ones in the config file
# Each dictionary keyword can be one of the sextractor config file keywords
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'NONE'
    #config['CHECKIMAGE_TYPE'] = 'SEGMENTATION'



# Extra functions the user has defined. Could be used with sampling type Function
def rand(minimum, maximum, ngal):
    return np.random.uniform( minimum, maximum, ngal )

def SampleFunction(x, y, xmax=1000, ymax=1000):
    dist = np.sqrt(x*x + y*y)
    max = np.sqrt(xmax*xmax + ymax*ymax)
    return dist/max

def StupidSize(size):
    return size * 1.0e-3

