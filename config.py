#!/usr/bin/env python

import os
import numpy as np
from model_class import *
import array


### In this function you can define your own command line arguments.
### Google python argparse for help if the syntax is unclear.
def CustomArgs(parser):
    # Catlog to sample simulation parameters from
    parser.add_argument( "-cs", "--catalogsample", help="Catalog used to sample simulated galaxy parameter distriubtions from", type=str, default=None)
    parser.add_argument( "-ext", "--ext", help="Index of the data extension for sampling catalog", type=int, default=1)
    parser.add_argument( "-reff", "--reff", help="Column name when drawing half light radius from catalog", type=str, default="HALF_LIGHT_RADIUS")
    parser.add_argument( "-nsersic", "--sersicindex", help="Column name when drawing sersic index catalog", type=str, default="SERSIC_INDEX")
    parser.add_argument( "-mag", "--mag", help="Column name when drawing magnitude from catalog", type=str, default=None)
    parser.add_argument( "-b", "--band", help="Which filter band to choose from COSMOS catalog. Only relevant if --mag is not given and using COSMOS catlalog.", type=str, default='i', choices=['g','r','i','z'])


### Throughout the remainder of the file, you have access to your custom command line arguments in the attributes of args.
### The three functions below execute in the order they appear. Attributes changed in ealier functions will propogate downstream.
### Your custom command line arguments will be logged as they exist at the end of the three functions.

### You also have local access to the native Balrog command line arguments. 
### However, to avoid accidentally breaking Balrog, any changes you make to these native args do not propagate outside this file.
###         e.g. You could change args.xmin for convenience if you wanted, but this would have no effect on the minimum x for your simulation area.
### Any native arguments you didn't specficy from the command line have already assumed their default values when these functions are called.

### How you want to parse your command line arguments
def CustomParseArgs(args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    if args.catalogsample==None:
        args.catalogsample = os.path.join(thisdir, 'cosmos.fits')
    if args.mag==None:
        args.mag = '%sMAG' %(args.band.upper())
    

### Each galaxy has a single {x, y, g1, g2, magnification}.
### Each galaxy can have as many sersic components as you want,
###    i.e {beta, axisratio, sersicindex, halflightradius, magnitude} 
###    can be arrays. This example shows the simplest case, using the
###    default of a single component.

### You can simulate parameters in 5 different ways:
###    The same constant for each galaxy --> rules.g1 = 0
###    An array with an element for each galaxy --> rules.g1 = np.zeros(args.ngal)
###    The same as another parameter --> rules.g2 = sampled.g1
###    Sample drawing from a FITS catalog --> rules.magnitude = Catalog(file='somefitscatalog.fits', ext=1, col='IMAG')
###    A function you defined --> rules.x = Function(function=rand, args=[args.xmin, args.xmax, args.ngal])

### Functions used with sampling type Function MUST return an array of length ngal
### Function accepts *args or **kwargs along the lines of usual python syntax,
###    where args is required arguments as a list and kwargs is a dictionary of keyword arguments.

### rules attributes are to be reassigned to define the simulation.
### sampled gives you access to what the parameters ultimately sample to. Sampled attributes cannot be reassigned.
### sampled attributes can be used as arguments to Function.
### Catalog arguments may also be used as Function arguments.
### Nesting a Function argument within another Function agument is also permissible.

### How you want to simulate your galaxies
def SimulationRules(args, rules, sampled):
    cat = args.catalogsample
    ext = args.ext

    rules.x = Function(function=rand, args=[args.xmin, args.xmax, args.ngal])
    rules.y = Function(function=rand, args=[args.ymin, args.ymax, args.ngal])

    rules.g1 = 0
    rules.g2 = 0
    rules.magnification = np.ones(args.ngal)
    
    # Being precise, halflightradius is along the major axis
    rules.halflightradius = Catalog(file=cat,ext=ext,col=args.reff)
    #rules.magnitude = Catalog(cat,ext,args.mag)
    rules.magnitude = 14
    rules.sersicindex = Catalog(cat,ext,args.sersicindex)

    rules.beta = Function(function=rand, args=[-90, 90, args.ngal])
    rules.axisratio = Function(function=rand, args=[0.05, 1, args.ngal])
    #rules.axisratio = Function(function=SampleFunction, args=[sampled.x, sampled.y, args.xmax, args.ymax])
    #rules.axisratio = Function(function=SampleFunction, args=[sampled.x, sampled.y], kwargs={'xmax':args.xmax, 'ymax':args.ymax})


### Adjust the galsim GSParams
def GalsimParams(args, gsparams, galaxies):
    gsparams.alias_threshold = 1e-3
    #gsparams.alias_threshold = Function( function=StupidSize, args=[galaxies.halflightradius] )



# These are extra configurations to give to sextractor which will override the ones in the config file
# Each dictionary keyword can be one of the sextractor config file keywords
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'NONE'



# Extra functions the user has defined. Could be used with sampling type Function
def rand(minimum, maximum, ngal):
    return np.random.uniform( minimum, maximum, ngal )

def SampleFunction(x, y, xmax=1000, ymax=1000):
    dist = np.sqrt(x*x + y*y)
    max = np.sqrt(xmax*xmax + ymax*ymax)
    return dist/max

def StupidSize(size):
    return size * 1.0e-3
