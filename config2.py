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
### Only changes to args made in CustomParseArgs persist throughout this file.

### You also have local access to the native Balrog command line arguments in args.
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
###    can be arrays. This example shows the case of 2 components.
###    Multi-component models require InitializeSeric to be called
###    prior to using {beta, axisratio, sersicindex, halflightradius, magnitude}.

### You can simulate parameters in 4 (5 depending how you count) different ways:
###    1)   The same constant for each galaxy --> rules.g1 = 0
###    2)   An array with an element for each galaxy --> rules.g1 = np.zeros(args.ngal)
###    2.5) The same as another parameter --> rules.g2 = sampled.g1
###    3)   Sample drawing from a FITS catalog --> rules.magnitude[0] = Catalog(file='somefitscatalog.fits', ext=1, col='IMAG')
###    4)   A function you defined --> rules.x = Function(function=rand, args=[args.xmin, args.xmax, args.ngal])

### Functions used with sampling type Function MUST return an array of length ngal
### Function accepts *args or **kwargs along the lines of usual python syntax,
###    where args is required arguments as a list and kwargs is a dictionary of keyword arguments.

### rules attributes are to be reassigned to define the simulation.
### Either the entire array or individual components can be reassigned for sersic components:
###    rules.beta[0] = sampled.beta[1]
###    rules.beta[1] = 0
###    rules.magnitude = [Catalog(cat,ext,args.mag), sampled.magnitude[0]]

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
    rules.magnification = 1
    
    InitializeSersic(rules, sampled, nProfiles=2)

    rules.beta[0] = sampled.beta[1]
    rules.beta[1] = 0
    pos = [sampled.x, sampled.y]
    pos_max = [args.xmax, args.ymax]
    axisratio = Function(function=SampleFunction, args=[pos, pos_max])
    rules.axisratio = [axisratio, sampled.axisratio[0]]

    rules.halflightradius = [Catalog(cat,ext,args.reff), sampled.halflightradius[0]]
    rules.magnitude = [Catalog(cat,ext,args.mag), sampled.magnitude[0]]

    # Both of these work
    '''
    nc = Catalog(file='cosmos_n=1.fits',ext=1,col='SERSIC_INDEX')
    n = Function(function=g, args=[4, 0.05, args.ngal, nc])
    rules.sersicindex = [1, n]
    '''
    ns = Function(function=f, args=[np.ones(args.ngal)])
    n = Function(function=g, args=[4, 0.05, args.ngal, ns])
    rules.sersicindex = [1, n]


# These are extra configurations to give to sextractor which will override the ones in the config file
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'NONE'


# Extra functions the user has defined. Could be used with sampling type Function
def rand(minimum, maximum, ngal):
    return np.random.uniform( minimum, maximum, ngal )

def f(item):
    return item

def g(avg, std, ngal, other):
    gg = gaussian(avg, std, ngal)
    return gg-other

def gaussian(avg, std, ngal):
    return np.random.normal( avg, std, ngal )

def SampleFunction(pos, pos_max=[1000,1000]):
    dist = np.sqrt(pos[0]*pos[0] + pos[1]*pos[1])
    max = np.sqrt(pos_max[0]*pos_max[0] + pos_max[1]*pos_max[1])
    return dist/max
