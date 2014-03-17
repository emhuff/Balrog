#!/usr/bin/env python

import os
import numpy as np
from model_class import *


### In this function you can define, your own command line arguments.
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
### Your custom command line arguments will be logged as they exist at the end of this file.

### You also have local access to the native Balrog command line arguments. 
### However, to avoid accidentally breaking Balrog, any changes you make to these native args do not propagate outside this file.
###         e.g. You could change args.xmin for convenience if you wanted, but this would have no effect on the minimum x for your simulation area.
### Any default arguments you didn't specficy from the command line have already assumed their default values when these functions are called.


### How you want to parse your command line arguments
def CustomParseArgs(args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )

    if args.catalogsample==None:
        #args.catalogsample = os.path.join(thisdir, 'cosmos_n=1.fits')
        #args.catalogsample = os.path.join(thisdir, 'cosmos_n=4.fits')
        args.catalogsample = os.path.join(thisdir, 'cosmos.fits')

    if args.mag==None:
        args.mag = '%sMAG' %(args.band.upper())
    
    '''
    args.catalogsample = os.path.join(thisdir, 'cosmos.fits')
    args.mag = 'IMAG'
    args.sersicindex = 'SERSIC_INDEX'
    args.reff = 'HALF_LIGHT_RADIUS'
    args.ext = 1
    '''


### How you want to simulate your galaxies
def SimulationRules(args, rules, sampled):
    '''
    thisdir = os.path.dirname( os.path.realpath(__file__) )
    args.catalogsample = os.path.join(thisdir, 'cosmos.fits')
    args.mag = 'IMAG'
    args.sersicindex = 'SERSIC_INDEX'
    args.reff = 'HALF_LIGHT_RADIUS'
    args.ext = 1
    '''

    cat = args.catalogsample
    ext = args.ext


    # Simulated galaxies only have one of each of these
    rules.x = Function(function=rand, args=(args.xmin, args.xmax, args.ngal))
    rules.y = Function(function=rand, args=(args.ymin, args.ymax, args.ngal))
    rules.g1 = 0
    rules.g2 = 0
    rules.magnification = 1
    
    # Simulated galaxies can have as many Sersic Profiles as you want. Make an array element for each.
    # Being precise, halflightradius is along the major axis (this is what sextractor measurses...I think)
    # All 3 of these blocks (2 of which are commented out) work fine.
    '''
    InitializeSersic(rules, sampled, nProfiles=1)
    rules.beta[0] = Function(function=rand, args=(-90, 90, args.ngal))
    rules.halflightradius[0] = Catalog(file=cat,ext=ext,col=args.reff)
    rules.magnitude[0] = Catalog(cat,ext,args.mag)
    rules.sersicindex[0] = Catalog(cat,ext,args.sersicindex)
    rules.axisratio[0] = Function(function=rand, args=(0.01, 1.0, args.ngal))
    #rules.axisratio[0] = Function(function=SampleFunction, args=(sampled.x, sampled.y, args.xmax, args.ymax))
    '''

    # If you have a single Sersic component, you don't need to use the indexing if you don't want. With more than one component the indexing is mandatory
    InitializeSersic(rules, sampled, nProfiles=1)
    rules.beta = Function(function=rand, args=(-90, 90, args.ngal))
    rules.halflightradius = Catalog(file=cat,ext=ext,col=args.reff)
    rules.magnitude = Catalog(cat,ext,args.mag)
    rules.sersicindex = Catalog(cat,ext,args.sersicindex)
    #rules.axisratio = Function(function=rand, args=(0.01, 1.0, args.ngal))
    rules.axisratio = Function(function=SampleFunction, args=(sampled.x, sampled.y, args.xmax, args.ymax))
    #rules.halflightradius = 2
    #rules.magnitude = 17
    #rules.sersicindex = 1
    #rules.axisratio = Function(function=SampleFunction, args=(sampled.x, sampled.y, args.xmax, args.ymax))

    '''
    InitializeSersic(rules, sampled, nProfiles=2)
    rules.beta[1] = 0
    rules.beta[0] = sampled.beta[1]
    rules.halflightradius = [Catalog(cat,ext,args.reff), sampled.halflightradius[0]]
    rules.magnitude = [Catalog(cat,ext,args.mag), sampled.magnitude[0]]
    ns = Function(function=f, args=(np.ones(args.ngal)))
    nc = Catalog(file='cosmos_n=1.fits',ext=1,col='SERSIC_INDEX')
    # Both of these work
    #n = Function(function=g, args=(4, 0.05, args.ngal, ns))
    n = Function(function=g, args=(4, 0.05, args.ngal, nc))
    rules.sersicindex = [1, n]
    #rules.sersicindex = [1]
    axisratio = Function(function=SampleFunction, args=(sampled.x, sampled.y, args.xmax, args.ymax))
    rules.axisratio = [axisratio, sampled.axisratio[0]]
    '''


# These are extra configurations to give to sextractor which will override the ones in the config file
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'NONE'



# Extra functions the user has defined. Could be used with sampling type Function
def f(item):
    return item

def g(avg, std, ngal, other):
    gg = gaussian(avg, std, ngal)
    return gg-other

def rand(minimum, maximum, ngal):
    return np.random.uniform( minimum, maximum, ngal )

def gaussian(avg, std, ngal):
    return np.random.normal( avg, std, ngal )

def SampleFunction(x,y, xmax,ymax):
    dist = np.sqrt(x*x + y*y)
    max = np.sqrt(xmax*xmax + ymax*ymax)
    return dist/max
