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


### How you want to simulate your galaxies
def SimulationRules(args, rules, sampled):
    cat = args.catalogsample
    ext = args.ext

    # Simulated galaxies only have one of each of these
    rules.x = Random(args.xmin, args.xmax)
    rules.y = Random(min=args.ymin, max=args.ymax)
    rules.g1 = 0
    #rules.g1 = Value(0)  # This line does exactly the same thing as the one above
    rules.g2 = sampled.g1
    #rules.g2 = Same('g1')# This line does exactly the same thing as the one above
    rules.magnification = np.ones(args.ngal)
    #rules.magnification = Array( np.ones(args.ngal) ) # This line does exactly the same thing as the one above
    
    # Simulated galaxies can have as many Sersic Profiles as you want. Make an array element for each.
    # Being precise, halflightradius is along the major axis (this is what sextractor measurses...I think)
    '''
    rules.nProfiles = 1
    rules.beta = [Random(-90, 90) ]
    rules.halflightradius = [Catalog(cat,ext,args.reff)]
    rules.magnitude = [Catalog(cat,ext,args.mag)]
    rules.sersicindex = [Catalog(cat,ext,args.sersicindex)]
    rules.axisratio = [Random(0.01, 1.0)]
    '''

    rules.nProfiles = 2
    rules.beta = [sampled.beta[1], 0]
    rules.halflightradius = [Gaussian(1.0, 0.1), Gaussian(avg=0.5, std=0.05)]
    rules.magnitude = [Catalog(cat,ext,args.mag), sampled.magnitude[0]]
    #rules.magnitude = [13, sampled.magnitude[0]]
    #rules.magnitude = [20, Same(0)] # This line does exactly the same thing as the one above
    rules.sersicindex = [1, 4]
    axisratio = Function(function=SampleFunction, args=(sampled.x, sampled.y, args.xmax, args.ymax))
    rules.axisratio = [axisratio, sampled.axisratio[0]]
    #rules.axisratio = [axisratio, Same(0,'axisratio')]


def SampleFunction(x,y, xmax,ymax):
    dist = np.sqrt(x*x + y*y)
    max = np.sqrt(xmax*xmax + ymax*ymax)
    return dist/max


# These are extra configurations to give to sextractor which will override the ones in the config file
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'NONE'

