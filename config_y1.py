#!/usr/bin/env python

import os
import numpy as np
import pywcs
import pyfits
from model_class import *

thisdir = os.path.dirname( os.path.realpath(__file__) )

### In this function you can define your own command line arguments.
### Google python argparse for help if the syntax is unclear.
def CustomArgs(parser):
    # Catlog to sample simulation parameters from
    parser.add_argument( "-cs", "--catalog", help="Catalog used to sample simulated galaxy parameter distriubtions from",
                         type=str, default=os.path.join(thisdir,"CMC_allband_upsample.fits"))
    parser.add_argument( "-ext", "--ext", help="Index of the data extension for sampling catalog", type=int, default=1)
    parser.add_argument( "-mag", "--mag", help="Column name when drawing magnitude from catalog", type=str, default=None)


def CustomParseArgs(args):
    if args.mag==None:
        args.mag = 'MAPP_%s_SUBARU' %(args.band.upper())
    

def SimulationRules(args, rules, sampled, TruthCat):

    rules.x = Function(function=rand, args=[args.xmin, args.xmax, args.ngal])
    rules.y = Function(function=rand, args=[args.ymin, args.ymax, args.ngal])
    rules.g1 = 0.
    rules.g2 = 0.
    rules.magnification = np.ones(args.ngal)
    
    # Being precise, halflightradius is along the major axis
    tab = Table(file=args.catalog, ext=args.ext)
    rules.halflightradius = tab.Column("halflightradius")
    rules.magnitude = tab.Column(args.mag)
    rules.sersicindex = 1. #tab.Column("sersicindex")
    rules.beta = Function(function=rand, args=[0., 360., args.ngal])
    rules.axisratio = tab.Column("axisratio")

### Adjust the galsim GSParams
def GalsimParams(args, gsparams, galaxies):
    gsparams.alias_threshold = 1e-3
    gsparams.maximum_fft_size = 8192

# Each dictionary keyword can be one of the sextractor config file keywords
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'SEGMENTATION'

# Extra functions the user has defined. Could be used with sampling type Function
def rand(minimum, maximum, ngal):
    return np.random.uniform( minimum, maximum, ngal )

def SampleFunction(x, y, xmax=1000, ymax=1000):
    dist = np.sqrt(x*x + y*y)
    max = np.sqrt(xmax*xmax + ymax*ymax)
    return dist/max

