#!/usr/bin/env python

import os
import numpy as np
from model_class import *


###In lieu of having time for proper documentation, for now I'll give you some commented examples, which you can extrapolate from.


### Your own command line arguments
### Google python argparse for help
def CustomArguments(parser):
    # Catlog to sample simulation parameters from
    parser.add_argument( "-cs", "--catalogsample", help="Catalog used to sample simulated galaxy parameter distriubtions from", type=str, default=None)
    parser.add_argument( "-ext", "--ext", help="Index of the data extension for sampling catalog", type=int, default=1)
    parser.add_argument( "-reff", "--reff", help="Column name when drawing half light radius from catalog", type=str, default="HALF_LIGHT_RADIUS")
    parser.add_argument( "-nsersic", "--sersicindex", help="Column name when drawing sersic index catalog", type=str, default="SERSIC_INDEX")
    parser.add_argument( "-mag", "--mag", help="Column name when drawing magnitude from catalog", type=str, default=None)
    parser.add_argument( "-b", "--band", help="Which filter band to choose from COSMOS catalog. Only relevant if --mag is not given and using COSMOS catlalog.", type=str, default='i', choices=['g','r','i','z'])


### How you want to parse your command line arguments
def CustomParseArgs(args):
    thisdir = os.path.dirname( os.path.realpath(__file__) )

    if args.catalogsample==None:
        #args.catalogsample = os.path.join(thisdir, 'cosmos_n=1.fits')
        args.catalogsample = os.path.join(thisdir, 'cosmos_n=4.fits')

    if args.mag==None:
        args.mag = '%sMAG' %(args.band.upper())


def SimulationRules(args, rules):
    cat = args.catalogsample
    ext = args.ext

    # Simulated galaxies only have one of each of these
    rules.x = Random(args.xmin, args.ymax)
    rules.y = Random(args.xmin, args.ymax)
    rules.g1 = Value(0)
    rules.g2 = Same('g1')
    rules.magnification = Array( np.ones(args.ngal) )
    #rules.magnification = Gaussian( 1, 0.05 )
    
    # Simulated galaxies can have as many Sersic Profiles as you want. Make an array element for each.
    # Being precise, halflightradius is along the major axis (this is what sextractor measurses...I think)
    rules.nProfiles = 1
    rules.axisratio = [Random(0.33, 1)]
    rules.beta = [Random(-90, 90) ]
    rules.halflightradius = [Catalog(cat,ext,args.reff)]
    #rules.magnitude = [Catalog(cat,ext,args.mag)]
    rules.magnitude = [Value(18)]
    rules.sersicindex = [Catalog(cat,ext,args.sersicindex)]

    '''
    cat = 'SomethingYouMade.fits'
    ext = 0
    nProfiles = 2
    axisratio = [Random(0.33, 1), Same(0)]
    #axisratio = [Random(0.33, 1), Same( (0,'axisratio') )]
    beta = [Random(-90, 90), Same(0)]
    halflightradius = [Catalog(cat,ext,'DISKSIZE'), Catalog(cat,ext,'BULGESIZE')]
    halflightradius = [Same(1), Catalog(cat,ext,'BULGESIZE')]
    magnitude = [Catalog(cat,ext,'DISKMAG'), Catalog(cat,ext,'BULGEMAG')]
    sersicindex = [Value(1), Value(4)]
    '''

# These are extra configurations to give to sextractor which will override the ones in the config file
def SextractorConfigs(args, config):
    config['CHECKIMAGE_TYPE'] = 'NONE'

