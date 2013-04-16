#!/usr/bin/env python

import subprocess
import os
import sys
import pyfits
import math
import argparse

class DEFAULT:
    BALROG_DIR = '/'.join( os.path.realpath(__file__).split('/')[:-1] )

    ASTRO_DIR = os.path.join( BALROG_DIR, 'astro_config' )
    SEX_CONFIG = os.path.join( ASTRO_DIR, 'sex.config' )
    SEX_PARAM = os.path.join( ASTRO_DIR, 'sex.param' )
    SEX_NNW = os.path.join( ASTRO_DIR, 'sex.nnw' )
    SEX_CONV = os.path.join( ASTRO_DIR, 'sex.conv' )

    LABEL = 'intermediate'
    DIR_LABEL = 'default'
    WORKING_DIR = os.path.join( BALROG_DIR, 'intermediate_files', DIR_LABEL )
    #OUTPUT_DIR = os.getcwd()


class Parser():

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument( "-i", "--image", help="Real image to insert simulated galaxies into", type=str, default="None" )
        self.parser.add_argument( "-w", "--weight", help="Real weight map to insert simulated galaxies into", type=str, default="None" )
        #self.parser.add_argument( "-o", "--output_dir", help="Output directory where final results are written", type=str, default=DEFAULT.OUTPUT_DIR )
        self.parser.add_argument( "-wd", "--working_dir", help="Working directory where intermdediate files are written", type=str, default=DEFAULT.WORKING_DIR )
        self.parser.add_argument( "-l", "--label", help="Label prepended to intermediate files", type=str, default=DEFAULT.LABEL )

        self.parser.add_argument( "-bd", "--balrog_dir", help="Balrog HOME directory", type=str, default=DEFAULT.BALROG_DIR )
        self.parser.add_argument( "-ad", "--astro_dir", help="Astronomy software configuration directory", type=str, default=DEFAULT.ASTRO_DIR )
        self.parser.add_argument( "-conf", "--sex_config", help="Sextractor configuration file", type=str, default=DEFAULT.SEX_CONFIG )
        self.parser.add_argument( "-param", "--sex_param", help="Sextractor parameters file", type=str, default=DEFAULT.SEX_PARAM )
        self.parser.add_argument( "-nnw", "--sex_nnw", help="Sextractor star neural network file", type=str, default=DEFAULT.SEX_NNW)
        self.parser.add_argument( "-conv", "--sex_conv", help="Sextractor convolution file", type=str, default=DEFAULT.SEX_CONV)

        self.args = self.parser.parse_args()
        self.HandleErrors()

    class OptionError(Exception):
        def __init__(self,value):
            self.value = value
        def __str__(self):
            return repr(self.value)

    def HandleErrors(self):
        if self.args.image=='None':
            raise self.OptionError('No image file specfied.')
        if self.args.weight=='None':
            raise self.OptionError('No weight map file specfied.')



if __name__=='__main__':

    parser = Parser()
    args = parser.args

    im_str = args.image + ',' + args.image
    weight_str = args.weight + ',' + args.weight
    out_cat = os.path.join( args.working_dir, '%s_cat.fits' %args.label )
    out_seg = os.path.join( args.working_dir, '%s_seg.fits' %args.label )

    chain = args.working_dir.split('/')
    working_path = '/'
    for sub in chain:
        working_path = os.path.join(working_path, sub)
        if not os.path.exists(working_path):
            subprocess.call( ['mkdir', working_path] )
        
    #if not os.path.exists(args.working_dir):
        #subprocess.call( ['mkdir', args.working_dir] )
   
    log_file = os.path.join( args.working_dir, '%s_log.txt' %args.label )
    log = open(log_file, 'w')

    subprocess.call( ['sex', im_str, 
                      '-c', args.sex_config, 
                      '-WEIGHT_IMAGE', weight_str, 
                      '-CATALOG_NAME', out_cat, 
                      '-CHECKIMAGE_NAME', out_seg, 
                      '-PARAMETERS_NAME', args.sex_param, 
                      '-FILTER_NAME', args.sex_conv, 
                      '-STARNNW_NAME', args.sex_nnw],
                    stdout=log, stderr=log )
    log.close()
    #sys.exit('DONE')

    """
    if filter!='Y':
        im_hdus = pyfits.open(fits_file)
        im_header = im_hdus[0].header
        fzp = im_header['FID_ZP']
        im_hdus.close()
        zp = -fzp - 25.0 + ( 2.5 * math.log10(90.0) )

    else:
        zp = 3.05
        if (clus=='elgordo') and (code in ['1','3','5','11']):
            zp = zp + ( 2.5 * math.log10(90.0/50.0) )

    cat_hdus = pyfits.open(link_cat, mode='update')
    cat_header = cat_hdus[2].header
    cat_header['ZPOFFSET'] = zp
    cat_hdus.flush()
    cat_hdus.close()
    """
