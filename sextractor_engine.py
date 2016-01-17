#!/usr/bin/env python

#import subprocess
#import sys
#import time
import os
import balrog
import copy

#import utils
#from defaults import *

class DEFAULT:
    IMAGE = 'NONE'
    WEIGHT = 'NONE'
    CHECK = 'NONE'
    CHECK_NAME = 'NONE'

    BALROG_DIR = '/'.join( os.path.realpath(__file__).split('/')[:-1] )
    ASTRO_DIR = os.path.join( BALROG_DIR, 'astro_config' )
    SEX_CONFIG = os.path.join( ASTRO_DIR, 'sex.config' )
    SEX_PARAM = os.path.join( ASTRO_DIR, 'sex.param' )
    SEX_NNW = os.path.join( ASTRO_DIR, 'sex.nnw' )
    SEX_CONV = os.path.join( ASTRO_DIR, 'sex.conv' )


class SextractorEngine():
    """Overview of SextractorEngine class:
        This is an object for running sextractor. Most importantly, it contains a dictionary of the configuration parameters given to sextractor.
        This dictionary is named SextractorEngine().config.  Editing the dictionary directly is one way to change the sextractor configuration.
        Values in the dictionary override those specified in the config file itself.  A default call to SextractorEngine will load default
        files for the sextractor config, param, nnw, and filter filters. These defaults are located in the astro_config directory packed with
        this file.
        
        Technically SextractorEngine() can be called with no arguments, but the sextractor run will fail if you do not supply an image and a weight file.
        The fields that can be supplied in the constructor are the sextractor config parameters which take file paths (and CHECKIMAGE_TYPE, which is related to CHECKIMAGE_NAMES).
           IMAGE 
           WEIGHT_IMAGE 
           CATALOG_NAME 
           c
           PARAMETERS_NAME
           STARNNW_NAME
           FILTER_NAME
           CHECKIMAGE_TYPE
           CHECKIMAGE_NAME
        However, any other parameters available to configure sextractor can also be modified simply by adding them to the object's config dictionary. 
       
        The image(s) to give to sextractor are normally only given on the command line (i.e. not in the config file). The dictionary keyword for it here is 'IMAGE'.
        When editing the dictionary do not use any -'s in the keywords. 
        For example when running 'sex image.fits -c sex.config -WEIGHT_IMAGE weight.fits', 
        you would have SE.config['c'] = 'sex.config', SE.config['IMAGE'] = image.fits, SE.config['WEIGHT_IMAGE'] in the config dictionary

        In the sextractor documenation sextractor keyword, value pairs in the config file are primarily upper case. I stuck with the same case conventions as sextractor, 
        so formatting should be the same as in the sextractor documentation (even though this can be slightly annoying to type).

        Methods of SextractorEngine:
            Use the run() method to run sextractor. The optional logfile argument is where to write the sextractor logging output. By default, it goes
            to the same directory as the catalog, stripping off .fits type extensions and appending _log.txt

            auto_checkimage_name() is a convenience function for automatically populating the CHECKIMAGE_NAME field, based on the CHECKIMAGE_TYPE and current
            image or catalog name. By default the output directory for the check files is the same as the catalogs, and the check files are named with the same
            base naming scheme as the image, but with a _IMAGETYPE appended to the file, e.g. image_BACKGROUND.fits.  The optional named_by argument can either be
            'image' or 'catalog' to choose if the base name is from the image or catalog file. The optional dir argument changes the output directory for
            check files.  It can be a single directory or a list the same length as the number of check files.

        A few notes:
            If you change some config parameters in the dictionary, and then change the config file referenced in the dictionary, the values in the dictionary will still take
            precedence over those in the (new) file.  Anything in the dictionary will become an explicit command line argument given to sextractor.
            To avoid this, if you would like to start fresh, with a new config file, call the constructor, e.g.
            SE = SextractorEngine( IMAGE=image.fits, WEIGHT_IMAGE=weight.fits, c=new_sex.config )
            where you could also specificy different params file, etc. if you wanted.

            Currently I do not have a way for changing the params (which will be output in the catalog), other than chaning the 'PARAMTERS_NAME' keyword.  I may add another way.

    """

    
    ############ You should not need to used any of these functions. They are things which helping under the hood.
    
    def _catalog_name(self, image, catalog_name):
        if catalog_name is None:
            self.config['CATALOG_NAME'] = os.path.join( os.getcwd(), image.strip().split('/')[-1].replace('.fits','.cat.fits') )
        else:
            self.config['CATALOG_NAME'] = catalog_name

    def _checkimage(self, checkimage_type, checkimage_name):        
        if checkimage_type!='NONE':
            self.config['CHECKIMAGE_TYPE'] = checkimage_type
            if checkimage_name!='NONE':
                self.config['CHECKIMAGE_NAME'] = checkimage_name
            else:
                self.auto_checkimage_name()


    def _strip(self, key, endings):
        name = self.config[key]
        for ending in endings:
            if self.config[key].endswith(ending):
                name = name.rstrip(ending)
        return name

    ################################################################################################3


    def __init__(self, IMAGE=DEFAULT.IMAGE, WEIGHT_IMAGE=DEFAULT.WEIGHT, CATALOG_NAME=None, c=DEFAULT.SEX_CONFIG, PARAMETERS_NAME=DEFAULT.SEX_PARAM, STARNNW_NAME=DEFAULT.SEX_NNW, FILTER_NAME=DEFAULT.SEX_CONV, CHECKIMAGE_TYPE=DEFAULT.CHECK, CHECKIMAGE_NAME=DEFAULT.CHECK_NAME, setup=None):
        #print locals()

        self.config = {}
        self.config['IMAGE'] = IMAGE
        #self.config['WEIGHT_IMAGE'] = WEIGHT_IMAGE
        self.config['c'] = c
        self.config['PARAMETERS_NAME'] = PARAMETERS_NAME
        self.config['STARNNW_NAME'] = STARNNW_NAME
        self.config['FILTER_NAME'] = FILTER_NAME

        self._catalog_name(IMAGE,CATALOG_NAME)
        self._checkimage(CHECKIMAGE_TYPE, CHECKIMAGE_NAME)

        self.path = 'sex'
        self.setup = setup


    def Path(self, path):
        self.path = path


    def auto_checkimage_name(self, dir=None, named_by='image'):
        tlist = self.config['CHECKIMAGE_TYPE'].split(',')
        for i in range(len(tlist)):
            tlist[i] = '_' + tlist[i].strip().replace('-','m') + '.fits'

        stripped_cat = self._strip('CATALOG_NAME', ['.cat.fits', '.fits'] )[-1]
        stripped_image = self._strip('IMAGE', ['.fits']).split('/')[-1]
        cat_dir = '/'.join( self.config['CATALOG_NAME'].split('/')[:-1] )
        
        if named_by=='image':
            stripped = stripped_image
        elif named_by=='catalog':
            stripped = stripped_cat

        for i in range(len(tlist)):
            file = stripped + tlist[i]
            if dir is None:
                d = cat_dir
            elif type(dir)==str:
                d = dir
            elif type(dir)==list:
                d = dir[i]
            tlist[i] = os.path.join(d,file)

        nstr = ','.join(tlist)
        self.config['CHECKIMAGE_NAME'] = nstr
    

    def run(self, msg=None):
        args = [self.path, self.config['IMAGE']]
        for key in self.config.keys():
            if key=='IMAGE':
                continue
            else:
                args.append( '-%s' %key )
            args.append( str(self.config[key]) )
   
        #cmd = ' '.join(args)
        cmd = args

        if msg is not None:
            balrog.SysInfoPrint(self.setup, msg, level='info')
        balrog.SystemCall(cmd, setup=self.setup)
