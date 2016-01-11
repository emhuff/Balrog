#!/usr/bin/env python

#import sys
#import logging
#import traceback
#import numpy.lib.recfunctions

import numpy as np
import astropy.io.fits as pyfits
import galsim
import copy
import os
from balrogexcept import *


class nComponentSersic(object):

    def __init__(self, ngal=100, ncomp=2, galkeys=[], compkeys=[]):
        self.ngal = ngal
        self.component = [None]*ncomp
        self.componentrule = [None]*ncomp

        self.galaxy = self._InitNone(galkeys)
        self.galaxyrule = self._InitRule(self.galaxy)

        for i in range(ncomp):
            self.component[i] = self._InitNone(compkeys)
            self.componentrule[i] = self._InitRule(self.component[i])

    
    def _InitNone(self, names):
        dict = {}
        for name in names: 
            dict[name] = None
        return dict


    def _InitRule(self,comp):
        dict = {}
        for key in comp.keys():
            dict[key] = Rule()
        return dict


    def ComponentRule(self, component=0, key=None, rule=None):
        
        if key is None:
            raise Exception('must indicate which field to give a rule for')

        if rule is None:
            raise Exception('must give a rule')

        self.componentrule[component][key] = rule


    def GalaxyRule(self, key=None, rule=None):

        if key is None:
            raise Exception('must indicate which field to give a rule for')
        if rule is None:
            raise Exception('must give a rule')

        self.galaxyrule[key] = rule

   
    def ReturnValue(self, val):
        return np.array( [val]*self.ngal )

    def DoValue(self, value):
        for i in range(len(value)):
            index, key, val = value[i]
            if index > -1:
                self.component[index][key] = self.ReturnValue(val)
            else:
                self.galaxy[key] = self.ReturnValue(val)


    def DoArray(self,array):
        for i in range(len(array)):
            index, key, arr = array[i]
            if index > -1:
                self.component[index][key] = np.array( arr )
            else:
                self.galaxy[key] = np.array( arr )
    

    def DoCatalog(self, catalogs, BalrogSetup, used=None):
        '''
        if used is None:
            used = []
        if len(used) > 0:
            pass
        '''
        if BalrogSetup.ngal==1:
            dtype = [ ('file',np.object), ('ext',np.int), ('rows', '(1,)i8') ]
        else:
            dtype = [ ('file',np.object), ('ext',np.int), ('rows', (np.int,BalrogSetup.ngal)) ]
        if IsNone(used):
            used = np.array( [], dtype=dtype )

        for catalog in catalogs:
            file = catalog[0][0]
            ext = catalog[0][1]
            hdus = pyfits.open(file)
            data = hdus[ext].data
            
            if (file not in used['file']) or (ext not in used['ext']):
                size = len(data)
                randints = np.random.randint(0,high=size, size=self.ngal)
                selected = data[randints]
                newused = np.array( [(file,ext,randints)], dtype=dtype )
                used = np.concatenate( (used,newused), axis=0 )
            else:
                cut = (used['file']==file) & (used['ext']==ext)
                randints = used[cut]['rows'][0]
                selected = data[randints]

            for tup in catalog[1]:
                index, key, column = tup
                if index > -1:
                    self.component[index][key] = selected[column]
                else:
                    self.galaxy[key] = selected[column]
        return used

    
    def ReturnComponent(self, mkey, mindex=-1):
        if mindex!=-1:
            return self.component[mindex][mkey]
        else:
            return self.galaxy[mkey]


    def DoComponent(self, component, comp=[]):
        completed = copy.copy(comp)
        for i in range(len(component)):
            if i in completed:
                continue

            index, key, mindex, mkey = component[i]
            c = self.ReturnComponent(mkey, mindex=mindex)
            if c is None:
                continue

            if index > -1:
                self.component[index][key] = c
            else:
                self.galaxy[key] = c
            completed.append(i)
        return completed


    def TryRule(self, arg, notready, used):
        if arg.type=='value':
            a = self.ReturnValue(arg.param[0])
        elif arg.type=='array':
            a = np.array( arg.param[0] )
        elif arg.type=='catalog':
            a = self.FunctionCatalog(used, arg.param)

        elif arg.type=='component':
            if type(arg.param)==str:
                arg.param = (-1,arg.param)
            a = self.ReturnComponent(arg.param[1],mindex=arg.param[0])
            if IsNone(a):
                notready = True

        elif arg.type=='function':
            aa, aaa, notready = self.OneFunction(arg.param[0], arg.param[1], arg.param[2], used)
            if notready==False:
                a = arg.param[0](*aa, **aaa)
            else:
                a = None

        return a, notready


    def TryArg(self, arg, notready, used):
        if type(arg).__name__=='CompResult':
            if arg.nProfiles==1:
                arg = arg[0]

        if type(arg).__name__=='Rule':
            a, notready = self.TryRule(arg, notready, used)

        else:
            try:
                length = len(arg)
                ok = True
            except:
                ok = False
            
            a = arg
            if ok:
                if length > 0:
                    if type(arg)==str:
                        return a, notready

                    if type(arg)==tuple:
                        aa = [None]*length
                    else:
                        aa = copy.copy(arg)

                    for i in range(length):
                        aa[i], notready = self.TryArg(arg[i], notready, used)
                        if notready:
                            break
                    a = aa 

        return a, notready


    def OneFunction(self, func, args, kwargs, used):
        notready = False

        arguments = []
        for arg in args:
            if notready:
                break
            a, notready = self.TryArg(arg, notready, used)
            arguments.append(a)

        kwarguments = {}
        for key in kwargs.keys():
            arg = kwargs[key]
            if notready:
                break
            a, notready = self.TryArg(arg, notready, used)
            kwarguments[key] = a

        return [arguments, kwarguments, notready]


    def DoFunction(self, function, used, comp=[]):
        completed = copy.copy(comp)
        for i in range(len(function)):
            if i in completed:
                continue

            index = function[i][0]
            key = function[i][1]
            func = function[i][2]
            args = function[i][3]
            kwargs = function[i][4]
            arguments, kwarguments, notready = self.OneFunction(func, args, kwargs, used)
            
            if notready:
                continue

            if index > -1:
                self.component[index][key] = func(*arguments, **kwarguments)
                try:
                    size = len(self.component[index][key])
                    arr = np.array(self.component[index][key])
                    if arr.ndim!=1:
                        raise FunctionReturnError(501, func.__name__)
                except:
                    raise FunctionReturnError(501, func.__name__)
            else:
                self.galaxy[key] = func(*arguments, **kwarguments)
                try:
                    size = len(self.galaxy[key])
                    arr = np.array(self.galaxy[key])
                    if arr.ndim!=1:
                        raise FunctionReturnError(501, func.__name__)
                except:
                    raise FunctionReturnError(501, func.__name__)
            
            if size!=self.ngal:
                raise FunctionReturnError(501, func.__name__)

            completed.append(i)

        return completed


    def FunctionCatalog(self, used, params):
        get_cat = params[0]
        get_ext = params[1]
        get_col = params[2]
        get_rows = None

        for u in used:
            used_cat, used_ext, used_rows = u
            if used_cat==get_cat and used_ext==get_ext:
                get_cat = used_cat
                get_ext = used_ext
                get_rows = used_rows
                break
       
        cat = pyfits.open(get_cat)[get_ext].data
        size = len(cat)
        if get_rows is None:
            get_rows = np.random.randint(0,high=size, size=self.ngal)
            used.append((get_cat, get_ext,get_rows))
        return cat[get_col][get_rows]

        

    def SortCatalog(self, catalog):
        tables = []
        for i in range(len(catalog)):
            ind, key, file, ext, col = catalog[i]
            table = (file, ext)
           
            index = None
            for j in range(len(tables)):
                if table == tables[j][0]:
                    index = j
            
            if index is None:
                tables.append( [ (file,ext), [(ind,key,col)] ] )
            else:
                tables[index][1].append( (ind,key,col) )
        return tables


    def ChoicesSample(self,rtype, param, i, key, value,array,catalog,component,function):
        if rtype=='value':
            value.append( (i, key, param) )
        elif rtype=='array':
            array.append( (i, key, param) )
        elif rtype=='catalog':
            catalog.append( (i, key, param[0], param[1], param[2]) )
        elif rtype=='component':
            component.append( (i, key, param[0], param[1]) )
        elif rtype=='function':
            function.append( (i, key, param[0], param[1], param[2]) )


    def TryFunctionComponent(self, function, component, used):
        all_f = np.arange(len(function)) 
        last_completed_f = None
        completed_f = []

        all_c = np.arange(len(component))
        last_completed_c = None
        completed_c = []

        while (last_completed_f!=completed_f) or (last_completed_c!=completed_c):
            last_completed_f = completed_f
            last_completed_c = completed_c
            
            completed_f = self.DoFunction(function, used, last_completed_f)
            completed_c = self.DoComponent(component, last_completed_c)


    def GetCompDefault(self, key, BalrogSetup, used, i):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'cosmos.fits')

        if key == 'beta':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of 0.' %(i, key))
            return np.zeros(self.ngal)
        if key=='axisratio':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of 1.' %(i, key))
            return np.ones(self.ngal)
        if key=='halflightradius':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of sampling from the supplied COSMOS catalog.' %(i, key))
            return self.FunctionCatalog(used, [file,1,'HALF_LIGHT_RADIUS'])
        if key=='sersicindex':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of sampling from the supplied COSMOS catalog.' %(i, key))
            return self.FunctionCatalog(used, [file,1,'SERSIC_INDEX'])
        if key=='magnitude' or key=='flux':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of sampling from the supplied COSMOS catalog.' %(i, key))
            return self.FunctionCatalog(used, [file,1,'MAPP_I_SUBARU'])

    def GetGalaxyDefault(self, key, used, BalrogSetup):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'cosmos.fits')

        if key in ['g1', 'g2']:
            BalrogSetup.runlogger.warning('A user-defined rule was not found for %s. Balrog will use the default of 0.' %(key))
            return np.zeros(self.ngal)
        if key == 'magnification':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for %s. Balrog will use the default of 1.' %(key))
            return np.ones(self.ngal)
        if key=='x':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for %s. Balrog will use the default of random positions.' %(key))
            return np.random.uniform( BalrogSetup.xmin, BalrogSetup.xmax, self.ngal )
        if key=='y':
            BalrogSetup.runlogger.warning('A user-defined rule was not found for %s. Balrog will use the default of random positions.' %(key))
            return np.random.uniform( BalrogSetup.ymin, BalrogSetup.ymax, self.ngal )
      

    def SimpleSample(self, BalrogSetup, used=None):
        value = []
        array = []
        catalog = []
        component = []
        function = []
        for i in range(len(self.componentrule)):
            for key in self.componentrule[i].keys():
                rtype = self.componentrule[i][key].type
                param = self.componentrule[i][key].param
                self.ChoicesSample(rtype, param, i, key, value,array,catalog,component,function)
        for key in self.galaxyrule.keys():
            rtype = self.galaxyrule[key].type
            param = self.galaxyrule[key].param
            self.ChoicesSample(rtype, param, -1, key, value,array,catalog,component,function)
        cat = self.SortCatalog(catalog)
       
        np.random.seed(BalrogSetup.seed)
        self.DoValue(value)
        self.DoArray(array)
        used = self.DoCatalog(cat, BalrogSetup, used=used)
        self.TryFunctionComponent(function, component, used)
        return used


    def Sample(self, BalrogSetup):
        used = self.SimpleSample(BalrogSetup)
        for i in range(len(self.component)):
            for key in self.component[i].keys():
                if IsNone(self.component[i][key]):
                    self.component[i][key] = self.GetCompDefault(key, BalrogSetup, used, i)
                    
        for key in self.galaxy.keys():
            if IsNone(self.galaxy[key]):
                self.galaxy[key] = self.GetGalaxyDefault(key, used, BalrogSetup)
        
        np.seterr(over='ignore')
        for i in range(len(self.component)):
            
            bad = ( (BalrogSetup.zeropoint - self.component[i]['flux']) > 50 )
            self.component[i]['flux'][bad] = 0
            self.component[i]['flux'][-bad] = np.power(10.0, (BalrogSetup.zeropoint - self.component[i]['flux'][-bad]) / 2.5)
            if np.sum(bad) > 0:
                BalrogSetup.runlogger.warning('Magnitude value caused overflow when computing flux. Its flux was set to 0. This may not be a problem if a large negative value for magnitude means no detection, e.g. -100')

            '''
            self.component[i]['flux'] = np.power(10.0, (BalrogSetup.zeropoint - self.component[i]['flux']) / 2.5)
            zeros = (self.component[i]['flux']==np.inf)
            if np.sum(zeros) > 0:
                self.component[i]['flux'][zeros] = 0
                BalrogSetup.runlogger.warning('Magnitude value caused overflow when computing flux. Its flux was set to 0. This may not be a problem if a large negative value for magnitude means no detection, e.g. -100')
            '''

            self.component[i]['halflightradius'] = self.component[i]['halflightradius'] * np.sqrt(self.component[i]['axisratio'])
        np.seterr(over='warn')
        return used


    '''
    def GetPSizes(self, BalrogSetup, wcs):
        psizes = np.zeros(self.ngal)
        athresh = np.zeros(self.ngal)
        for i in range(self.ngal):
            psizes[i], athresh[i] = self.GetPSize(i, BalrogSetup, wcs)
        return psizes, athresh


    def GetPSize(self, i, BalrogSetup, wcs):
        fid_seeing = 1.5

        ncomp = len(self.component)
        flux_thresh = BalrogSetup.fluxthresh / float(ncomp)
        test_size = np.zeros(ncomp)
        test_flux = np.zeros(ncomp)
        total_flux = 0
        for j in range(ncomp):
            n = float(self.component[j]['sersicindex'][i])
            reff = float(self.component[j]['halflightradius'][i])
            flux = float(self.component[j]['flux'][i])
            q = float(self.component[j]['axisratio'][i])
            g = np.sqrt( self.galaxy['g1'][i]*self.galaxy['g1'][i] + self.galaxy['g2'][i]*self.galaxy['g2'][i] )
            gq = (1-g) / (1+g)
            k =  self.galaxy['magnification'][i]

            sersicObj = galsim.Sersic(n=n, half_light_radius=reff, flux=flux)
            re = reff / np.sqrt(q)
            re = re / np.sqrt(gq)
            fe = sersicObj.xValue(galsim.PositionD(reff,0))
            f0 = sersicObj.xValue(galsim.PositionD(0,0))
            b = b_n_estimate(n)

            re_frac =  np.power( 1.0 - (1.0 / b) * np.log(flux_thresh/fe), n )
            intrinsic = k * re * re_frac
            seeing = fid_seeing * np.sqrt(2 * np.log(f0/flux_thresh) )
            total = intrinsic + seeing

            test_size[j] = 2 * total
            total_flux += flux
            test_flux[j] = f0

        x = self.galaxy['x'][i]
        y = self.galaxy['y'][i]
        pos = galsim.PositionD(x,y)
        local = wcs.local(image_pos=pos)
        step = min([local.dudx, local.dvdy])
       
        f = max(test_flux)
        psize = max(test_size)
        psize = np.ceil( psize / step )
        if psize%2==0:
            psize += 1

        alias_thresh = BalrogSetup.fluxthresh / total_flux
        return psize, alias_thresh
    '''


    def GetConvolved(self, psfmodel, i, wcs, gsparams, BalrogSetup):
        for j in range(len(self.component)):
            reff = float(self.component[j]['halflightradius'][i])
            if reff < 0:
                reff = 0

            n = float(self.component[j]['sersicindex'][i])
            flux = float(self.component[j]['flux'][i])
            q = float(self.component[j]['axisratio'][i])
            beta = self.component[j]['beta'][i]*galsim.degrees 
            intrinsic_shear = galsim.Shear(q=q, beta=beta )

            sersicObj = galsim.Sersic(n=n, half_light_radius=reff, flux=flux, gsparams=gsparams)
            sersicObj.applyShear(intrinsic_shear)

            if j==0:
                combinedObj = sersicObj
            else:
                combinedObj = combinedObj + sersicObj
       
        lensing_shear = galsim.Shear(g1=self.galaxy['g1'][i], g2=self.galaxy['g2'][i])
        combinedObj.applyShear(lensing_shear)
        combinedObj.applyMagnification(self.galaxy['magnification'][i])

        x = float(self.galaxy['x'][i])
        y = float(self.galaxy['y'][i])
        ix = int(x) 
        iy = int(y)
        dx = x-ix
        dy = y-iy

        pos = galsim.PositionD(x,y)
        local = wcs.local(image_pos=pos)
        localscale = np.sqrt(local.dudx * local.dvdy)
        #combinedObj.applyShift(dx*local.dudx, dy*local.dvdy)
        combinedObj.applyShift(dx*localscale, dy*localscale)

        psf = psfmodel.getPSF(pos, gsparams=gsparams)
        psf.setFlux(1.)
        psf_centroid = psf.centroid()
        psf.applyShift(-psf_centroid.x, -psf_centroid.y)
        #psf.applyShift(psf_centroid.x, psf_centroid.y)

        #pix = galsim.Box(width=local.dudx, height=local.dvdy, gsparams=gsparams)
        #pix = galsim.Pixel(scale=localscale, gsparams=gsparams)
        #combinedObj = galsim.Convolve([psf,pix,combinedObj])
        combinedObj = galsim.Convolve([psf,combinedObj])

        '''
        f1 = os.path.join(BalrogSetup.outdir,'balrog_image','obj%i.fits'%i)
        i1 = combinedObj.draw(scale=localscale)
        galsim.fits.write(image=i1, file_name=f1)

        f2 = os.path.join(BalrogSetup.outdir,'balrog_image','psf%i.fits'%i)
        c = galsim.Convolve([psf,pix])
        c.applyShift(dx*localscale, dy*localscale)
        i2 = c.draw(scale=localscale)
        galsim.fits.write(image=i2, file_name=f2)
        '''
        
        return combinedObj


'''
def b_n_estimate(n):
    order0 = 2*n - 1.0/3.0
    order1 = 4.0 / (405.0 * n)
    order2 = 46.0 / (25515.0 * n*n)
    order3 = 131.0 / (1148175.0 * n*n*n)
    order4 = -2194697.0 / (30690717750 * n*n*n*n)
    sum = order0  + order1 + order2 + order3 + order4
    return sum
'''


class Rule(object):

    def __init__(self, type=None, average=None, sigma=None, joint=False, value=None, array=None, component=None, minimum=None, maximum=None, function=None, args=None, catalog=None, ext=None, column=None, kwargs=None ):

        if type=='catalog':
            if catalog is None:
                raise CatalogArgError(503, 'a catalog file (file)')
            if ext is None:
                raise CatalogArgError(503, 'a FITS extention index (ext)')
            if column is None:
                raise CatalogArgError(503, 'a column name (col)')
            
            try:
                hdus = pyfits.open(catalog)
            except:
                raise CatalogFileError(504, catalog)

            try:
                data = hdus[ext].data
            except:
                raise CatalogExtError(505, catalog, ext)

            try:
                col = data[column]
            except:
                raise CatalogColError(506, catalog, ext, column)

            self.param = [catalog,ext,column]

        elif type=='value':
            if value is None:
                raise Exception('must specify a value with sample type value')
            self.param = value

        elif type=='array':
            if IsNone(array):
                raise Exception('must specify an array with sample type array')
            self.param = array

        elif type=='component':
            if component is None:
                raise Exception('Same takes argument(s)')
            self.param = component

        elif type=='function':
            if function is None:
                raise FunctionArgError(502, 'a function name (function)')
            if args is None:
                raise FunctionArgError(502, 'the arguments to the function (args)')
            self.param = [function, args, kwargs]

        elif type is None:
            self.param = None
            type = None

        else:
            raise Exception('unknown smpling type')

        self.type = type

def IsNone(rule):
    try:
        len(rule)
    except:
        if rule is None:
            return True
    return False


def Value( val=None ):
    return Rule(type='value', value=val)

def Array( arr=None ):
    return Rule(type='array', array=arr)

def Catalog( file=None, ext=None, col=None ):
    return Rule(type='catalog', catalog=file, ext=ext, column=col)

def Column( file=None, ext=None, col=None ):
    return Rule(type='catalog', catalog=file, ext=ext, column=col)

def Same( comp ):
    return Rule(type='component', component=comp)

def Function(function=None, args=(), kwargs={}):
    return Rule(type='function', function=function, args=args, kwargs=kwargs)

class Table(object):
    def __init__(self, file=None, ext=None):
        super(Table, self).__setattr__('file', file)
        super(Table, self).__setattr__('ext',  ext)

    def Column(self, col=None):
        return Catalog(file=self.file, ext=self.ext, col=col)

    def __setattr__(self, name, value):
        raise TableAssignmentError(801)


def Tuplify(g, k):
    if type(g.param)==int:
        g.param = (g.param,k) 
    elif type(g.param)==str:
        g.param = (-1,g.param)
    return g


def MagFlux(g):
    if g.param[1]=='magnitude':
        g.param = (g.param[0], 'flux')
    return g



def HandleFunction(g, k):
    arguments = g.param[1]
    for arg in arguments:
        if type(arg).__name__=='Rule' and arg.type=='component':
            arg = Tuplify(arg,k)
        if type(arg).__name__=='Rule' and arg.type=='function':
            HandleFunction(arg, k)
    g = MagFlux(g)
    return g


def DefineRules(ngal, galkeys, galrules, compkeys, comprules, nProfiles):
    simulatedgals = nComponentSersic(ngal=ngal, ncomp=nProfiles, galkeys=galkeys, compkeys=compkeys)

    keys = galkeys
    for g,k in zip(galrules,keys):
        if g is not None:
            if g.type=='component':
                g = Tuplify(g,k)
                g = MagFlux(g)
            if g.type=='function':
                g = HandleFunction(g,k)
            simulatedgals.GalaxyRule(key=k, rule=g)

    keys = compkeys
    for j in range(len(comprules)):
        key = keys[j]
        if comprules[j] is not None:
            size = len(comprules[j])
            if size!=nProfiles:
                if key=='flux':
                    k = 'magnitude'
                else:
                    k = key
                raise Exception('rules.%s has %i array element(s). Must match rules.nProfiles = %i.' %(k,size,nProfiles))
        for i in range(nProfiles):
            if comprules[j] is not None:
                comp = comprules[j][i]
            else:
                comp = None

            if comp is not None:
                if comp.type=='component':
                    comp = Tuplify(comp, key)
                    comp = MagFlux(comp)
                if comp.type=='function':
                    arg = HandleFunction(comp, key)
                simulatedgals.ComponentRule(component=i, key=key, rule=comp)
                      
    return simulatedgals


def InitializeSersic(rules, sampled, TruthCat, nProfiles=1):
    rules.InitializeSersic(nProfiles=nProfiles)
    sampled.InitializeSersic(nProfiles=nProfiles)
    TruthCat.InitializeSersic(nProfiles=nProfiles) 
