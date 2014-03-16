#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import galsim
import copy
import sys
import os
import logging
import traceback
from balrogexcept import *


class nComponentSersic(object):

    def __init__(self, ngal=100, ncomp=2):
        self.ngal = ngal
        self.component = [None]*ncomp
        self.componentrule = [None]*ncomp

        self.galaxy = self._InitGalaxy(ngal)
        self.galaxyrule = self._InitRule(self.galaxy)

        for i in range(ncomp):
            self.component[i] = self._InitSersicComponent(ngal)
            self.componentrule[i] = self._InitRule(self.component[i])


    def _InitSersicComponent(self,ngal):
        dict = {}
        dict['flux'] = None
        dict['halflightradius'] = None
        dict['sersicindex'] = None
        dict['axisratio'] = None
        dict['beta'] = None
        return dict

    def _InitGalaxy(self,ngal):
        dict = {}
        dict['x'] = None
        dict['y'] = None
        dict['g1'] = None
        dict['g2'] = None
        dict['magnification'] = None
        return dict


    def _InitRule(self,comp):
        dict = {}
        for key in comp.keys():
            dict[key] = Rule()
        return dict


    def ComponentRule(self, component=0, key=None, rule=None):
        
        if key==None:
            raise Exception('must indicate which field to give a rule for')
        if key=='mag':
            key = 'flux'

        if rule==None:
            raise Exception('must give a rule')

        ruleset = self.componentrule[component]
        ruleset[key] = rule


    def GalaxyRule(self, key=None, rule=None):

        if key==None:
            raise Exception('must indicate which field to give a rule for')
        if rule==None:
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
    

    def DoCatalog(self, catalogs):
        used = []
        for catalog in catalogs:
            file = catalog[0][0]
            ext = catalog[0][1]
            hdus = pyfits.open(file)
            data = hdus[ext].data
            size = len(data)
            randints = np.random.randint(0,high=size, size=self.ngal)
            selected = data[randints]
            used.append( (file, ext, randints) )

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
            index, key, mindex, mkey = component[i]
            c = self.ReturnComponent(mkey, mindex=mindex)
            if c==None:
                continue

            if index > -1:
                self.component[index][key] = c
            else:
                self.galaxy[key] = c
            completed.append(i)
        return completed
       

    def OneFunction(self, func, args, used):
        arguments = []
        notready = False
        for arg in args:
            if type(arg).__name__=='Rule':
                if arg.type=='value':
                    a = self.ReturnValue(arg.param[0])
                elif arg.type=='array':
                    a = np.array( arg.param[0] )
                elif arg.type=='catalog':
                    a = self.FunctionCatalog(used, arg.param)

                elif arg.type=='component':
                    a = self.ReturnComponent(arg.param[1],mindex=arg.param[0])
                    if a==None:
                        notready = True

                elif arg.type=='function':
                    aa, notready = self.OneFunction(arg.param[0], arg.param[1], used)
                    if notready==False:
                        a = arg.param[0](aa)

                if notready:
                    break

            else:
                a = arg

            arguments.append(a)
        return [arguments, notready]


    def DoFunction(self, function, used, comp=[]):
        completed = copy.copy(comp)
        for i in range(len(function)):
            if i in completed:
                continue

            index = function[i][0]
            key = function[i][1]
            func = function[i][2]
            args = function[i][3]
            arguments, notready = self.OneFunction(func, args, used)

            if notready:
                continue

            if index > -1:
                self.component[index][key] = func(*arguments)
            else:
                self.galaxy[key] = func(*arguments)

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
        if get_rows==None:
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
            
            if index==None:
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
            function.append( (i, key, param[0], param[1]) )


    def TryFunctionComponent(self, function, component, used):
        all_f = np.arange(len(function)) 
        last_completed_f = None
        completed_f = []

        all_c = np.arange(len(component))
        last_completed_c = None
        completed_c = []

        while (last_completed_f!=completed_f) and (last_completed_c!=completed_c):
            last_completed_f = completed_f
            last_completed_c = completed_c

            completed_f = self.DoFunction(function, used, last_completed_f)
            completed_c = self.DoComponent(component, last_completed_c)
            

    def GetCompDefault(self, key, BalrogSetup, used, i):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'cosmos.fits')

        if key == 'beta':
            BalrogSetup.logger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of 0.' %(i, key))
            return np.zeros(self.ngal)
        if key=='axisratio':
            BalrogSetup.logger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of 1.' %(i, key))
            return np.ones(self.ngal)
        if key=='halflightradius':
            BalrogSetup.logger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of sampling from the supplied COSMOS catalog.' %(i, key))
            return self.FunctionCatalog(used, [file,1,'HALF_LIGHT_RADIUS'])
        if key=='sersicindex':
            BalrogSetup.logger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of sampling from the supplied COSMOS catalog.' %(i, key))
            return self.FunctionCatalog(used, [file,1,'SERSIC_INDEX'])
        if key=='magnitude' or key=='flux':
            BalrogSetup.logger.warning('A user-defined rule was not found for component %i of %s. Balrog will use the default of sampling from the supplied COSMOS catalog.' %(i, key))
            return self.FunctionCatalog(used, [file,1,'IMAG'])

    def GetGalaxyDefault(self, key, used, BalrogSetup):
        thisdir = os.path.dirname( os.path.realpath(__file__) )
        file = os.path.join(thisdir, 'cosmos.fits')

        if key in ['g1', 'g2']:
            BalrogSetup.logger.warning('A user-defined rule was not found for %s. Balrog will use the default of 0.' %(key))
            return np.zeros(self.ngal)
        if key == 'magnification':
            BalrogSetup.logger.warning('A user-defined rule was not found for %s. Balrog will use the default of 1.' %(key))
            return np.ones(self.ngal)
        if key=='x':
            BalrogSetup.logger.warning('A user-defined rule was not found for %s. Balrog will use the default of random positions.' %(key))
            return np.random.uniform( BalrogSetup.xmin, BalrogSetup.xmax, self.ngal )
        if key=='y':
            BalrogSetup.logger.warning('A user-defined rule was not found for %s. Balrog will use the default of random positions.' %(key))
            return np.random.uniform( BalrogSetup.ymin, BalrogSetup.ymax, self.ngal )
        

    def Sample(self, BalrogSetup):
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
        used = self.DoCatalog(cat)
        self.TryFunctionComponent(function, component, used)

        for i in range(len(self.component)):
            for key in self.component[i].keys():
                if self.component[i][key]==None:
                    self.component[i][key] = self.GetCompDefault(key, BalrogSetup, used, i)
                    
        for key in self.galaxy.keys():
            if self.galaxy[key]==None:
                self.galaxy[key] = self.GetGalaxyDefault(key, used, BalrogSetup)
         
        for i in range(len(self.component)):
            self.component[i]['flux'] = np.power(10.0, (BalrogSetup.zeropoint - self.component[i]['flux']) / 2.5)
            self.component[i]['halflightradius'] = self.component[i]['halflightradius'] * np.sqrt(self.component[i]['axisratio'])


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
          


    def GetPSFConvolved(self, psfmodel, i, wcs, athresh):
        #gsparams = galsim.GSParams(alias_threshold=athresh[i])#, kvalue_accuracy=athresh[i])
        
        for j in range(len(self.component)):
            n = float(self.component[j]['sersicindex'][i])
            reff = float(self.component[j]['halflightradius'][i])
            flux = float(self.component[j]['flux'][i])
            q = float(self.component[j]['axisratio'][i])
            beta = self.component[j]['beta'][i]*galsim.degrees 
            intrinsic_shear = galsim.Shear(q=q, beta=beta )

            #sersicObj = galsim.Sersic(n=n, half_light_radius=reff, flux=flux, gsparams=gsparams)
            sersicObj = galsim.Sersic(n=n, half_light_radius=reff, flux=flux)
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
        combinedObj.applyShift(dx*local.dudx, dy*local.dvdy)

        #psf = psfmodel.getPSF(pos, gsparams=gsparams)
        psf = psfmodel.getPSF(pos)
        psf.setFlux(1.)
        psf_centroid = psf.centroid()
        psf.applyShift(-psf_centroid.x, -psf_centroid.y)
        combinedObj = galsim.Convolve([psf,combinedObj])

        #pix = galsim.Box(width=local.dudx, height=local.dvdy, gsparams=gsparams)
        pix = galsim.Box(width=local.dudx, height=local.dvdy)
        combinedObj = galsim.Convolve([pix,combinedObj])
        
        return combinedObj


def b_n_estimate(n):
    order0 = 2*n - 1.0/3.0
    order1 = 4.0 / (405.0 * n)
    order2 = 46.0 / (25515.0 * n*n)
    order3 = 131.0 / (1148175.0 * n*n*n)
    order4 = -2194697.0 / (30690717750 * n*n*n*n)
    sum = order0  + order1 + order2 + order3 + order4
    return sum


class Rule(object):

    def __init__(self, type=None, average=None, sigma=None, joint=False, value=None, array=None, component=None, minimum=None, maximum=None, function=None, args=None, catalog=None, ext=None, column=None, ):


        if type=='catalog':
            if catalog==None:
                raise Exception('must specify a catalog file when sampling from catalog')
            if ext==None:
                raise Exception('must specify an extenstion when sampling from catalog')
            if column==None:
                raise Exception('must specify a column when sampling from catalog')
            self.param = [catalog,ext,column]

        elif type=='value':
            if value==None:
                raise Exception('must specify a value with sample type value')
            self.param = value

        elif type=='array':
            if array==None:
                raise Exception('must specify an array with sample type array')
            self.param = array

        elif type=='component':
            if component==None:
                raise Exception('Same takes argument(s)')
            self.param = component

        elif type=='function':
            if function==None:
                raise Exception('must specify a function with sample type function')
            if args==None:
                raise Exception('must specify args with sample type function')
            self.param = [function, args]

        elif type==None:
            self.param = None
            type = None

        else:
            raise Exception('unknown smpling type')

        self.type = type


def Value( val=None ):
    return Rule(type='value', value=val)

def Array( arr=None ):
    return Rule(type='array', array=arr)

def Catalog( file=None, ext=None, col=None ):
    return Rule(type='catalog', catalog=file, ext=ext, column=col)

def Same( comp ):
    return Rule(type='component', component=comp)

def Function(function=None, args=()):
    return Rule(type='function', function=function, args=args)


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


def DefineRules(opts, x=None, y=None, g1=None, g2=None, magnification=None, nProfiles=1, axisratio=None, beta=None, halflightradius=None, magnitude=None, sersicindex=None ):
    simulatedgals = nComponentSersic(ngal=opts.ngal, ncomp=nProfiles)
    out = open(opts.simruleslog, 'w')

    galrules = [x, y, g1, g2, magnification]
    keys = ['x', 'y', 'g1', 'g2', 'magnification']
    for g,k in zip(galrules,keys):
        if g!=None:
            if g.type=='component':
                g = Tuplify(g,k)
                g = MagFlux(g)
            if g.type=='function':
                g = HandleFunction(g,k)
            simulatedgals.GalaxyRule(key=k, rule=g)

        if g!=None:
            out.write('%s %s %s\n' %(k, g.type, str(g.param)) )
        else:
            out.write('%s None\n' %(k))
   
    out.write('\n')
    keys = ['axisratio', 'beta', 'halflightradius', 'flux', 'sersicindex']
    comprules = [axisratio, beta, halflightradius, magnitude, sersicindex]
    for j in range(len(comprules)):
        key = keys[j]
        if comprules[j]!=None:
            size = len(comprules[j])
            if size!=nProfiles:
                if key=='flux':
                    k = 'magnitude'
                else:
                    k = key
                raise Exception('rules.%s has %i array element(s). Must match rules.nProfiles = %i.' %(k,size,nProfiles))
        for i in range(nProfiles):
            if comprules[j]!=None:
                comp = comprules[j][i]
            else:
                comp = None

            if comp!=None:
                if comp.type=='component':
                    comp = Tuplify(comp, key)
                    comp = MagFlux(comp)
                if comp.type=='function':
                    arg = HandleFunction(comp, key)
                simulatedgals.ComponentRule(component=i, key=key, rule=comp)
                      
            if comp!=None:
                out.write('%s %s %s %s\n' %(str(i), key, comp.type, str(comp.param)) )
            else:
                out.write('%s %s None\n' %(str(i), key))

    out.close()
    return simulatedgals
