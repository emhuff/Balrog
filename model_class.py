#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import galsim


"""
deprecated stuff 
"""


class BaseSampler():


    def __init__(self, ngal, ncomp, label, catalog):
        self.init(ngal,ncomp,label,catalog)


    def initrules(self):
        self.ComponentRule(0, 'beta', Rule(type='uniform', minimum=-90, maximum=90) )
        self.ComponentRule(0, 'axisratio', Rule(type='uniform', minimum=0, maximum=1) )
        for i in range(len(self.component)-1):
            self.ComponentRule(i+1, 'beta', Rule(type='component', component=0) )
            self.ComponentRule(i+1, 'axisratio', Rule(type='component', component=0) )
        for i in range(len(self.component)):
            self.ComponentRule(i, 'flux', Rule(type='uniform', minimum=25, maximum=30) )
            self.ComponentRule(i, 'halflightradius', Rule(type='uniform', minimum=0.3, maximum=1.0) )
            self.ComponentRule(i, 'sersicindex', Rule(type='value', value=i+1) )

        self.GalaxyRule('g1', Rule(type='value', value=0) )
        self.GalaxyRule('g2', Rule(type='value', value=0) )
        self.GalaxyRule('magnification', Rule(type='value', value=1) )


    def init(self, ngal, ncomp, label, catalog):
        self.ngal = ngal
        self.component = [None]*ncomp
        self.rule = [None]*ncomp
        self.catalog = catalog

        self.galaxy = self._DefaultGalaxy(ngal)
        self.galaxyrule = self._DefaultRule(self.galaxy)

        for i in range(ncomp):
            self.component[i] = self._DefaultSersicComponent(ngal)
            self.rule[i] = self._DefaultRule(self.component[i])

        if label!=None:
            if len(label)!=len(self.component):
                raise Exception('labels for for components must have same length = ncomp')
            self.ruledict = {}
            for i in range(len(label)):
                exec 'self.%s = self.component[%i]' %(label[i],i)
                exec 'self.ruledict["%s"] = self.rule[%i]' %(label[i],i)

        self.initrules()


    def _DefaultSersicComponent(self,ngal):
        dict = {}
        dict['flux'] = np.zeros(ngal)
        dict['halflightradius'] = np.zeros(ngal)
        dict['sersicindex'] = np.zeros(ngal)
        dict['axisratio'] = np.zeros(ngal)
        dict['beta'] = np.zeros(ngal)
        return dict

    def _DefaultGalaxy(self,ngal):
        dict = {}
        dict['x'] = np.zeros(ngal)
        dict['y'] = np.zeros(ngal)
        dict['g1'] = np.zeros(ngal)
        dict['g2'] = np.zeros(ngal)
        dict['magnification'] = np.zeros(ngal)
        return dict


    def _DefaultRule(self,comp):
        dict = {}
        for key in comp.keys():
            dict[key] = Rule()
        return dict

    
    def ComponentRules(self, component, rules):
        if type(component)==str:
            ruleset = self.ruledict[component]
        elif type(component)==int:
            ruleset = self.rule[component]

        for key in rules:
            ruleset[key] = rules[key]


    def ComponentRule(self, component=0, key=None, rule=None):
        
        if key==None:
            raise Exception('must indicate which field to give a rule for')
        if key=='mag':
            key = 'flux'

        if rule==None:
            raise Exception('must give a rule')

        if type(component)==str:
            ruleset = self.ruledict[component]
        elif type(component)==int:
            ruleset = self.rule[component]

        ruleset[key] = rule

    def GalaxyRule(self, key=None, rule=None):

        if key==None:
            raise Exception('must indicate which field to give a rule for')
        if rule==None:
            raise Exception('must give a rule')
        """
        if rule.type=='component':
            raise Exception('rule type component is not valid with galaxy')
        """

        self.galaxyrule[key] = rule

    

    def DoGaussian(self, gaussian, seed):
        np.random.seed(seed)
        for i in range(len(gaussian)):
            index, key, avg, std = gaussian[i]
            if index > -1:
                self.component[index][key][:] = np.random.normal( avg,std, len(self.component[index][key]) )
            else:
                self.galaxy[key][:] = np.random.normal( avg,std, len(self.galaxy[key]) )


    def DoUniform(self, uniform, seed):
        np.random.seed(seed)
        for i in range(len(uniform)):
            index, key, min, max = uniform[i]
            if index > -1:
                self.component[index][key][:] = np.random.uniform( min,max, len(self.component[index][key]) )
            else:
                self.galaxy[key][:] = np.random.uniform( min,max, len(self.galaxy[key]) )


    def DoValue(self, value):
        for i in range(len(value)):
            index, key, val = value[i]
            if index > -1:
                self.component[index][key][:] = val
            else:
                self.galaxy[key][:] = val


    def DoArray(self,array):
        for i in range(len(array)):
            index, key, arr = array[i]
            if index > -1:
                self.component[index][key][:] = arr[:]
            else:
                self.galaxy[key][:] = arr[:]


    def DoCatalog(self, catalogs, seed):
        np.random.seed(seed)
        for catalog in catalogs:
            file = catalog[0][0]
            ext = catalog[0][1]
            hdus = pyfits.open(file)
            data = hdus[ext].data
            size = len(data)
            randints = np.random.randint(0,high=size, size=self.ngal)
            selected = data[randints]

            for tup in catalog[1]:
                index, key, column = tup
                if index > -1:
                    self.component[index][key] = selected[column]
                else:
                    self.galaxy[key] = selected[column]


    def DoComponent(self, component):
        for i in range(len(component)):
            index, key, mindex, mkey = component[i]
            if index > -1:
                if mindex > -1:
                    self.component[index][key] = self.component[mindex][mkey]
                else:
                    self.component[index][key] = self.galaxy[mkey]
            else:
                if mindex > -1:
                    self.galaxy[key] = self.component[mindex][mkey]
                else:
                    self.galaxy[key] = self.galaxy[mkey]


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


    def ChoicesSample(self,rtype, param, i, key, gaussian,uniform,value,array,catalog,component):
        if rtype=='gaussian':
            gaussian.append( (i, key, param[0], param[1]) )
        elif rtype=='uniform':
            uniform.append( (i, key, param[0], param[1]) )
        elif rtype=='value':
            value.append( (i, key, param) )
        elif rtype=='array':
            array.append( (i, key, param) )
        elif rtype=='catalog':
            catalog.append( (i, key, param[0], param[1], param[2]) )
        elif rtype=='component':
            component.append( (i, key, param[0], param[1]) )


    def BetterSample(self, BalrogSetup):
        gaussian = []
        uniform = []
        value = []
        array = []
        catalog = []
        component = []
        for i in range(len(self.rule)):
            for key in self.rule[i].keys():
                rtype = self.rule[i][key].type
                param = self.rule[i][key].param
                self.ChoicesSample(rtype, param, i, key, gaussian,uniform,value,array,catalog,component)
        for key in self.galaxyrule.keys():
            rtype = self.galaxyrule[key].type
            param = self.galaxyrule[key].param
            self.ChoicesSample(rtype, param, -1, key, gaussian,uniform,value,array,catalog,component)
        cat = self.SortCatalog(catalog)
       
        self.DoGaussian(gaussian, BalrogSetup.seed)
        self.DoUniform(uniform, BalrogSetup.seed)
        self.DoValue(value)
        self.DoArray(array)
        self.DoCatalog(cat, BalrogSetup.seed)
        self.DoComponent(component)

    
    def Sample(self, BalrogSetup):
        self.BetterSample(BalrogSetup)

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
        fid_seeing = 1.0

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


class Rule():

    def __init__(self, type=None, average=None, sigma=None, joint=False, value=None, array=None, component=None, minimum=None, maximum=None, function=None, args=None, catalog=None, ext=None, column=None, ):

        if type=='gaussian':
            if average==None:
                raise Exception('must specifiy an average with sample type gaussian')
            if sigma==None:
                raise Exception('must specify a sigma with sample type gaussian')
            self.param = [average,sigma]

        elif type=='uniform':
            if minimum==None:
                raise Exception('must specify a minimum with sample type uniform')
            if maximum==None:
                raise Exception('must sepcify a maximum with sample type uniform')
            self.param = [minimum, maximum]

        elif type=='catalog':
            if column==None:
                raise Exception('must specify a column with sample type catalog')
            #self.param = [column,joint]
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
                raise Exception('must specify a component wth sample type component')
            self.param = component

        elif type==None:
            self.param = 0
            type = 'value'

        else:
            raise Exception('unknown type')

        '''
        elif type=='function':
            if component==None:
                raise Exception('must specify a component wth sample type function')
            if function==None:
                raise Exception('must specify a function wth sample type function')
            if args==None:
                raise Exception('must specify args wth sample type function')
            self.param = [component, function, args]
        '''

        self.type = type



class nComponentSersic(BaseSampler):
    def __init__(self, ngal=100, ncomp=2, label=None, catalog=None):
        self.init(ngal, ncomp, label, catalog)
    

"""
class Disk(BaseSampler):
    def __init__(self, ngal=100, catalog=None):
        self.init(ngal, 1, ['profile'], catalog)


class Bulge(BaseSampler):
    def __init__(self, ngal=100, catalog=None) :
        self.init(ngal, 1, ['profile'], catalog) 


class BulgeDisk(BaseSampler):
    def __init__(self, ngal=100, catalog=None):
        self.init(ngal, 2, ['bulge','disk'], catalog)


class DiskBulge(BaseSampler):
    def __init__(self, ngal=100, catalog=None):
        self.init(ngal, 2, ['disk','bulge'], catalog)


class deVaucouleur(BaseSampler):
    def __init__(self, ngal=100, catalog=None):
        self.init(ngal, 1, ['profile'], catalog) 
        self.ComponentRule(0, 'sersicindex', Rule(type='value',value=4))


class Exponential(BaseSampler):
    def __init__(self, ngal=100, catalog=None): 
        self.init(ngal, 1, ['profile'], catalog) 
        self.ComponentRule(0, 'sersicindex', Rule(type='value',value=1))


class im3shape(BaseSampler):
    def __init__(self, ngal=100, catalog=None):
        self.init(ngal, 2, ['bulge','disk'], catalog)
        self.ComponentRule('bulge', 'sersicindex', Rule(type='value', value=4))
        self.ComponentRule('disk', 'sersicindex', Rule(type='value', value=1))
        self.ComponentRule(1, 'halflightradius', Rule(type='component', component=0))
"""


def Random( min,max ):
    return Rule(type='uniform', minimum=min, maximum=max)

def Value( val ):
    return Rule(type='value', value=val)

def Gaussian( avg,std ):
    return Rule(type='gaussian', average=avg, sigma=std)

def Array( arr ):
    return Rule(type='array', array=arr)

def Catalog( file,ext,col ):
    return Rule(type='catalog', catalog=file, ext=ext, column=col)

def Same( comp ):
    return Rule(type='component', component=comp)


def DefineRules(args, x=None, y=None, g1=None, g2=None, magnification=None, nProfiles=1, axisratio=None, beta=None, halflightradius=None, magnitude=None, sersicindex=None ):
    simulatedgals = nComponentSersic(ngal=args.ngal, catalog=None, ncomp=nProfiles)
    out = open(args.simruleslog, 'w')

    galrules = [x, y, g1, g2, magnification]
    keys = ['x', 'y', 'g1', 'g2', 'magnification']
    for g,k in zip(galrules,keys):
        if g!=None:
            if g.type=='component':
                if type(g.param)==int:
                    g.param = (g.param,k) 
                elif type(g.param)==str:
                    g.param = (-1,k)
            simulatedgals.GalaxyRule(key=k, rule=g)
        out.write('%s %s %s\n' %(k, g.type, str(g.param)) )
   
    out.write('\n')
    keys = ['axisratio', 'beta', 'halflightradius', 'flux', 'sersicindex']
    comprules = [axisratio, beta, halflightradius, magnitude, sersicindex]
    for j in range(len(comprules)):
        key = keys[j]
        size = len(comprules[j])
        if size!=nProfiles:
            raise Exception('%s has %i elements. Must match nProfiles = %i' %(key,size,nProfiles))
        for i in range(nProfiles):
            comp = comprules[j][i]
            if comp!=None:
                if comp.type=='component':
                    if type(comp.param)==int:
                        comp.param = (comp.param,key) 
                    elif type(comp.param)==str:
                        comp.param = (-1,key)
                simulatedgals.ComponentRule(component=i, key=key, rule=comp)
            out.write('%s %s %s %s\n' %(str(i), key, comp.type, str(comp.param)) )

    out.close()
    return simulatedgals
