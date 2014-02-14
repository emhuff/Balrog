#!/usr/bin/env python

import numpy as np
import pyfits
import galsim


"""
deprecated stuff 
"""


class BaseSampler():

    """
    def ReadCatalogFrom(self,cat):
        self.catalog = cat


    def PrintRules(self):
        for i in range(len(self.rule)):
            print '%i' %i
           
            for key in self.rule[i].keys():
                type = self.rule[i][key].type
                param = self.rule[i][key].param
                print '%s:  %s --> %s' %(key,type,param)
            
            print ''


    def SampleComponent(self, opts):
        value = []
        gaussian = []
        uniform = []
        array = []
        component = []

        single = []
        joint = []
        for i in range(len(self.rule)):
            for key in self.rule[i].keys():
                rtype = self.rule[i][key].type
                param = self.rule[i][key].param
                if rtype=='gaussian':
                    gaussian.append( (i, key, param[0], param[1]) )
                elif rtype=='uniform':
                    uniform.append( (i, key, param[0], param[1]) )
                elif rtype=='value':
                    value.append( (i, key, param) )
                elif rtype=='array':
                    array.append( (i, key, param) )
                elif rtype=='component':
                    component.append( (i, key, param) )
                elif rtype=='catalog':
                    if param[1]:
                        joint.append( (i, key, param[0]) )
                    else:
                        single.append( (i, key, param[0]) )
        
        for i in range(len(value)):
            index, key, param = self.index_key(value[i])
            val = param[0]
            self.component[index][key][:] = val

        for i in range(len(gaussian)):
            index, key, param = self.index_key(gaussian[i])
            avg, std = param
            self.component[index][key][:] = np.random.normal( avg,std, len(self.component[index][key]) )

        for i in range(len(uniform)):
            index, key, param = self.index_key(uniform[i])
            min, max = param
            self.component[index][key][:] = np.random.uniform( min,max, len(self.component[index][key]) )

        for i in range(len(array)):
            index, key, param = self.index_key(array[i])
            arr = param[0]
            self.component[index][key][:] = arr[:]

        if len(single)>0 or len(joint)>0:
            self.FromCatalog(single,joint, opts.ext)

        for i in range(len(component)):
            index, key, param = self.index_key(component[i])
            comp = param[0]
            if type(comp)==int:
                self.component[index][key] = self.component[comp][key]
            elif type(comp)==str:
                exec 'self.component[%i]["%s"] = self.%s["%s"]' %(index,key, comp,key)


    def SampleGalaxy(self, opts, redo=None):
        # redo = ('x', [1,8,11])
        value = []
        gaussian = []
        uniform = []
        array = []

        single = []
        joint = []

        if redo==none:
            okkeys = self.galaxyrule.keys()
            okindecies = np.arange(self.ngal)
        else:
            okkeys = redo[0]
            okindecies = redo[1]

        #for key in self.galaxyrule.keys():
        for key in okkeys:
            rtype = self.galaxyrule[key].type
            param = self.galaxyrule[key].param
            if rtype=='gaussian':
                gaussian.append( (key, param[0], param[1]) )
            elif rtype=='uniform':
                uniform.append( (key, param[0], param[1]) )
            elif rtype=='value':
                value.append( (key, param) )
            elif rtype=='array':
                array.append( (key, param) )
            elif rtype=='catalog':
                if param[1]:
                    joint.append( (key, param[0]) )
                else:
                    single.append( (key, param[0]) )
        
        for i in range(len(value)):
            key, param = self.noindex_key(value[i])
            val = param[0]
            self.galaxy[key][okindecies] = val

        for i in range(len(gaussian)):
            key, param = self.noindex_key(gaussian[i])
            avg, std = param
            self.galaxy[key][okindecies] = np.random.normal( avg,std, len(okindecies) )

        for i in range(len(uniform)):
            key, param = self.noindex_key(uniform[i])
            min, max = param
            self.galaxy[key][okindecies] = np.random.uniform( min,max, len(okindecies) )

        for i in range(len(array)):
            key, param = self.noindex_key(array[i])
            arr = param[0]
            self.galaxy[key][okindecies] = arr[okindecies]

        if len(single)>0 or len(joint)>0:
            self.FromCatalogGalaxy(single,joint, okindecies, opts.ext)


    def OutsidePostageStamp(self, psizes, opts, index=None):
        if index==None:
            index = np.arange(len(psizes))

        inds = []
        for ii in range(len(psizes)):
            bad = False
            i = index[ii] 
            psize = psizes[ii]
            x = int(self.galaxy['x'][i])
            y = int(self.galaxy['y'][i])
        
            if (x - (psize/2+opts.frame)) < opts.xmin:
                bad = True
            if (x + (psize/2+opts.frame)) > opts.xmax:
                bad = True
            if (y - (psize/2+opts.frame)) < opts.ymin:
                bad = True
            if (y + (psize/2+opts.frame)) > opts.ymax:
                bad = True

            if bad:
                inds.append(i)

        return np.array(inds)


    def FromCatalog(self, single, joint, ext):
        hdus = pyfits.open(self.catalog)
        data = hdus[ext].data

        #sanitycuts = (data['FLUX_AUTO'] > 0) & (data['FLUX_RADIUS'] > 0) & (data['FWHM_IMAGE'] > 0) & (data['FLAGS']==0)
        #data = data[sanitycuts]

        size = len(data)
        
        for i in range(len(single)):
            index, key, param = self.index_key(single[i])
            randints = np.random.randint(0,high=size, size=self.ngal)
            column = param[0]
            self.component[index][key] = data[randints][column]

        if len(joint)>0:
            randints = np.random.randint(0,high=size, size=self.ngal)
            selected = data[randints]
            for i in range(len(joint)):
                index, key, param = self.index_key(joint[i])
                column = param[0]
                self.component[index][key] = selected[column]


    def FromCatalogGalaxy(self, single, joint, okindecies, ext):
        hdus = pyfits.open(self.catalog)
        data = hdus[ext].data
        size = len(data)
        
        for i in range(len(single)):
            key, param = self.noindex_key(single[i])
            randints = np.random.randint(0,high=size, size=len(okindecies))
            column = param[0]
            self.galaxy[key][okindecies] = data[randints][column]

        if len(joint)>0:
            randints = np.random.randint(0,high=size, size=len(okindecies))
            selected = data[randints]
            for i in range(len(joint)):
                key, param = self.noindex_key(joint[i])
                column = param[0]
                self.galaxy[key][okindecies] = selected[column]


    def index_key(self, arr):
        return [arr[0],arr[1],arr[2:]]

    def noindex_key(self, arr):
        return [arr[0],arr[1:]]

    """

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

    

    def DoGaussian(self, gaussian):
        np.random.seed()
        for i in range(len(gaussian)):
            index, key, avg, std = gaussian[i]
            if index > -1:
                self.component[index][key][:] = np.random.normal( avg,std, len(self.component[index][key]) )
            else:
                self.galaxy[key][:] = np.random.normal( avg,std, len(self.galaxy[key]) )


    def DoUniform(self, uniform):
        np.random.seed()
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


    def DoCatalog(self, catalogs):
        np.random.seed()
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


    def BetterSample(self):
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
       
        self.DoGaussian(gaussian)
        self.DoUniform(uniform)
        self.DoValue(value)
        self.DoArray(array)
        self.DoCatalog(cat)
        self.DoComponent(component)

    
    def Sample(self, opts, psfmodel):
        """
        self.SampleComponent(opts)
        self.SampleGalaxy(opts)
        """
        self.BetterSample()

        for i in range(len(self.component)):
            self.component[i]['flux'] = np.power(10.0, (opts.zeropoint - self.component[i]['flux']) / 2.5)
            self.component[i]['halflightradius'] = self.component[i]['halflightradius'] * np.sqrt(self.component[i]['axisratio'])
        psizes = self.GetPostageStampSizes(opts, psfmodel, stamp=True)
        
        """
        outside = self.OutsidePostageStamp(psizes,opts)
        while len(outside) > 0:
            self.SampleGalaxy(opts, redo=('x',outside))
            self.SampleGalaxy(opts, redo=('y',outside))
            ps = self.GetPostageStampSizes(opts, psfmodel, index=outside, stamp=True )
            psizes[outside] = ps
            outside = self.OutsidePostageStamp(ps,opts, index=outside
        """
           
        '''
        for i in range(len(self.component)):
            self.component[i]['flux'] = np.power(10.0, (opts.zeropoint - self.component[i]['flux']) / 2.5)
        '''

        return psizes




    def GetPSFConvolved(self, opts, psfmodel, i, stamp=False):
        for j in range(len(self.component)):
            n = float(self.component[j]['sersicindex'][i])
            reff = float(self.component[j]['halflightradius'][i])
            flux = float(self.component[j]['flux'][i])
            q = float(self.component[j]['axisratio'][i])
            if stamp==True:
                beta = 0 * galsim.degrees
            else:
                beta = self.component[j]['beta'][i]*galsim.degrees 
            intrinsic_shear = galsim.Shear(q=q, beta=beta )
            
            sersicObj = galsim.Sersic(n=n, half_light_radius=reff, flux=flux)
            sersicObj.applyShear(intrinsic_shear)
            if j==0:
                combinedObj = sersicObj
            else:
                combinedObj = combinedObj + sersicObj
        
        if stamp==True:
            gmag = np.sqrt( self.galaxy['g1'][i]*self.galaxy['g1'][i] + self.galaxy['g2'][i]*self.galaxy['g2'][i] )
            lensing_shear = galsim.Shear(g1=float(gmag), g2=0)
        else:
            lensing_shear = galsim.Shear(g1=self.galaxy['g1'][i], g2=self.galaxy['g2'][i])
        combinedObj.applyShear(lensing_shear)
        combinedObj.applyMagnification(self.galaxy['magnification'][i])

        x = float(self.galaxy['x'][i])
        y = float(self.galaxy['y'][i])
        ix = int(np.floor(x)) 
        iy = int(np.floor(y))
        dx = x-ix
        dy = y-iy
        pos = galsim.PositionD(x,y)
        combinedObj.applyShift(dx*opts.pixscale,dy*opts.pixscale)
        
        psf = psfmodel.getPSF(pos,opts.pixscale)
        psf.setFlux(1.)
        combinedObj = galsim.Convolve([psf,combinedObj])

        return combinedObj


    def GetPostageStampSizes(self,opts, psfmodel, index=None, stamp=True):
        if index==None:
            index = np.arange( self.ngal )
        psizes = []

        for i in index:
            combinedObj = self.GetPSFConvolved(opts, psfmodel, i, stamp=stamp)
            size = opts.minsize
            while (combinedObj.xValue(galsim.PositionD(-size*opts.pixscale/2,0))>opts.fluxthresh) or (combinedObj.xValue(galsim.PositionD(size*opts.pixscale/2,0))>opts.fluxthresh):
                size += opts.inc
            if size%2==0:
                size += 1
            psizes.append(size)

        return np.array(psizes)



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

    return simulatedgals
