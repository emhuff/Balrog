#!/usr/bin/env python

import numpy as np
import pyfits


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
            self.ComponentRule(i, 'flux', Rule(type='uniform', minimum=1.0e2, maximum=1.0e4) )
            self.ComponentRule(i, 'halflightradius', Rule(type='uniform', minimum=2, maximum=10) )
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

    
    def ReadCatalogFrom(self,cat):
        self.catalog = cat

    
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
        if rule.type=='component':
            raise Exception('rule type component is not valid with galaxy')

        self.galaxyrule[key] = rule

    
    def PrintRules(self):
        for i in range(len(self.rule)):
            print '%i' %i
           
            for key in self.rule[i].keys():
                type = self.rule[i][key].type
                param = self.rule[i][key].param
                print '%s:  %s --> %s' %(key,type,param)
            
            print ''


    def SampleComponent(self):
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
            self.FromCatalog(single,joint)

        for i in range(len(component)):
            index, key, param = self.index_key(component[i])
            comp = param[0]
            if type(comp)==int:
                self.component[index][key] = self.component[comp][key]
            elif type(comp)==str:
                exec 'self.component[%i]["%s"] = self.%s["%s"]' %(index,key, comp,key)


    def SampleGalaxy(self):
        value = []
        gaussian = []
        uniform = []
        array = []

        single = []
        joint = []

        for key in self.galaxyrule.keys():
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
            self.galaxy[key][:] = val

        for i in range(len(gaussian)):
            key, param = self.noindex_key(gaussian[i])
            avg, std = param
            self.galaxy[key][:] = np.random.normal( avg,std, len(self.galaxy[key]) )

        for i in range(len(uniform)):
            key, param = self.noindex_key(uniform[i])
            min, max = param
            self.galaxy[key][:] = np.random.uniform( min,max, len(self.galaxy[key]) )

        for i in range(len(array)):
            key, param = self.noindex_key(array[i])
            arr = param[0]
            self.galaxy[key][:] = arr[:]

        if len(single)>0 or len(joint)>0:
            self.FromCatalogGalaxy(single,joint)


    def Sample(self):
        self.SampleComponent()
        self.SampleGalaxy()


    def FromCatalog(self, single, joint):
        hdus = pyfits.open(self.catalog)
        data = hdus[2].data

        sanitycuts = (data['FLUX_AUTO'] > 0) & (data['FLUX_RADIUS'] > 0) & (data['FWHM_IMAGE'] > 0) & (data['FLAGS']==0)
        data = data[sanitycuts]

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


    def FromCatalogGalaxy(self, single, joint):
        hdus = pyfits.open(self.catalog)
        data = hdus[2].data
        size = len(data)
        
        for i in range(len(single)):
            key, param = self.noindex_key(single[i])
            randints = np.random.randint(0,high=size, size=self.ngal)
            column = param[0]
            self.galaxy[key] = data[randints][column]

        if len(joint)>0:
            randints = np.random.randint(0,high=size, size=self.ngal)
            selected = data[randints]
            for i in range(len(joint)):
                key, param = self.noindex_key(joint[i])
                column = param[0]
                self.galaxy[key] = selected[column]


    def index_key(self, arr):
        return [arr[0],arr[1],arr[2:]]

    def noindex_key(self, arr):
        return [arr[0],arr[1:]]


class Rule():

    def __init__(self, type=None, average=None, sigma=None, column=None, joint=False, value=None, array=None, component=None, minimum=None, maximum=None):

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
            self.param = [column,joint]

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

        self.type = type



class nComponentSersic(BaseSampler):
    def __init__(self, ngal=100, ncomp=2, label=None, catalog=None):
        self.init(ngal, ncomp, label, catalog)
    

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


if __name__=='__main__':

    ones = np.ones( 20 )
    disksamp = Exponential(ngal=20, catalog='/n/des/suchyta.1/des/SV/SV_clusters_project_home/coadd_products/rxj/catalogs/rxj_i.44_det_rriizz.43.42.44.40.44.40_nomodel.cat.fits')
    #disksamp.ComponentRule(key='flux', rule=Rule(type='value', value=175) )
    #disksamp.ComponentRule(key='flux', rule=Rule(type='uniform', minimum=160, maximum=180) )
    #disksamp.ComponentRule(key='flux', rule=Rule(type='catalog', column='FLUX_AUTO', joint=False) )
    #disksamp.ComponentRule(key='halflightradius', Rule(type='array', array=ones ) )
    disksamp.ComponentRule(key='flux', rule=Rule(type='catalog', column='FLUX_AUTO', joint=True) )
    disksamp.ComponentRule(key='halflightradius', rule=Rule(type='catalog', column='FLUX_RADIUS', joint=True) )
    disksamp.GalaxyRule(key='x', rule=Rule(type='uniform', minimum=0, maximum=28000) )
    disksamp.GalaxyRule(key='y', rule=Rule(type='uniform', minimum=0, maximum=28000) )
    disksamp.Sample()
    print disksamp.profile
    print disksamp.galaxy

    '''
    bdsamp = BulgeDisk(ngal=20)
    bdsamp.ComponentRule(component='bulge', key='flux', rule=Rule(type='gaussian', average=170, sigma=10) )
    bdsamp.ComponentRule(component='disk', key='flux', rule=Rule(type='component', component='bulge') )
    bdsamp.Sample()
    '''
