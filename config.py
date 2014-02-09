#!/usr/bin/env python

from model_class import *


def SimulationRules(opts):
    simulatedgals = nComponentSersic(ngal=opts.ngal, catalog=opts.catalogsample, ncomp=1)
    simulatedgals.GalaxyRule(key='x', rule=Rule(type='uniform', minimum=opts.xmin, maximum=opts.xmax) )
    simulatedgals.GalaxyRule(key='y', rule=Rule(type='uniform', minimum=opts.ymin, maximum=opts.ymax) )
    simulatedgals.ComponentRule(key='axisratio',rule=Rule(type='uniform', minimum=0.33, maximum=1) )
    simulatedgals.ComponentRule(key='beta',rule=Rule(type='uniform', minimum=-90, maximum=90 ) )
    simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='catalog', column=opts.reff, joint=True) )
    simulatedgals.ComponentRule(key='mag', rule=Rule(type='catalog', column=opts.mag, joint=True) )
    simulatedgals.ComponentRule(key='sersicindex', rule=Rule(type='catalog', column=opts.sersicindex, joint=True) )
    return simulatedgals


def SextractorConfigs():
    config = {}
    config['CHECKIMAGE_TYPE'] = 'NONE'
    return config
