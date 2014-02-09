#!/usr/bin/env python

from model_class import *


###In lieu of having time for proper documentation, for now I'll give you some commented examples, which you can extrapolate from.

'''
This function controls how your galaxies will be simulated. If you're using a catalog to sample from, e.g. the COMSOS
catalog, that's specfied in a command line option (--catalogsample) in balrog.

The overall "default" behavior (based on my command line arguments) is to sample from the COSMOS catalog. If that's what
you are doing you won't need to change much. However if you do want to simulate otherwise, these classes
make that very straightforward for you to do. You could possibly even add extra command line arguments if you wanted
and then use those in the function below instead of the command line options I'm using.
There's a function GetArgs() in balrog.py. I'm not going to copy that into here, because don't change it if you don't need to.
There's less chance of screwing up that way.

GalaxyRule is for things the galaxy only has one of: ['x', 'y', 'g1', 'g2', 'magnification'].
    g1, g2 are lensing shears
    The default is no lensing:  g1 = g2 = 0,  magnification = 1.
    Default positions are random.

ComponentRule is for aspects of the Sersic profile models. You can model a galaxy with as many superimposed sersic profiles as you want.
    ncomp --> The number of sersic profiles, e.g. bulge + disk would have ncomp = 2
    key --> Which parameter of the sersic profile: ['mag', 'halflightradius', 'sersicindex', 'axisratio', 'beta']
            They all come with defaults, and you overwrite them however you want.  
            You won't want the default mag, halflightradius, or sersicindex since there aren't obvious defaults to use here.
            (axisratio = b/a, default is uniform between 1/3 and 1/1 (galsim didn't seem to like to draw very elliptial))
            (beta is the angle the major axis makes with x-direction, default is uniform between -90 and 90 degrees)
    type --> How you want to sample your galaxies. The other parameters should make sense then give this type).
             I've given an example of all the types I've implemented so far (uniform prob., gaussian, one value, sampling from
             a catalog, or the catch all where the user gives an array. Catalog sampling can be done jointly or independtly
             between different parameters, but more than likely you want joint)

There are shortcut classes is you want, e.g. see the end of model_class.py for all of them:
    deVaucouleur()
    Exponential()
    BulgeDisk()

nComponentSerisic() can do any of the above
'''

def SimulationRules(opts):
    # If there's multiple components, you can either give an integer for which component you mean or by defined keywords.
    # For BulgeDisk 0 = 'bulge', 1 = 'disk'
    # The default is 0, so if it's a single component Serisc you don't need the component keyword

    #simulatedgals = BulgeDisk(ngal=5, catalog='sample_from_here.fits')
    #simulatedgals.ComponentRule(component='bulge', key='mag', rule=Rule(type='catalog', column='IMAG', joint=True) )
    #simulatedgals.ComponentRule(component=0, key='halflightradius', rule=Rule(type='catalog', column='R50', joint=True) )
    #simulatedgals.ComponentRule(component='bulge', key='sersicindex', rule=Rule(type='value', value=4) )
    #simulatedgals.ComponentRule(component='disk', key='mag', rule=Rule(type='gaussian', average=23, sigma=0.5) )
    #simulatedgals.ComponentRule(component=1, key='halflightradius', rule=Rule(type='array', array=np.arange(5)/2.0) )
    #simulatedgals.ComponentRule(component='disk', key='sersicindex', rule=Rule(type='value', value=1) )


    simulatedgals = nComponentSersic(ngal=opts.ngal, catalog=opts.catalogsample, ncomp=1)
    simulatedgals.GalaxyRule(key='x', rule=Rule(type='uniform', minimum=opts.xmin, maximum=opts.xmax) )
    simulatedgals.GalaxyRule(key='y', rule=Rule(type='uniform', minimum=opts.ymin, maximum=opts.ymax) )
    simulatedgals.ComponentRule(key='axisratio',rule=Rule(type='uniform', minimum=0.33, maximum=1) )
    simulatedgals.ComponentRule(key='beta',rule=Rule(type='uniform', minimum=-90, maximum=90 ) )
    simulatedgals.ComponentRule(key='halflightradius', rule=Rule(type='catalog', column=opts.reff, joint=True) )
    simulatedgals.ComponentRule(key='mag', rule=Rule(type='catalog', column=opts.mag, joint=True) )
    simulatedgals.ComponentRule(key='sersicindex', rule=Rule(type='catalog', column=opts.sersicindex, joint=True) )
    return simulatedgals


'''
This is extra configurations to give to sextractor other than the config file in astro_config/
'''
def SextractorConfigs():
    config = {}
    config['CHECKIMAGE_TYPE'] = 'NONE'
    return config


