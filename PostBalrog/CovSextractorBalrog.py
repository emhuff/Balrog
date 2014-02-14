#!/usr/bin/env python

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from plotfunction import *


def getmodel(m, data):
    model = {}
    for param in m:
        model[param[0]] = data[param[1]]
    return model


def GetTruthModelBulgeDisk(file, zp=30):
    hdus = pyfits.open(file)
    head = hdus[1].header
    data = hdus[1].data
    
    truth = {}
    index = 0
    assockeys = {}
    while True:
        key = 'V%i' %index
        if key not in head.keys():
            break
        assockeys[head[key]] = index
        index += 1
    for key in assockeys:
        pos = assockeys[key]
        k = key.replace('_0', '')
        truth[k] = data['VECTOR_ASSOC'][:,pos]
    truth['halflightradius'] = truth['halflightradius'] / np.sqrt(truth['axisratio'])
    truth['mag'] = -2.5 * np.log10(truth['flux']) + zp
    truth['surfacebrightness'] = truth['mag'] / (np.pi * np.power(truth['halflightradius'],2.0))

    pb = [('flux','FLUX_SPHEROID'),
          ('fluxerr', 'FLUXERR_SPHEROID'),
          ('mag', 'MAG_SPHEROID'),
          ('halflightradius', 'SPHEROID_REFF_WORLD'),
          ('beta', 'SPHEROID_THETA_IMAGE'),
          ('axisratio', 'SPHEROID_ASPECT_IMAGE'),
          ('sersicindex', 'SPHEROID_SERSICN')]
    bulge = getmodel(pb, data)
    bulge['halflightradius'] = bulge['halflightradius'] * 3600.0
    bulge['surfacebrightness'] = bulge['mag'] / (np.pi * np.power(bulge['halflightradius'],2.0))
    ok = (bulge['fluxerr'] > 0)
    bulge['s/n'] = np.zeros(len(ok))
    bulge['s/n'][ok] = bulge['flux'][ok] / bulge['fluxerr'][ok]

    pd = [('flux','FLUX_DISK'),
          ('fluxerr', 'FLUXERR_DISK'),
          ('mag', 'MAG_DISK'),
          ('halflightradius', 'DISK_SCALE_WORLD'),
          ('beta', 'DISK_THETA_IMAGE'),
          ('axisratio', 'DISK_ASPECT_IMAGE')]
    disk = getmodel(pd, data)
    disk['halflightradius'] = disk['halflightradius'] * 1.678 * 3600.0
    disk['surfacebrightness'] = disk['mag'] / (np.pi * np.power(disk['halflightradius'],2.0))
    disk['sersicindex'] = np.ones( len(disk['flux']) )
    ok = (disk['fluxerr'] > 0)
    disk['s/n'] = np.zeros(len(ok))
    disk['s/n'][ok] = disk['flux'][ok] / disk['fluxerr'][ok]

    pm = [('x','XMODEL_IMAGE'),
          ('y','YMODEL_IMAGE'),
          ('ra','ALPHAMODEL_J2000'),
          ('dec','DELTAMODEL_J2000'),
          ('flux', 'FLUX_MODEL'),
          ('fluxerr','FLUXERR_MODEL'),
          ('mag', 'MAG_MODEL'),
          ('a','AMODEL_IMAGE'),
          ('b','BMODEL_IMAGE'),
          ('e1','ELLIP1MODEL_IMAGE'),
          ('e2','ELLIP2MODEL_IMAGE'),
          ('flags', 'FLAGS'),
          ('wflags', 'FLAGS_WEIGHT'),
          ('dflags', 'FLAGS_DETMODEL'),
          ('flux','FLUX_MODEL'),
          ('beta','THETAMODEL_IMAGE'),
          ('mumax', 'MU_MAX_MODEL'),
          ('mueff', 'MU_EFF_MODEL'),
          ('mumean', 'MU_MEAN_MODEL'),
          ('R50', 'FLUX_RADIUS')]
    model = getmodel(pm, data)
    model['axisratio'] = model['b'] / model['a']
    ok = (model['fluxerr'] > 0)
    model['s/n'] = np.zeros(len(ok))
    model['s/n'][ok] = model['flux'][ok] / model['fluxerr'][ok]
    model['b/d'] = bulge['flux'] / disk['flux']
    
    return truth, model, bulge, disk


def single_covariance(m1, m2, t1, t2):
    num = len(m1)
    f1 = m1 - t1
    f2 = m2 - t2
    #s = np.sum(f1*f2)
    s = np.sum(f1/t1 * f2/t2)
    cov = s / num
    return cov


def compute_covariance( models_truth, sn, cut=None, sn_low=30, sn_high=40 ):
    sncut = (sn > sn_low) & (sn < sn_high)
    cuts = (cut) & (sncut)

    size = len(models_truth)
    cov = np.empty( (size,size) )
    reg = np.zeros( (size,size) )
    labels = []
    for i in range(size):
        for j in range(i,size):
            m1 = models_truth[i][0][cuts]
            m2 = models_truth[j][0][cuts]
            t1 = models_truth[i][1][cuts]
            t2 = models_truth[j][1][cuts]

            c = single_covariance(m1,m2,t1,t2)
            cov[i][j] = c
            reg[i][j] = 1
            cov[j][i] = c
            if j!=i:
                reg[j][i] = 0

        labels.append( models_truth[i][2] )

    return cov, labels, cuts, np.bool_(reg)


def CovPlot(*args, **kwargs):
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    if kwargs['hist']!=None:
        plt.figure(kwargs['fig']+1)
        bins = np.arange(0,8.5, 0.15)
        plt.hist(kwargs['type']['sersicindex'][kwargs['cut']], bins=bins)
    
    fig = plt.figure(kwargs['fig'], figsize=(10,8))
    fig = plt.gcf()
    fig.suptitle(kwargs['title'])

    plt.subplot(221)
    MakePlot(cut=kwargs['cut'], tcat=kwargs['truth'], mcat=kwargs['type'], y='halflightradius', ylabel='R^e', xcat=kwargs['truth'], x='halflightradius', xlabel='R^e', fig=kwargs['fig'], err=True, logx=False, logy=False, color=kwargs['color'], clabel=kwargs['clabel'], size=5, low=0.05, meas=kwargs['meas'], xunit='arcsec', yunit='', ylarger=True)

    plt.subplot(222)
    MakePlot(cut=kwargs['cut'], tcat=kwargs['truth'], mcat=kwargs['type'], y='mag', ylabel='M^i', xcat=kwargs['truth'], x='mag', xlabel='M^i', fig=kwargs['fig'], err=True, logx=False, logy=False, color=kwargs['color'], clabel=kwargs['clabel'], size=5, low=0.05, meas=kwargs['meas'], xunit='mag', yunit='', ylarger=True)

    plt.subplot(223)
    MakePlot(cut=kwargs['cut'], tcat=kwargs['truth'], mcat=kwargs['type'], y='surfacebrightness', ylabel='\\mu', xcat=kwargs['truth'], x='surfacebrightness', xlabel='\\mu', fig=kwargs['fig'], err=True, logx=True, logy=False, color=kwargs['color'], clabel=kwargs['clabel'], size=5, low=0.05, meas=kwargs['meas'], xunit='\\frac{mag}{arcsec^2}',yunit='', ylarger=True )
    #MakePlot(cut=kwargs['cut'], tcat=kwargs['truth'], mcat=kwargs['type'], y='sersicindex', ylabel='n', xcat=kwargs['truth'], x='sersicindex', xlabel='n', fig=kwargs['fig'], err=True, logx=False, logy=False, color=sn, clabel='S/N', size=5, xlim=[3.8, 4.2], meas=kwargs['meas'], xunit='',yunit='', ylarger=True )

    plt.subplot(224)
    MakePlot(cut=kwargs['cut'], tcat=kwargs['truth'], mcat=kwargs['type'], y='axisratio', ylabel='(b/a)', xcat=kwargs['truth'], x='axisratio', xlabel='(b/a)', fig=kwargs['fig'], err=False, logx=False, logy=False, color=kwargs['color'], clabel=kwargs['clabel'], size=5, xlim=[0.33,1],ylim=[0.33,1], meas=kwargs['meas'], xunit='',yunit='', ylarger=False )
    '''
    plt.xticks( np.arange(0.5,len(kwargs['labels'])+0.5,1) )
    plt.yticks( np.arange(0.5,len(kwargs['labels'])+0.5,1) )
    plt.gca().set_xticklabels(kwargs['labels'], fontsize=16)
    plt.gca().set_yticklabels(kwargs['labels'], fontsize=16)
    max = np.amax(np.abs(kwargs['cov'].flatten()))
    #cov[-reg] = 0
    plt.pcolor(kwargs['cov'], cmap=plt.get_cmap('RdBu'), vmin=-max, vmax=max)
    pc = plt.colorbar()
    pc.set_label('Covariance')
    '''

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()
    #plt.savefig('n1_example.png')


if __name__ == "__main__":

    lowsn = 700
    highsn = 1500
    rcut = 0.3
    mcut = 26

    cat = sys.argv[1]
    image = sys.argv[2]
    zp = pyfits.open(image)[0].header['AVG_ZP']
    truth, model, bulge, disk = GetTruthModelBulgeDisk(cat, zp=zp)
    flagcut = (model['flags']==0) & (model['wflags']==0) & (model['dflags']==0) & (model['R50']*0.263 > rcut)
 

    t = 'disk'
    if t=='disk':
        type = disk
        nt = 1
        histo = None
    elif t=='bulge':
        type = bulge
        nt = 4
        histo = nt

    typecut = (truth['sersicindex']==nt)
    sizecut = (type['halflightradius'] > 0*rcut)
    magcut = (type['mag'] < mcut )
    datacuts = flagcut & sizecut & magcut
    typecuts = datacuts & typecut
    sn = model['s/n']
    typesn = type['s/n']

    models_truth = [(model['mag'],truth['mag'], r'$M^i_{disk+bulge}$'),
                    (type['mag'],truth['mag'], '$M^i_{%s}$'%(t)),
                    (type['halflightradius'],truth['halflightradius'], '$R^e_{%s}$'%(t)),
                    (type['surfacebrightness'],truth['surfacebrightness'], '$\mu_{%s}$'%(t))]
    #cov, labels, sncuts, reg = compute_covariance( models_truth, sn, sn_low=lowsn, sn_high=highsn, cut=typecuts )
    cov, labels, sncuts, reg = compute_covariance( models_truth, typesn, sn_low=lowsn, sn_high=highsn, cut=typecuts )
    allcuts = typecuts&sncuts

    title = r'$n_{truth}=%i$'%(nt) + '\n' + r'$R^e_{%s} > %.1f~arcsec$'%(t,rcut) + '\n' + r'$M^i_{%s} < %.1f$'%(t,mcut) + '\n' + r'$FLAGS==0$' + '\n' + r'$%i < S/N < %i$'%(lowsn,highsn)
    CovPlot(title=title, cut=allcuts, truth=truth, type=type, cov=cov, labels=labels, fig=1, hist=histo, color=typesn, clabel='S/N', meas=t)
    #CovPlot(title=title, cut=allcuts, truth=truth, type=type, cov=cov, labels=labels, fig=1, hist=histo, color=np.log10(bulge['flux']/disk['flux']), clabel='B/D', meas=t)



    ''' 
    unspec = {}
    for key in disk.keys():
        unspec[key] = np.append(disk[key],bulge[key])
    for key in truth.keys():
        truth[key] = np.append(truth[key],truth[key])
    usn = unspec['s/n']

    sizecut = (unspec['halflightradius'] > rcut*0)
    magcut = (unspec['mag'] < mcut ) & (unspec['mag'] > 19)
    fcut = np.append(flagcut,flagcut)
    datacuts = fcut & sizecut & magcut

    models_truth = [(unspec['mag'],truth['mag'], '$M^i$'),
                    (unspec['halflightradius'],truth['halflightradius'], '$R^e$'),
                    (unspec['surfacebrightness'],truth['surfacebrightness'], '$\mu$')]
    cov, labels, sncuts, reg = compute_covariance( models_truth, usn, sn_low=lowsn, sn_high=highsn, cut=datacuts )
    allcuts = datacuts&sncuts

    title = r'$R^e > %.1f~arcsec$'%(rcut) + '\n' + r'$M^i < %.1f$'%(mcut) + '\n' + r'$FLAGS==0$' + '\n' + r'$%i < S/N < %i$'%(lowsn,highsn)
    print np.amin(usn)
    CovPlot(title=title, cut=allcuts, truth=truth, type=unspec, cov=cov, labels=labels, fig=1, hist=None, color=usn, clabel='S/N', meas=None)
    '''


    '''
    models_truth = [#(disk['axisratio'],truth['axisratio']),
                    #(disk['beta'],truth['beta']),
                    #(disk['sersicindex'],truth['sersicindex']),
                    
                    #(bulge['axisratio'],truth['axisratio']),
                    #(bulge['beta'],truth['beta']),
                    #(bulge['sersicindex'],truth['sersicindex']),

                    #(model['beta'],truth['beta']),
                    #(model['mumax'],truth['surfacebrightness']),
                    #(model['mumean'],truth['surfacebrightness']),
                    #(model['mueff'],truth['surfacebrightness']),
                    #(model['axisratio'],truth['axisratio']),
    '''
