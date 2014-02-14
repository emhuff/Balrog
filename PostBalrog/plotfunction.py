#!/usr/bin/env python

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pyfits


def MakeScatter(*arguments, **keywords):
    plt.figure(keywords['fig'])

    if keywords['cut']!=None:
        keywords['measured'] = keywords['measured'][keywords['cut']]
        keywords['truth'] = keywords['truth'][keywords['cut']]
        if keywords['color']!=None:
            keywords['color'] = keywords['color'][keywords['cut']]
    length = len(keywords['truth'])

    if 'weight' not in keywords.keys():
        sorted = np.zeros( length, dtype=[('x','f8'), ('y','f8')] )
        sorted['x'] = keywords['truth']
        sorted['y'] = keywords['measured']
    else:
        sorted = np.zeros( length, dtype=[('x','f8'), ('y','f8'), ('w','f8')] )
        sorted['x'] = keywords['truth']
        sorted['y'] = keywords['measured']
        sorted['w'] = keywords['weight']

    high = 1 - keywords['low']
    if keywords['xlim']==None:
        sorted = np.sort(sorted, order=['x'])
        imin = int( length*keywords['low'] )
        imax = int( length*high-1 )
        xmin = sorted['x'][imin]
        xmax = sorted['x'][imax]
        keywords['xlim'] = [xmin,xmax]
    if keywords['ylim']==None:
        sorted = np.sort(sorted, order=['y'])
        imin = int( length*keywords['low'] )
        imax = int( length*high-1 )
        ymin = sorted['y'][imin]
        ymax = sorted['y'][imax]
        keywords['ylim'] = [ymin,ymax]

    if keywords['bin']:
        sorted = np.sort(sorted, order=['x'])
        xcut = (sorted['x'] >= xmin) & (sorted['x'] <= xmax)
        cut = sorted[xcut]
        newlength = len(cut)
        perbin = 1000

        xbins = []
        ybins = []
        elows = []
        ehighs = []
        index = 0
        while index <= newlength-perbin:
            xs = cut['x'][index:(index+perbin)]
            ys = cut['y'][index:(index+perbin)]
            if 'weight' not in keywords.keys():
                xbin = np.average( xs )
                ybin = np.average( ys )
                '''
                xbin = np.median( xs )
                ybin = np.median( ys )
                '''
            else:
                ws = cut['w'][index:(index+perbin)]
                xbin = np.average( xs, weights=ws )
                ybin = np.average( ys, weights=ws )
            upcut = (ys >= ybin)
            lowcut = (ys < ybin)
            '''
            ehigh = np.std( ys[upcut] )
            elow = np.std( ys[lowcut] )
            '''
            ehigh = np.std( ys ) / np.sqrt(perbin)
            elow = np.std( ys ) / np.sqrt(perbin)
            xbins.append(xbin)
            ybins.append(ybin)
            elows.append(elow)
            ehighs.append(ehigh)
            index += perbin
        ebins = np.array([elows,ehighs])

    nx = 100
    if keywords['line']==0:
        linex = np.linspace( keywords['xlim'][0], keywords['xlim'][1], num=nx )
        liney = np.zeros( len(linex) )
    else:
        linex = np.linspace( keywords['xlim'][0], keywords['xlim'][1], num=nx )
        liney = linex

    if keywords['logx']:
        keywords['truth'] = np.log10(keywords['truth'])
        linex = np.log10(linex)
        keywords['xlim'] = np.log10(keywords['xlim'])
        if keywords['bin']:
            xbins = np.log10(xbins)
    if keywords['logy']:
        keywords['measured'] = np.log10(keywords['measured'])
        liney = np.log10(liney)
        keywords['ylim'] = np.log10(keywords['ylim'])
        if keywords['bin']:
            ebins = (1.0 / np.log(10)) * (ebins / ybins)
            ybins = np.log10(ybins)


    if keywords['color']==None:
        plt.scatter(keywords['truth'], keywords['measured'], s=keywords['size'], lw=0, c='red')
    else:
        plt.scatter(keywords['truth'], keywords['measured'], s=keywords['size'], lw=0, c=keywords['color'])
        pc = plt.colorbar()
        pc.set_label(keywords['clabel'])

    if keywords['bin']:
        plt.scatter(xbins, ybins, c='blue', s=5)
        plt.errorbar(xbins, ybins, yerr=ebins, fmt='o')
    plt.plot(linex,liney, c='black', lw=0.5)
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel(keywords['xlabel'], fontsize=16)
    if keywords['ylarger']:
        plt.ylabel(keywords['ylabel'], fontsize=22)
    else:
        plt.ylabel(keywords['ylabel'], fontsize=16)
    plt.xlim(keywords['xlim'])
    plt.ylim(keywords['ylim'])
    plt.figtext(0.75, 0.8, keywords['text'])
    plt.title(keywords['title'])


def PrepScatter(*arguments, **keywords):
    keywords['truth'] = keywords['tcat'][keywords['tkey']]
    if keywords['err']:
        keywords['measured'] = (keywords['mcat'][keywords['mkey']] - keywords['truth']) / keywords['truth']
    else:
        keywords['measured'] = keywords['mcat'][keywords['mkey']]

    keywords['truth'] = keywords['xcat'][keywords['tx']]

    if keywords['weightkey']!=None:
        keywords['weight'] = keywords['mcat'][keywords['weightkey']]
    MakeScatter(**keywords)


def DefaultKwargs():
    defaults = {}
    defaults['tcat']=None
    defaults['mcat']=None
    defaults['xcat']=None

    defaults['x']='flux'
    defaults['xlabel']='F'
    defaults['xunit']=''
    defaults['logx']=False
    defaults['xlim']=None

    defaults['y']='flux'
    defaults['ylabel']='F'
    defaults['yunit']=''
    defaults['logy']=False
    defaults['ylim']=None
    defaults['ylarger']=False

    defaults['color']=None
    defaults['clabel']=''
    defaults['size']=0.5

    defaults['err']=False
    defaults['bin']=False
    defaults['weightkey']=None
    defaults['cut']=None
    defaults['low']=0.005

    defaults['fig']=1
    defaults['text']=''
    defaults['title']=''
    defaults['meas']='meas'

    return defaults


def GetKeywords(*arguments, **keywords):
    default = DefaultKwargs()
    for key in default.keys():
        if key not in keywords:
            keywords[key] = default[key]
    return keywords


def MakePlot(**keywords):
    keywords = GetKeywords(**keywords) 
    if keywords['err']:
        n = "%s_{%s} - %s_{truth}" %(keywords['ylabel'],keywords['meas'],keywords['ylabel'])
        d = "%s_{truth}" %(keywords['ylabel'])
        keywords['ylabel'] = "\\frac{%s}{%s}" %(n,d)
    else:
        keywords['ylabel'] = "%s_{%s}" %(keywords['ylabel'], keywords['meas'])
    if not keywords['yunit']=='':
        keywords['ylabel'] = "%s/[%s]" %(keywords['ylabel'], keywords['yunit'])
    if keywords['logy']:
        keywords['ylabel'] = "\\log_{10} \\left( %s \\right)" %(keywords['ylabel'])
    keywords['ylabel'] = r"$%s$" %(keywords['ylabel'])

    if keywords['x']!='s/n':
        keywords['xlabel'] = "%s_{truth}" %(keywords['xlabel'])
    if not keywords['xunit']=='':
        keywords['xlabel'] = "%s/[%s]" %(keywords['xlabel'], keywords['xunit'])
    if keywords['logx']:
        keywords['xlabel'] = "\\log_{10} \\left( %s \\right)" %(keywords['xlabel'])
    keywords['xlabel'] = r"$%s$" %(keywords['xlabel'])

    if keywords['err']:
        '''
        if ylim==None:
            ylim = [-1,2]
        '''
        keywords['line'] = 0
    else:
        keywords['line'] = 1

    keywords['mkey'] = keywords['y']
    keywords['tkey'] = keywords['y']
    keywords['tx'] = keywords['x']
    PrepScatter(**keywords)


