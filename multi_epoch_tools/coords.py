import astropy.wcs as apywcs
import galsim
import numpy as np
"""
def get_wcoords(x, y, header):
    wcs = pywcs.WCS(header)
    pcoords = np.dstack((x,y))[0]
    wcoords = wcs.wcs_pix2sky(pcoords, 1)
    return wcoords.T
"""
def get_wcoords(xs, ys, wcs):
    wcoords=np.zeros((2,xs.shape[0]))
    for i,(x,y) in enumerate(zip(xs,ys)):
        im_pos=galsim.PositionD(np.float(x),np.float(y))
        wcoords_galsim=wcs.toWorld(im_pos)
        wcoords[0,i],wcoords[1,i]=wcoords_galsim.ra/galsim.degrees,wcoords_galsim.dec/galsim.degrees
    return wcoords



def get_pcoords(ra, dec, wcs):
    #wcs=galsim.wcs.readFromFitsHeader(header)[0]                                                                              
    pcoords=np.zeros((2,ra.shape[0]))
    #this is pretty lame....but pywcs/astropy don't handle TPV, so using galsim. Is there a better way to                      
    #do this in galsim?                                                                                                        
    for i,(r,d) in enumerate(zip(ra,dec)):
        world_pos=galsim.celestial.CelestialCoord(r*galsim.degrees,d*galsim.degrees)
        pcoords_galsim=wcs.toImage(world_pos)
        pcoords[0,i],pcoords[1,i]=pcoords_galsim.x,pcoords_galsim.y
    return pcoords
