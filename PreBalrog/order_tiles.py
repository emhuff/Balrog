#!/usr/bin/env python

import utils
import numpy as np


if __name__=='__main__':

    cur = utils.get_cursor()
    cur.execute("select ra, dec, tilename from sva1_coadd_file where band='i'")
    all = np.array( cur.fetchall() )
    ra = all[:,0]
    dec = all[:,1]
    tile = all[:,2]

    arr = np.zeros( len(ra), dtype=[('ra','f4'),('dec','f4'),('tile','a12')])
    arr['ra'] = ra
    arr['dec'] = dec
    arr['tile'] = tile

    sorted = np.sort(arr, order=['dec','ra'])
    for sort in sorted:
        if sort['tile'].find('-4831')!=-1:
            if sort['ra'] > 60:
                print sort['tile']
