#!/usr/bin/env python

import utils
import numpy as np
import subprocess


if __name__=='__main__':

    tiles = ['DES0419-4831',
             'DES0423-4831',
             'DES0427-4831',
             'DES0432-4831',
             'DES0436-4831',
             'DES0440-4831',
             'DES0445-4831',
             'DES0449-4831',
             'DES0453-4831',
             'DES0458-4831',
             'DES0502-4831']

    for tile in tiles:
        subprocess.call( ['./fetch_tile.py', '--band', 'i', '--tile', tile] )
