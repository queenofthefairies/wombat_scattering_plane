# the virtual environment to run this script in is openopt-ubmatrix

import sys
import numpy as np
import pandas as pd

# Add DerApproximator to path for openopt compatibility
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/DerApproximator')
# Add openopt to path for openopt compatibility
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/OpenOpt/openopt')
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/OpenOpt')

# Add UBmatrix to path
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/ubmatrix')

from oo import NLSP
import ubmatrix

# import wombat_scattering_plane module
import wombat_scattering_plane 

# Unit cell params for Y2SiO5
unit_cell_params = [14.406, 6.728, 10.42, 90.00, 122.1938, 90.00]

wavelength = 2.41 # in Angstrom
wom_stth = 13 # first 2theta angle of wombat detector

sample_name_prefix = 'Y2SiO5_test'
hkl_limits = [-4, 4, -4, 4, -4, 4] # [h_min, h_max, k_min, k_max, l_min, l_max]
UB_matrix = np.array([[-0.07128796726465, -0.00493851210922, -0.00514440611005],
                      [0.04053531959653, -0.00255302991718,  0.11326986551285],
                      [-0.00167354044970,  0.14852857589722,  0.00177592528053]])

# calculate star (a.k.a. reciprocal lattice params)
star = ubmatrix.star(*unit_cell_params)
star = dict(zip(('astar','bstar','cstar','alphastar','betastar','gammastar'),
                 star))
print('star (reciprocal lattice params)')
print(star)

# calculate B matrix
B_matrix = ubmatrix.calcB(star['astar'],star['astar'],star['astar'],
                          star['alphastar'],star['betastar'],star['gammastar'],
                          unit_cell_params[2], unit_cell_params[3])

# calculate 2theta of a reflection 
hkl_list = [[0,0,2],
            [-6,0,6],
            [-2,0,4],
            [-4,0,6],
            [0,0,4],
            [2,0,0],
            [0,0,2],
            [-8,0,2], 
            [-10, 0, 2], 
            [-8, 0, 4], 
            [-4, 0, 4]
            ]
for hkl in hkl_list:
    twotheta = ubmatrix.calcTwoTheta(hkl, star, wavelength)
    print('({0}, {1}, {2}) 2theta = {3:.3f}'.format(*hkl, twotheta))

# calculate angle between two hkls
#interplanar_angle = np.arccos(d1 dot d2/ absd1 absd2)
for i in range(len(hkl_list)-1):
    hkl1 = hkl_list[i]
    hkl2 = hkl_list[i+1]
    interplanar_angle = np.arccos(ubmatrix.scalar(*hkl1, *hkl2, star)/(ubmatrix.modvec(*hkl1, star)*ubmatrix.modvec(*hkl2, star)))*180/np.pi
    print('interplanar angle between ({0}, {1}, {2}) and ({3}, {4}, {5}): {6:.3f} degrees'.format(*hkl1, *hkl2, interplanar_angle))

# calculate angle between two hkls
#interplanar_angle = np.arccos(d1 dot d2/ absd1 absd2)
#hkl1 = [-6,0,6]
##hkl2 = [-2,0,4]
#interplanar_angle = np.arccos(ubmatrix.scalar(*hkl1, *hkl2, star)/(ubmatrix.modvec(*hkl1, star)*ubmatrix.modvec(*hkl2, star)))*180/np.pi
#print('interplanar angle between ({0}, {1}, {2}) and ({3}, {4}, {5}): {6:.3f} degrees'.format(*hkl1, *hkl2, interplanar_angle))
#print('interplanar angle - 180 {0}'.format(180-interplanar_angle))