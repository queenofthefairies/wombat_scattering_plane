# the virtual environment to run this script in is wombatsinglecrystal

import sys
import numpy as np

# Add DerApproximator to path for openopt compatibility
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/DerApproximator')
# Add openopt to path for openopt compatibility
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/OpenOpt/openopt')
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/OpenOpt')

# Add UBmatrix to path
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/ubmatrix')

# Add wombat_scattering_plane to path
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/siobhan_scattering_plane')

# add location of Wombat DMCpy scripts
sys.path.append('J:\wombat_instrument_work\DMCpy_for_wombat\Wombat_DMCpy')

from wombatDMCpy._tools import readCryFileFromInt3D
from oo import NLSP
import ubmatrix
import wombat_scattering_plane 

# sample name
sample_number = 2
sample_name_prefix = 'Er2SiO5_April_crystal_{0}'.format(sample_number)

# Read in unit cell params for Er2SiO5 and UB matrix found from Int3D
unit_cell_params, UB_matrix = readCryFileFromInt3D('UB7.cry')

wavelength = 2.41 # in Angstrom

# calculate star (a.k.a. reciprocal lattice params)
star = ubmatrix.star(*unit_cell_params)
star = dict(zip(('astar','bstar','cstar','alphastar','betastar','gammastar'),
                star))

# calculate B matrix
B_matrix = ubmatrix.calcB(star['astar'],star['astar'],star['astar'],
                          star['alphastar'],star['betastar'],star['gammastar'],
                          unit_cell_params[2], unit_cell_params[3])

# have a look at which scattering planes are accessible given the UB matrix
wombat_scattering_plane.evaluate_possible_scattering_planes(sample_name_prefix, UB_matrix, wavelength, star)

# calculate 2theta of a reflection 
hkl_to_calc_list = [[5,3,0],
                    [5,1,0],
                    [6,2,0],
                    [2,0,0],
                    [4,0,0],
                    [5,3,0]
                    ]
for hkl_peak in hkl_to_calc_list:
     print('hkl: ({0:.2f}, {1:.2f}, {2:.2f})   two theta: {3:.2f}'.format(hkl_peak[0],hkl_peak[1],hkl_peak[2], ubmatrix.calcTwoTheta(hkl_peak, star, wavelength)))
# given echi = 14.19 and ephi = 137.56 puts sample in the hk0 plane, find the 2thetas & omegas of the 
# reflections in hkl_to_calc_list
echi_in_plane = 14.19
ephi_in_plane = 137.56
for hkl in hkl_to_calc_list:
    wombat_scattering_plane.hkl_in_plane_omega_twotheta(hkl, UB_matrix, 
                                                        echi_in_plane, 
                                                        ephi_in_plane, 
                                                        wavelength, star)

# for the hk0 scattering plane (which is apparently accessible), 
# see which reflections are visible on Wombat given the Q limits
# (+-3, +-3, 0) 
# i.e. hkl_max_component_val = 3
# plane_name = 'hk0'
# hkl1 = [1, 0, 0]
# hkl2 = [0, 1, 0]
# hkl_max_component_val = 6 
# wom_stth = 13
# wombat_scattering_plane.accessible_hkl_in_scattering_plane(sample_name_prefix, plane_name, 
#                                                            hkl1, hkl2, hkl_max_component_val,
#                                                            UB_matrix, wavelength, star, wom_stth)

