# the virtual environment to run this script in is wombatsinglecrystal

import sys
import numpy as np

# Add DerApproximator to path for openopt compatibility
sys.path.append('C:/Users/wombat/OOSuite/DerApproximator')
# Add openopt to path for openopt compatibility
sys.path.append('C:/Users/wombat/OOSuite/OpenOpt/openopt')
sys.path.append('C:/Users/wombat/OOSuite/OpenOpt')

# Add UBmatrix to path
sys.path.append('C:/Users/wombat/ubmatrix')

# Add wombat_scattering_plane to path
sys.path.append('C:/Users/wombat/wombat_scattering_plane')

# add location of Wombat DMCpy scripts
sys.path.append('C:/Users/wombat/WombatDMCpy')

from wombatDMCpy._tools import readCryFileFromInt3D
from oo import NLSP
import ubmatrix
import wombat_scattering_plane 

# sample name
sample_name_prefix = 'Y2SiO5_crystal'

# Read in unit cell params for Y2SiO5 and UB matrix found from Int3D
unit_cell_params, UB_matrix = readCryFileFromInt3D('example_UB_fromInt3D_Y2SiO5.cry')

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
hkl_to_calc_list = [[-5,0,1],
                    [-5,0,3],
                    [-6,0,2],
                    [-2,0,2],
                    [-4,0,4],
                    [-4,0,2]
                    ]
# given echi = 2.14 deg and ephi = 27.33 deg puts sample in the h0l plane, find the 2thetas & omegas of the 
# reflections in hkl_to_calc_list
echi_in_plane = 2.14
ephi_in_plane = 27.33
for hkl in hkl_to_calc_list:
    wombat_scattering_plane.hkl_in_plane_omega_twotheta(hkl, UB_matrix, 
                                                        echi_in_plane, 
                                                        ephi_in_plane, 
                                                        wavelength, star)

# for the h0l scattering plane (which is apparently accessible), 
# see which reflections are visible on Wombat given the Q limits
# (+-3, +-3, 0) 
# i.e. hkl_max_component_val = 3
# plane_name = 'h0l'
# hkl1 = [1, 0, 0]
# hkl2 = [0, 0, 1]
# hkl_max_component_val = 6 
# wom_stth = 13
# wombat_scattering_plane.accessible_hkl_in_scattering_plane(sample_name_prefix, plane_name, 
#                                                            hkl1, hkl2, hkl_max_component_val,
#                                                            UB_matrix, wavelength, star, wom_stth)

