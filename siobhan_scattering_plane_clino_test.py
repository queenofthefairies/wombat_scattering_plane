# the virtual environment to run this script in is openopt-ubmatrix

import sys
import numpy as np

# Add DerApproximator to path for openopt compatibility
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/DerApproximator')
# Add openopt to path for openopt compatibility
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/OpenOpt/openopt')
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/OOSuite/OpenOpt')

# Add UBmatrix to path
sys.path.append('J:/wombat_instrument_work/scattering_plane_software/ubmatrix')

#import oo as openopt
#from openopt import NLSP
from oo import NLSP
import ubmatrix

import wombat_scattering_plane 

# sample name
sample_name_prefix = 'clino'

# Unit cell params for clinoatacamite
unit_cell_params = [6.144, 6.805, 9.112, 90.00, 99.55, 90.00]
wavelength = 2.41 # in Angstrom
# UB matrix found with Int3D
UB_matrix = np.array([[-0.04263322800398, 0.13622714579105, 0.02506458014250],
                      [0.15944631397724, 0.03636400774121, 0.02596857771277],
                      [0.00023406882246, 0.04140481352806, -0.10527275502682]])

# calculate star (a.k.a. reciprocal lattice params)
star = ubmatrix.star(*unit_cell_params)
star = dict(zip(('astar','bstar','cstar','alphastar','betastar','gammastar'),
                star))
print('star')
print(star)

# calculate B matrix
B_matrix = ubmatrix.calcB(star['astar'],star['astar'],star['astar'],
                          star['alphastar'],star['betastar'],star['gammastar'],
                          unit_cell_params[2], unit_cell_params[3])

# calculate 2theta of a reflection 
hkl = [-1,0,5]
print('hkl: ({0}, {1}, {2})'.format(hkl[0],hkl[1],hkl[2]))
print(ubmatrix.calcTwoTheta(hkl, star, wavelength))

# have a look at which scattering planes are accessible given the UB matrix
wombat_scattering_plane.evaluate_possible_scattering_planes(sample_name_prefix, UB_matrix, wavelength, star)

# for the hk0 scattering plane (which is apparently accessible), 
# see which reflections are visible on Wombat given the Q limits
# (+-3, +-3, 0) 
# i.e. hkl_max_component_val = 3
plane_name = 'hk0'
hkl1 = [1, 0, 0]
hkl2 = [0, 1, 0]
hkl_max_component_val = 3
wom_stth = 13
wombat_scattering_plane.accessible_hkl_in_scattering_plane(sample_name_prefix, plane_name, 
                                                           hkl1, hkl2, hkl_max_component_val,
                                                           UB_matrix, wavelength, star, wom_stth)

plane_name = 'hkk'
hkl1 = [1, 0, 0]
hkl2 = [0, 1, 1]
hkl_max_component_val = 3
wom_stth = 13
wombat_scattering_plane.accessible_hkl_in_scattering_plane(sample_name_prefix, plane_name, 
                                                           hkl1, hkl2, hkl_max_component_val,
                                                           UB_matrix, wavelength, star, wom_stth)

# calculate angle between two hkls
#interplanar_angle = np.arccos(d1 dot d2/ absd1 absd2)
hkl1 = [1,2,2]
hkl2 = [2,3,3]
interplanar_angle = np.arccos(ubmatrix.scalar(*hkl1, *hkl2, star)/(ubmatrix.modvec(*hkl1, star)*ubmatrix.modvec(*hkl2, star)))*180/np.pi
print('interplanar angle between ({0}, {1}, {2}) and ({3}, {4}, {5}): {6:.3f} degrees'.format(*hkl1, *hkl2, interplanar_angle))

# calculate angle between two hkls
#interplanar_angle = np.arccos(d1 dot d2/ absd1 absd2)
hkl1 = [-1,2,2]
hkl2 = [0,2,2]
interplanar_angle = np.arccos(ubmatrix.scalar(*hkl1, *hkl2, star)/(ubmatrix.modvec(*hkl1, star)*ubmatrix.modvec(*hkl2, star)))*180/np.pi
print('interplanar angle between ({0}, {1}, {2}) and ({3}, {4}, {5}): {6:.3f} degrees'.format(*hkl1, *hkl2, interplanar_angle))
