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
hkl = [-2,0,4]
twotheta = ubmatrix.calcTwoTheta(hkl, star, wavelength)
print('({0}, {1}, {2}) 2theta = {3:.3f}'.format(*hkl, twotheta))

# # calculate scattering plane
hkl1 = [2,0,0]
hkl2 = [0,0,2]
#chi, phi = ubmatrix.calcScatteringPlane(hkl1, hkl2, UB_matrix, wavelength, star)
#print('scattering plane of ({0}, {1}, {2}) and ({3}, {4}, {5}) chi: {6:.3f}    phi: {7:.3f}'.format(*hkl1, *hkl2, chi,phi))

plane_name = 'h0l'
hkl_max_component_val = 2
wombat_scattering_plane.accessible_hkl_in_scattering_plane(sample_name_prefix, plane_name, hkl1, hkl2, hkl_max_component_val, UB_matrix, wavelength, star, wom_stth)

# generate list of hkls
#hkl_list_to_test = wombat_scattering_plane.generate_hkl_list(hkl_limits)

# test if hkls allowed by space group
# space_group_number = 15
# hkl_allowed_list = wombat_scattering_plane.hkl_allowed(hkl_list_to_test, space_group_number)
# print('allowed reflections')
# print(hkl_allowed_list)

# calculate list of accessible hkl reflections, keeping eom fixed at zero degrees
# accessible_reflections_df = wombat_scattering_plane.accessible_hkl_omega_zero_list(sample_name_prefix, hkl_list_to_test, 
#                                                                                    UB_matrix, B_matrix, wavelength, 
#                                                                                    star, wom_stth)

# # if importing a different excel spreadsheet of reflections do so here
# # reflections_df = pd.read_excel('my_spreadsheet.xlsx')
                                                                                  
# # create script to do a radcollscan in eom at each hkl in the accessible reflections dataframe
# eom_min = -1
# eom_step = 0.2 
# num_steps = 21
# oscillations_per_step = 1
# wombat_scattering_plane.generate_hkl_eom_scan_script(accessible_reflections_df, sample_name_prefix, eom_min, eom_step, 
#                                                      num_steps, oscillations_per_step)

# # calculate ideal angles (with omega = 0)
# twotheta, theta, omega, chi, phi = ubmatrix.calcIdealAngles(hkl1, UB_matrix, B_matrix, wavelength, star)
# print('for ({0}, {1}, {2}) omega: {3:.3f}   chi: {4:.3f}    phi: {5:.3f}   2theta: {6:.3f}'.format(*hkl1, omega,chi,phi,twotheta))

# # calculate ideal angles (with omega = 0)
# twotheta, theta, omega, chi, phi = ubmatrix.calcIdealAngles(hkl2, UB_matrix, B_matrix, wavelength, star)
# print('for ({0}, {1}, {2}) omega: {3:.3f}   chi: {4:.3f}    phi: {5:.3f}   2theta: {6:.3f}'.format(*hkl2, omega,chi,phi,twotheta))

# # calculate angle between two hkls
# #interplanar_angle = np.arccos(d1 dot d2/ absd1 absd2)
# hkl1 = [-6,0,6]
# hkl2 = [0,0,2]
# interplanar_angle = np.arccos(ubmatrix.scalar(*hkl1, *hkl2, star)/(ubmatrix.modvec(*hkl1, star)*ubmatrix.modvec(*hkl2, star)))*180/np.pi
# print('interplanar angle between ({0}, {1}, {2}) and ({3}, {4}, {5}): {6:.3f} degrees'.format(*hkl1, *hkl2, interplanar_angle))

# # calculate angle between two hkls
# #interplanar_angle = np.arccos(d1 dot d2/ absd1 absd2)
# hkl1 = [-6,0,6]
# hkl2 = [-2,0,4]
# interplanar_angle = np.arccos(ubmatrix.scalar(*hkl1, *hkl2, star)/(ubmatrix.modvec(*hkl1, star)*ubmatrix.modvec(*hkl2, star)))*180/np.pi
# print('interplanar angle between ({0}, {1}, {2}) and ({3}, {4}, {5}): {6:.3f} degrees'.format(*hkl1, *hkl2, interplanar_angle))
# print('interplanar angle - 180 {0}'.format(180-interplanar_angle))