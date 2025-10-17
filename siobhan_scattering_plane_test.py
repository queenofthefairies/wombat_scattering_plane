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

# Unit cell params for clinoatacamite
clino_unit_cell_params = [6.144, 6.805, 9.112, 90.00, 99.55, 90.00]
wavelength = 2.41 # in Angstrom
# UB matrix found with Int3D
clino_UBmatrix = np.array([[-0.04263322800398, 0.13622714579105, 0.02506458014250],
                           [0.15944631397724, 0.03636400774121, 0.02596857771277],
                           [0.00023406882246, 0.04140481352806, -0.10527275502682]])

# calculate star (a.k.a. reciprocal lattice params)
clino_star = ubmatrix.star(*clino_unit_cell_params)
clino_star = dict(zip(('astar','bstar','cstar','alphastar','betastar','gammastar'),
                         clino_star))
print('clino star')
print(clino_star)

# calculate B matrix
clino_Bmatrix = ubmatrix.calcB(clino_star['astar'],clino_star['astar'],clino_star['astar'],
                               clino_star['alphastar'],clino_star['betastar'],clino_star['gammastar'],
                               clino_unit_cell_params[2], clino_unit_cell_params[3])

# calculate 2theta of a reflection 
hkl = [-1,0,5]
print('hkl: ({0}, {1}, {2})'.format(hkl[0],hkl[1],hkl[2]))
print(ubmatrix.calcTwoTheta(hkl, clino_star, wavelength))

# calculate scattering plane
h1 = [0,0,1]
h2 = [0,1,0]
chi, phi = ubmatrix.calcScatteringPlane(h1, h2, clino_UBmatrix, wavelength,clino_star)
print('chi: {0}    phi: {1}'.format(chi,phi))

# calculate scattering plane
h1 = [0,0,1]
h2 = [0,-1,0]
chi, phi = ubmatrix.calcScatteringPlane(h1, h2, clino_UBmatrix, wavelength,clino_star)
print('chi: {0}    phi: {1}'.format(chi,phi))

# calculate ideal angles (with omega = 0)
twotheta, theta, omega, chi, phi = ubmatrix.calcIdealAngles(h1, clino_UBmatrix, clino_Bmatrix, wavelength, clino_star)
print('omega: {0}   chi: {1}    phi: {2}   2theta: {3}'.format(omega,chi,phi,twotheta))