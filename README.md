# Scripts for orienting crystal to a specific scattering plane for Wombat
This utilises the modules `OpenOpt` and `UBmatrix` (courtesy of William Ratcliff)    

## Functions
`UBmatrix` functions that I am testing:

### star
- input: lattice params      
- output: reciprocal lattice params     

### calcTwoTheta
- input: hkl vector, wavelength, reciprocal lattice params  
- output: 2theta of hkl of interest  

### calcScatteringPlane      
- input: two hkl vectors, UB matrix, wavelength, reciprocal lattice params   
- output: chi and phi for scattering plane  
- use: you can do an omega scan to sweep the scattering plane

### calcIdealAngles
- input: hkl vector, UB matrix, B matrix, wavelength, reciprocal lattice params     
- output: 2theta, theta, omega, chi, phi     
- use: constrains omega == 0        

### calcIdealAngles2
- input: hkl vector, chi, phi, UB matrix, wavelength, reciprocal lattice params
- output: 2theta, theta, omega
- use: once you have your UB matrix and scattering plane sorted out, it tells you what omega the desired hkl should be at       

### calcIdealAngles3
- input: hkl vector, phi, UB matrix, wavelength, reciprocal lattice params
- output: 2theta, theta, omega, chi         
- use: constrains phi to a fixed value    

## Limits for Eulerian cradle
-30 < eom < 40
-30 < echi < 90
0 < ephi < 359




