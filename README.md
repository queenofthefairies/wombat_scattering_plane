# Scripts for orienting crystal to a specific scattering plane for Wombat
This utilises the modules `OpenOpt` and `UBmatrix` (courtesy of William Ratcliff)    

## UBmatrix functions
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

## wombat_scattering_plane functions

### generate_hkl_list
- input: hkl_limits
- output: hkl_list         
- use: generates a list of hkl reflections within [+- hkl_limit, +- hkl_limit, +- hkl_limit]    

### hkl_allowed
- input: hkl_list, space group number
- output: subset of hkl_list that are allowed by space group        
- use: not all space groups accommodated yet!    

### is_angle_accessible
- input: twotheta, omega, chi, phi, wom_stth
- output: boolean (1 = accessible)         
- use: given limits on Eulerian cradle and range of Wombat detector, determine whether reciprocal space location is accessible or not    

### accessible_hkl_omega_zero_list
- input: sample_name_prefix, hkl_list_to_test, UB_matrix, B_matrix, wavelength, reciprocal lattice params, wom_stth
- output: dataframe of accessible hkl with omega = 0 degrees         
- use: setting omega = 0 degrees, determine which hkl of hkl_list_to_test are accessible on Wombat, save these hkl to a dataframe    

### generate_hkl_eom_scan_script
- input: accessible_hkl_df, sample_name_prefix, eom_min, eom_step, num_steps, oscillations_per_step
- output: script file as .txt         
- use: generates part of a script that will do eom scans (radcollscans) of a list of hkl reflections    

### evaluate_possible_scattering_planes
- input: sample_name_prefix, UB_matrix, wavelength, reciprocal lattice params
- output: dataframe and .xlsx of accessible scattering planes and the echi, ephi angles that define them         
- use: goes through common scattering planes (e.g. hk0, h0l, hhl etc.) and evaluates whether these are accessible given the UB matrix and Eulerian cradle limits    

### accessible_hkl_in_scattering_plane
- input: sample_name_prefix, plane_name, hkl1, hkl2, hkl_max_component_val, UB_matrix, wavelength, reciprocal lattice params, wom_stth
- output: 2theta, theta, omega, chi         
- use: For a scattering plane defined by two reciprocal space vectors hkl1 and hkl2, make a list of all in-plane reflections accessible on Wombat given eom and twotheta limits. the in-plane reflections to test go up to a maximum Q of +- Y x hkl1 +- Y hkl2, where Y is hkl_max_component_val i.e. if plane is hk0 and Y = 4, then reflections from (-4, -4, 0) to (4, 4, 0) will be tested    

### hkl_in_plane_omega
- input: sample_name_prefix, plane_name, hkl1, hkl2, hkl_to_find, UB_matrix, wavelength, reciprocal lattice params, wom_stth
- output: omega of hkl_to_find if hkl_to_find is in the scattering plane defined by hkl1 and hkl2         
- use: in the scattering plane, find the omega that will give you hkl_to_find    

## Limits for Eulerian cradle
-30 < eom < 40      
-30 < echi < 90         
0 < ephi < 359          




