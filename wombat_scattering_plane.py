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

#import oo as openopt
#from openopt import NLSP
from oo import NLSP
import ubmatrix

################################################################################
def generate_hkl_list(hkl_limits):
    # create list of reflections to test
    h_list = range(hkl_limits[0],hkl_limits[1])
    k_list = range(hkl_limits[2],hkl_limits[3])
    l_list = range(hkl_limits[4],hkl_limits[5])

    hkl_list_to_test = []
    for h in h_list:
        for k in k_list:
            for l in l_list:
                hkl_to_add = [h,k,l]
                hkl_list_to_test.append(hkl_to_add)
    
    return hkl_list_to_test

################################################################################
def hkl_allowed(hkl_list_to_test, space_group_number):
    hkl_allowed_list = []
    for hkl in hkl_list_to_test:
        h = hkl[0]
        k = hkl[1]
        l = hkl[2]

        if space_group_number == 15:
            # hkl: h+k = 2n
            if h != 0 and k !=0 and l != 0:
                if (h+k)%2 == 0:
                    hkl_allowed_list.append(hkl)
            # h0l: h, l = 2n
            if h != 0 and k == 0 and l != 0:
                if h%2 == 0 and l%2 == 0:
                    hkl_allowed_list.append(hkl)
            # 0kl: k = 2n
            if h == 0 and k != 0 and l != 0:
                if k%2 == 0:
                    hkl_allowed_list.append(hkl)
            # hk0: h + k = 2n
            if h != 0 and k != 0 and l == 0:
                if (h+k)%2 == 0:
                    hkl_allowed_list.append(hkl)
            # h00: h = 2n
            if h != 0 and k == 0 and l == 0:
                if h%2 == 0:
                    hkl_allowed_list.append(hkl)
            # 0k0: k = 2n
            if h == 0 and k != 0 and l == 0:
                if k%2 == 0:
                    hkl_allowed_list.append(hkl)
            # 00l: l = 2n
            if h == 0 and k == 0 and l != 0:
                if l%2 == 0:
                    hkl_allowed_list.append(hkl)

    else: 
        print('that space group is not in our list, we should add it!')
        print('hkl allowed list will be empty')
    return hkl_allowed_list

################################################################################
def is_angle_accessible(twotheta, omega, chi, phi, wom_stth):
    angle_accessible_bool = False
    twotheta_accessible_bool = False
    omega_accessible_bool = False
    chi_accessible_bool = False
    phi_accessible_bool = False
    # calculate ideal angles (with omega = 0)
    wom_twotheta_min = wom_stth + 2
    wom_twotheta_max = wom_stth + 118
    wom_eom_min = -30
    wom_eom_max = 40
    wom_echi_min = -30
    wom_echi_max = 92.5
    wom_ephi_min = 0
    wom_ephi_max = 359.999
    if twotheta < wom_twotheta_max and twotheta > wom_twotheta_min:
        twotheta_accessible_bool = True

    if omega < wom_eom_max and omega > wom_eom_min:
        omega_accessible_bool = True

    if chi < wom_echi_max and chi > wom_echi_min:
        chi_accessible_bool = True
    
    if phi < wom_ephi_max and phi > wom_ephi_min:
        phi_accessible_bool = 1

    if twotheta_accessible_bool and omega_accessible_bool and chi_accessible_bool and phi_accessible_bool:
        angle_accessible_bool = True

    return angle_accessible_bool

################################################################################
def accessible_hkl_omega_zero_list(sample_name_prefix, hkl_list_to_test, UB_matrix, B_matrix, wavelength, star, wom_stth):
    # go through list of reflections, add those that are accessible to a new list accessible_hkl_list
    accessible_hkl_list = []
    for hkl in hkl_list_to_test:
        twotheta, theta, omega, chi, phi = ubmatrix.calcIdealAngles(hkl, UB_matrix, B_matrix, wavelength, star)
        is_accessible_bool = is_angle_accessible(twotheta, omega, chi, phi, wom_stth)
        if is_accessible_bool == 1:
            accessible_hkl_to_add = [hkl[0], hkl[1], hkl[2], twotheta, omega, chi, phi]
            accessible_hkl_list.append(accessible_hkl_to_add)
        else:
            pass
    
    # create a dataframe for accessible hkl and save dataframe to .csv file and .xlsx file
    accessible_hkl_array = np.array(accessible_hkl_list)
    accessible_hkl_df = pd.DataFrame(accessible_hkl_array, columns=['h','k','l','twotheta','eom','echi','ephi'])
    print('\ngiven wombat_stth = {0} degrees and wavelength = {1} Angstrom, there are {2} accessible reflections with eom = 0 degrees \n'.format(wom_stth, wavelength, len(accessible_hkl_list)))
    accessible_hkl_df = accessible_hkl_df.sort_values(by=['twotheta'])
    accessible_hkl_df.to_csv('{0}_hkl_accessible.csv'.format(sample_name_prefix), index=False)
    accessible_hkl_df.to_excel('{0}_hkl_accessible.xlsx'.format(sample_name_prefix), index=False)
    print('\naccessible hkl saved to {0}_hkl_accessible.csv and {0}_hkl_accessible.xlsx \n'.format(sample_name_prefix))
    
    return accessible_hkl_df

################################################################################
def generate_hkl_eom_scan_script(accessible_hkl_df, sample_name_prefix, eom_min, eom_step, num_steps, oscillations_per_step):
    ''' generates part of a script that will do eom scans (radcollscans) of a list of hkl reflections '''
    # sort accessible hkl by echi value 
    accessible_hkl_df = accessible_hkl_df.sort_values(by=['echi'])

    # write script
    # for each hkl, change sample name, drive to appropriate echi and ephi value, do an eom scan 
    script_string = ''
    radcollscan_string = 'radcollscan eom {0} {1} {2} {3} \n'.format(eom_min, eom_step, num_steps, oscillations_per_step)

    for index, hkl_info in accessible_hkl_df.iterrows():
        sample_name_str = 'samplename {0} hkl {1:.2f} {2:.2f} {3:.2f} \n'.format(sample_name_prefix, hkl_info['h'], hkl_info['k'], hkl_info['l'])
        echi_str = 'drive echi {0:.3f} \n'.format(hkl_info['echi'])
        ephi_str = 'drive ephi {0:.3f} \n'.format(hkl_info['ephi'])
        script_string = script_string + sample_name_str + echi_str + ephi_str + radcollscan_string + '\n' 

    # write script to file
    hkl_eom_scan_script_filename = '{0}_hkl_eom_scan_draft_script.txt'.format(sample_name_prefix)
    with open(hkl_eom_scan_script_filename, "w") as file:
        file.write(script_string)
    print()
    print('hkl eom scan draft script saved to {0}'.format(hkl_eom_scan_script_filename))
    print('note you will need to edit the script before running it in gumtree \n')

    return 

################################################################################
def evaluate_possible_scattering_planes(sample_name_prefix, UB_matrix, wavelength, star):
    possible_planes_dict = {'hk0': [[1, 0, 0], [0, 1, 0]], # hk0 family
                            'h-k0': [[1, 0, 0], [0, -1, 0]],
                            '-h-k0': [[-1, 0, 0], [0, -1, 0]],
                            '-hk0': [[-1, 0, 0], [0, 1, 0]],
                            'h0l': [[1, 0, 0], [0, 0, 1]], # h0l family
                            'h0-l': [[1, 0, 0], [0, 0, -1]],
                            '-h0-l': [[-1, 0, 0], [0, 0, -1]],
                            '-h0l': [[-1, 0, 0], [0, 0, 1]],
                            '0kl': [[0, 1, 0], [0, 0, 1]], # 0kl family
                            '0k-l': [[0, 1, 0], [0, 0, -1]],
                            '0-k-l': [[0, -1, 0], [0, 0, -1]],
                            '0-kl': [[0, -1, 0], [0, 0, 1]],
                            'hhl': [[1, 1, 0], [0, 0, 1]], # hhl family
                            '-h-hl': [[-1, -1, 0], [0, 0, 1]],
                            'hh-l': [[1, 1, 0], [0, 0, -1]],
                            '-h-h-l': [[-1, -1, 0], [0, 0, -1]],
                            'hkh': [[1, 0, 1], [0, 1, 0]], # hkh family
                            '-hk-h': [[-1, 0, -1], [0, 1, 0]],
                            'h-kh': [[1, 0, 1], [0, -1, 0]],
                            '-h-k-h': [[-1, 0, -1], [0, -1, 0]],
                            'hkk': [[1, 0, 0], [0, 1, 1]], # hkk family
                            '-hkk': [[-1, 0, 0], [0, 1, 1]],
                            'h-k-k': [[1, 0, 0], [0, -1, -1]],
                            '-h-k-k': [[-1, 0, 0], [0, -1, -1]],
                            'h-hl': [[1, -1, 0], [0, 0, 1]], # h-hl family
                            '-hhl': [[-1, 1, 0], [0, 0, 1]],
                            'h-h-l': [[1, -1, 0], [0, 0, -1]],
                            '-hh-l': [[-1, 1, 0], [0, 0, -1]],
                            'hk-h': [[1, 0, -1], [0, 1, 0]], # hk-h family
                            '-hkh': [[-1, 0, 1], [0, 1, 0]],
                            'h-k-h': [[1, 0, -1], [0, -1, 0]],
                            '-h-kh': [[-1, 0, 1], [0, -1, 0]],
                            'hk-k': [[1, 0, 0], [0, 1, -1]], # hk-k family
                            '-hk-k': [[-1, 0, 0], [0, 1, -1]],
                            'h-kk': [[1, 0, 0], [0, -1, 1]],
                            '-h-kk': [[-1, 0, 0], [0, -1, 1]]
                            }
    
    possible_planes_list = []

    # euler cradle limits
    wom_echi_min = -30
    wom_echi_max = 92.5
    wom_ephi_min = 0
    wom_ephi_max = 359.999

    # go through each plane and calculate the eulerian cradle angles required to 
    # put sample into that scattering plane
    # then evaluate whether those eulerian cradle angles are actually accessible
    # or not
    for planes, vectors in possible_planes_dict.items():
        echi, ephi = ubmatrix.calcScatteringPlane(vectors[0], vectors[1], UB_matrix, wavelength, star)
        wom_stth = 13 # dummy value
        wom_twotheta = 90 # dummy value
        eom = 0 
        if echi < wom_echi_min or echi > wom_echi_max:
            echi = echi%360
        if ephi < wom_ephi_min or ephi > wom_ephi_max:
            ephi = ephi%360
        angle_accessible_bool = is_angle_accessible(wom_twotheta, eom, echi, ephi, wom_stth)
        plane_info_to_add = [planes, *vectors[0], *vectors[1], echi, ephi, angle_accessible_bool]
        possible_planes_list.append(plane_info_to_add)
            
        if angle_accessible_bool:
            print('\ndone {0} plane calculation: accessible \n'.format(planes))
        else:
            print('\ndone {0} plane calculation: not accessible \n'.format(planes))
    
    # save scattering plane info to excel spreadsheet
    scattering_planes_df = pd.DataFrame(possible_planes_list, 
                                        columns=['Plane name', 'h1','k1','l1', 
                                                 'h2','k2','l2','echi','ephi', 
                                                 'Accessible?'])
    scattering_planes_df.to_excel('{0}_scattering_planes.xlsx'.format(sample_name_prefix), index=False)
    print('\nScattering plane info saved to {0}_scattering_planes.xlsx  \n'.format(sample_name_prefix))

    # save accessible scattering planes to another spreadsheet
    accessible_scattering_planes_df = scattering_planes_df[scattering_planes_df['Accessible?'] == True]
    accessible_scattering_planes_df.to_excel('{0}_scattering_planes_accessible.xlsx'.format(sample_name_prefix), index=False)
    print('\nAccessible planes also saved to {0}_scattering_planes_accessible.xlsx  \n'.format(sample_name_prefix))

    return

################################################################################
def where_in_the_omega_is_my_hkl(my_hkl, ref_hkl, ref_hkl_omega, wavelength, star):
    ''' Assuming you are in the scattering plane of my_hkl and ref_hkl, and ref_hkl is at ref_omega, 
    find the omega for my_hkl
    '''
    my_hkl_2theta = ubmatrix.calcTwoTheta(my_hkl, star, wavelength)
    ref_hkl_2theta = ubmatrix.calcTwoTheta(ref_hkl, star, wavelength)

    diff_2theta = my_hkl_2theta - ref_hkl_2theta
    if diff_2theta < 0:
        sign_diff_2theta = -1
    else:
        sign_diff_2theta = 1
    
    interplanar_angle = np.arccos(ubmatrix.scalar(*my_hkl, *ref_hkl, star)/
                                  (ubmatrix.modvec(*my_hkl, star)*ubmatrix.modvec(*ref_hkl, star)))*180/np.pi

    my_hkl_omega = -sign_diff_2theta*interplanar_angle + 0.5*(diff_2theta) + ref_hkl_omega

    print('given ({0} {1} {2}) at omega = {3:.2f}, hkl ({4} {5} {6}) is at omega = {7:.2f}'.format(*ref_hkl, ref_hkl_omega,
                                                                                                   *my_hkl, my_hkl_omega))
        
    return my_hkl_omega

################################################################################
def UB_matrix_in_plane(UB_matrix, echi, ephi):
    # echi, ephi in degrees

    chi_angle = echi*np.pi/180 # convert echi to radians
    # right-handed chi rotation matrix
    chi_rot_matrix = np.array([[np.cos(chi_angle), 0, np.sin(chi_angle)],
                                [0, 1, 0],
                                [-np.sin(chi_angle), 0, np.cos(chi_angle)]])


    phi_angle = ephi*np.pi/180 # convert ephi to radians
    # right-handed phi rotation matrix
    phi_rot_matrix = np.array([[np.cos(phi_angle), np.sin(phi_angle), 0],
                                [-np.sin(phi_angle), np.cos(phi_angle), 0],
                                [0, 0, +1]])

    UB_matrix_transformed = np.matmul(chi_rot_matrix, np.matmul(phi_rot_matrix, UB_matrix))
    return UB_matrix_transformed 

################################################################################
def hkl_in_plane_omega_twotheta(hkl_to_find, UB_matrix, echi, ephi, wavelength, star):
    ''' For a scattering plane defined by two angles echi and ephi, 
    find the vector hkl in the plane
    '''
    # apply echi, ephi rotation matrices to transform UB matrix
    UB_matrix_transformed = UB_matrix_in_plane(UB_matrix, echi, ephi)

    # convert hkl to Qx, Qy, Qz
    pos = np.array([*hkl_to_find])
    Qx, Qy, Qz = np.einsum('ij,j...->i...', UB_matrix_transformed,pos)

    # Now calculate theta, omega
    theta_sign = -1
    hkl_two_theta = ubmatrix.calcTwoTheta(hkl_to_find, star, wavelength)

    Qx = np.asarray(Qx)
    Qy = np.asarray(Qy)
    Qz = np.asarray(Qz)

    #QC = np.array([Qx,Qy,Qz*0.0])
    #q = np.linalg.norm(QC)

    U1V = np.array([Qx.flatten(),Qy.flatten(),Qz.flatten()],dtype=float).flatten()

    U1V/=np.linalg.norm(U1V)
    U2V = np.array([0.0,0.0,1.0],dtype=float)

    U1V /=np.linalg.norm(U1V)
    U2V /=np.linalg.norm(U2V)
    T3V = np.cross(U1V.T, U2V.T).T
    
    TV = np.array([U1V,U2V,T3V])
    R = np.linalg.inv(TV)
    
    ss = 1.0
    
    cossgl = np.sqrt(R[0,0]*R[0,0]+R[1,0]*R[1,0])
    om = np.rad2deg(np.arctan2(R[1,0]/cossgl,R[0,0]/cossgl))
    
    hkl_omega = +om - np.sign(theta_sign)*ss*0.5*hkl_two_theta #+ A3Off

    if np.isnan(hkl_omega) and np.isnan(hkl_two_theta):
        pass
    elif np.isnan(hkl_omega) or np.isnan(hkl_two_theta):
        pass
    else:
        print('hkl ({0:.2f} {1:.2f} {2:.2f})    2theta = {3:.2f}, omega = {4:.2f}'.format(*hkl_to_find, 
                                                                                        hkl_two_theta, 
                                                                                        hkl_omega))
    return hkl_omega, hkl_two_theta

################################################################################
def accessible_hkl_in_scattering_plane(sample_name_prefix, plane_name, hkl1, hkl2, hkl_max_component_val, UB_matrix, wavelength, star, wom_stth):
    ''' For a scattering plane defined by two reciprocal space vectors hkl1 and hkl2, make
    a list of all in-plane reflections accessible on Wombat given eom and twotheta limits
    
    the in-plane reflections to test go up to a maximum Q of +- Y x hkl1 +- Y hkl2,
    where Y is hkl_max_component_val
    i.e. if plane is hk0 and Y = 4, then reflections from (-4, -4, 0) to (4, 4, 0)
    will be tested'''
    echi, ephi = ubmatrix.calcScatteringPlane(hkl1, hkl2, UB_matrix, wavelength, star)
    # euler cradle limits
    wom_echi_min = -30
    wom_echi_max = 92.5
    wom_ephi_min = 0
    wom_ephi_max = 359.999
    if echi < wom_echi_min or echi > wom_echi_max:
        echi = echi%360
    if ephi < wom_ephi_min or ephi > wom_ephi_max:
        ephi = ephi%360

    # make a list of in-plane reflections to test up to a certain maximum Q
    # in order to add the vectors defining the scattering plane together, make them arrays rather than lists
    hkl1_array = np.array(hkl1)
    hkl2_array = np.array(hkl2)

    hkl_in_plane_to_test_list = []
    for hkl1_component_val in range(-hkl_max_component_val, hkl_max_component_val+1):
        for hkl2_component_val in range(-hkl_max_component_val, hkl_max_component_val+1):
            hkl_to_test = hkl1_component_val*hkl1_array + hkl2_component_val*hkl2_array
            hkl_in_plane_to_test_list.append(hkl_to_test)

    # test the reflections
    hkl_in_plane_info_list = []
    for hkl in hkl_in_plane_to_test_list:
        hkl_list = hkl.tolist()

        eom, twotheta = hkl_in_plane_omega_twotheta(hkl, UB_matrix, echi, ephi, wavelength, star)

        # euler cradle limits
        wom_eom_min = -179.99999
        wom_eom_max = 179.99999
        if eom < wom_eom_min or echi > wom_eom_max:
            eom = eom%360

        angle_accessible_bool = is_angle_accessible(twotheta, eom, echi, ephi, wom_stth)
        hkl_in_plane_info_list.append([*hkl_list, twotheta, eom, echi, ephi, wom_stth, angle_accessible_bool])

        if angle_accessible_bool:
            print('reflection {0} {1} {2} in {3} plane is accessible'.format(*hkl_list, plane_name))
        else:
            print('reflection {0} {1} {2} in {3} plane is not accessible'.format(*hkl_list, plane_name))

    # save scattering plane info to excel spreadsheet
    in_plane_hkl_df = pd.DataFrame(hkl_in_plane_info_list, columns=['h','k','l', 'twotheta', 'eom', 'echi','ephi', 'Wombat stth', 'Accessible?'])
    in_plane_hkl_df = in_plane_hkl_df.sort_values(by=['eom'])
    in_plane_hkl_df.to_excel('{0}_{1}_plane_hkl.xlsx'.format(sample_name_prefix, plane_name), index=False)
    print('\nIn-plane hkl info saved to {0}_{1}_plane_hkl.xlsx  \n'.format(sample_name_prefix, plane_name))

    # save accessible scattering planes to another spreadsheet
    accessible_in_plane_hkl_df = in_plane_hkl_df[in_plane_hkl_df['Accessible?'] == True]
    accessible_in_plane_hkl_df.to_excel('{0}_{1}_plane_hkl_accessible.xlsx'.format(sample_name_prefix, plane_name, index=False))
    print('\nAccessible in-plane hkl also saved to {0}_{1}_plane_hkl_accessible.xlsx  \n'.format(sample_name_prefix, plane_name))

    return