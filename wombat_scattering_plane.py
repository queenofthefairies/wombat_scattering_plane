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

def hkl_allowed(hkl_list_to_test, space_group_number):
    hkl_allowed_list = []
    for hkl in hkl_list_to_test:
        h = hkl[0]
        k = hkl[1]
        l = hkl[2]

        if space_group_number == 15:
            if h == 0 or k == 0 or l == 0:
                if h == 0:
                    if k != 0:
                        if k%2 == 0:
                            # 0kl: k = 2n
                            hkl_allowed_list.append(hkl)
                    elif k == 0:
                        if l%2 == 0:
                            # 00l: l = 2n
                            hkl_allowed_list.append(hkl)
                    if l == 0:
                        if k%2 == 0:
                            # 0k0: k = 2n
                            hkl_allowed_list.append(hkl)


            else:
                if (h + k)%2 == 0:
                    hkl_allowed_list.append(hkl)
                
    else: 
        print('that space group is not in our list, we should add it!')
        print('hkl allowed list will be empty')
    return hkl_allowed_list



def is_angle_accessible_omega_zero(hkl, twotheta, theta, omega, chi, phi, wom_stth):
    angle_accessible_omega_zero_bool = 0
    # calculate ideal angles (with omega = 0)
    wom_twotheta_min = wom_stth + 2
    wom_twotheta_max = wom_stth + 118
    wom_echi_min = -30
    wom_echi_max = 92.5
    wom_ephi_min = 0
    wom_ephi_max = 395
    if twotheta < wom_twotheta_max and twotheta > wom_twotheta_min:
        if chi < wom_echi_max and chi > wom_echi_min:
            if phi < wom_ephi_max and phi > wom_ephi_min:
                angle_accessible_omega_zero_bool = 1
    return angle_accessible_omega_zero_bool

def accessible_angle_omega_zero_list(sample_name_prefix, hkl_list_to_test, UB_matrix, B_matrix, wavelength, star, wom_stth):
    # go through list of reflections, add those that are accessible to a new list accessible_hkl_list
    accessible_hkl_list = []
    for hkl in hkl_list_to_test:
        twotheta, theta, omega, chi, phi = ubmatrix.calcIdealAngles(hkl, UB_matrix, B_matrix, wavelength, star)
        is_accessible_bool = is_angle_accessible_omega_zero(hkl, twotheta, theta, omega, chi, phi, wom_stth)
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
    accessible_hkl_df.to_csv('{0}_accessible_hkl.csv'.format(sample_name_prefix), index=False)
    accessible_hkl_df.to_excel('{0}_accessible_hkl.xlsx'.format(sample_name_prefix), index=False)
    print('\naccessible hkl saved to {0}_accessible_hkl.csv and {0}_accessible_hkl.xlsx \n'.format(sample_name_prefix))
    
    return accessible_hkl_df

def generate_hkl_eom_scan_script(accessible_hkl_df, sample_name_prefix, eom_min, eom_step, num_steps, oscillations_per_step):
    # sort accessible hkl by echi value 
    accessible_hkl_df = accessible_hkl_df.sort_values(by=['echi'])

    # write script
    # for each hkl, change sample name, drive to appropriate echi and ephi value, do an eom scan 
    script_string = ''
    radcollscan_string = 'RadCollScan eom {0} {1} {2} {3} \n'.format(eom_min, eom_step, num_steps, oscillations_per_step)

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
