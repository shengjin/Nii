from constant import year
from math import radians,cos

one_planet = True #False


#target star:
ra0 = 5.0
dec0 = 14.0

#field of view:
w_ra_dec = 0.44

# number reference stars:
n_ref_stars = 8
# pm, mu as/year
ref_stars_pm_ra = [0.000, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.00] 
# pm, mu as/year
ref_stars_pm_dec = [0.000, 0.0, 0.00, 0.00, -0.00, -0.0, -0.0, -0.00] 
#ref_stars_pm_dec = [0.001, 0.1, 0.01, 0.02, -0.02, -0.01, -0.1, -0.001] 


# observation period ( time in days )
t0 = 0
t1 = year*5
# observation frequency
N_time_init = 300.0
drop_ratio = 0.3 # drop_ratio*100%
drop_seg_min = 4
drop_seg_max = 8

# noise simulation
noise_mean = 0.0
noise_std = 1.0



# proxima A,B
# Proxima C
#mp_Mearth = 312.0
#ms_Msun = 0.1221
#d_pc = 1.3012
#a_AU = 5.2

ms_Msun = 1.0
d_pc = 3.0
#
mp_Mearth = 8.1
a_AU = 1.3
e_orbit = 0.1
# i_orbit: degree: 0-90 or 0-180 ?
i_orbit = 30.1
cos_i_orbit = cos(radians(i_orbit))
# ascend_node_Omega: phase
ascend_node_Omega =  30.0 # 0-360
# periapsis_omega: phase
periapsis_omega = 200.0  # 0-360
# M0: phase
M0 =  100.0 # 0-360
#
mp_Mearth2 = 30.0
a_AU2 = 2.0
e_orbit2 = 0.1
i_orbit2 = 20.5
cos_i_orbit2 = cos(radians(i_orbit2))
ascend_node_Omega2 =  100.0 # 0-360
periapsis_omega2 = 20.0  # 0-360
M02 =  200.0 # 0-360
#


