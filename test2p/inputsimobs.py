from constant import year

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
mp_Mearth = 2.9808009550481276
a_AU = 1.3023467300523999
e_orbit = 0.340198255282122
# i_orbit: degree: 0-90 or 0-180 ?
i_orbit = 65.84364682095553
# ascend_node_Omega: phase
ascend_node_Omega =  30.0 # 0-360
# periapsis_omega: phase
periapsis_omega = 200.0  # 0-360
# M0: phase
M0 =  100.0 # 0-360
#
mp_Mearth2 = 24.633853751287443
a_AU2 = 1.80123370931395
#882.9453243126169
e_orbit2 = 0.1483524402330199
i_orbit2 = 105.49297579977537
ascend_node_Omega2 =  100
periapsis_omega2 = 20
M02 =  200.0

