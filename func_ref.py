
# create some random reference in the field
# coordinate: 0,0 at the center of the field of view 
def create_random_ref_star(n_ref_stars, w_ra_dec): #ra0, dec0, w_ra_dec):
    import numpy as np
    rand_2d_array = np.random.rand(n_ref_stars, 2)
    ref_stars = np.zeros((n_ref_stars, 2), dtype=np.float64)
    ref_stars[:,0] = (rand_2d_array[:,0]-0.5)*w_ra_dec # + ra0
    ref_stars[:,1] = (rand_2d_array[:,1]-0.5)*w_ra_dec # + dec0
    return ref_stars


# calculate the dx dy between target and the reference stars
def simulate_dx_dy(ra, dec, ref_stars, ref_stars_pm_ra, ref_stars_pm_dec, t_in_days): #, mu_per_pixel, resolution):
    # ra, dec: coordinate of the target start
    import numpy as np
    from constant import year
    n_ref_stars = ref_stars.shape[0]
    dra_ddec_pm = np.zeros((n_ref_stars, 2), dtype=np.float64)
    #
    # simulate pm of ref_stars
    ref_stars_pm_ra = np.array(ref_stars_pm_ra)
    ref_stars_pm_dec = np.array(ref_stars_pm_dec)
    #
    dra_ddec_pm[:,0] = ref_stars[:,0] - ra +  ref_stars_pm_ra*t_in_days/year
    dra_ddec_pm[:,1] = ref_stars[:,1] - dec + ref_stars_pm_dec*t_in_days/year
    #
    return dra_ddec_pm


def calc_dx_dy(delta_dx_dy_pm, i, dra_ddec_old_pm, dra_ddec_new_pm, ref_stars_pm_ra, ref_stars_pm_dec, t_in_days, noise_mean, noise_std, n_ref_stars): 
    import numpy as np
    from constant import year
    #
    ref_stars_pm_ra = np.array(ref_stars_pm_ra)
    ref_stars_pm_dec = np.array(ref_stars_pm_dec)
    ##
    # ra
    noise = np.random.normal(noise_mean,noise_std,n_ref_stars)
    delta_dx_dy_pm[i,:,0] = dra_ddec_new_pm[:,0] - dra_ddec_old_pm[:,0] + noise - ref_stars_pm_ra*t_in_days/year
    ##
    # dec
    noise = np.random.normal(noise_mean,noise_std,n_ref_stars)
    delta_dx_dy_pm[i,:,1] = dra_ddec_new_pm[:,1] - dra_ddec_old_pm[:,1] + noise - ref_stars_pm_dec*t_in_days/year
    #
    return delta_dx_dy_pm 



def create_delta_dxdy(time_con, n_ref_stars, w_ra_dec, as_mu, noise_mean, noise_std, ref_stars_pm_ra, ref_stars_pm_dec):
    import numpy as np
    from constant import degree1_mu
    ##
    N_time = len(time_con)
    ##
    ref_stars_coor_in_mu = create_random_ref_star(n_ref_stars, w_ra_dec)*degree1_mu
    ##
    delta_dx_dy_pm = np.zeros((int(N_time)-1, n_ref_stars, 2), dtype=np.float64)
    ##
    dra_ddec_old_pm = simulate_dx_dy(as_mu[0,0], as_mu[0,1], ref_stars_coor_in_mu, ref_stars_pm_ra, ref_stars_pm_dec, time_con[0])
    #
    for i in range(int(N_time)-1):
        #
        dra_ddec_new_pm = simulate_dx_dy(as_mu[i+1,0], as_mu[i+1,1], ref_stars_coor_in_mu, ref_stars_pm_ra, ref_stars_pm_dec, time_con[i+1])
        #
        calc_dx_dy(delta_dx_dy_pm, i, dra_ddec_old_pm, dra_ddec_new_pm, ref_stars_pm_ra, ref_stars_pm_dec, time_con[i+1], noise_mean, noise_std, n_ref_stars)

    return delta_dx_dy_pm


def create_delta_dxdy_pureAsMu(time_con,  as_mu):
    import numpy as np
    from constant import degree1_mu
    ##
    N_time = len(time_con)
    ##
    ##
    delta_dx_dy = np.zeros((int(N_time)-1, 2), dtype=np.float64)
    ##
    for i in range(int(N_time)-1):
        #
        # (ref_i-as_i) - (ref_0-as_0) = (ref_i-as_i-ref_0+as_0)
        # NOTE: i+1 
        delta_dx_dy[i,0] = as_mu[0,0] - as_mu[i+1,0]
        delta_dx_dy[i,1] = as_mu[0,1] - as_mu[i+1,1]
        
    return delta_dx_dy



