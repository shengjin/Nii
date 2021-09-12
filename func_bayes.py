from inputbayes import *
#from inputsimobs import one_planet
import copy

########################################
# prior: eccentricity
def prior_e(ecc):
    if(ecc<ecc_max) and (ecc>=ecc_min):
        return 1/(ecc_max-ecc_min)
    else:
        return 0

########################################
# prior: inclination in degree
def prior_cos_i(cos_incl):
    if(cos_incl<cos_incl_max) and (cos_incl>=cos_incl_min):
        return 1/(cos_incl_max-cos_incl_min)
    else:
        return 0

########################################
# prior: ascend node Omega in degree
def prior_anO(an_Omega):
    if(an_Omega<an_Omega_max) and (an_Omega>=an_Omega_min):
        return 1/(an_Omega_max-an_Omega_min)
    else:
        return 0

########################################
# prior: periapsis omega in degree
def prior_po(p_omega):
    if(p_omega<p_omega_max) and (p_omega>=p_omega_min):
        return 1/(p_omega_max-p_omega_min)
    else:
        return 0

########################################
# prior: anomaly M0 in degree
def prior_M0(M0):
    if(M0<M0_max) and (M0>=M0_min):
        return 1/(M0_max-M0_min)
    else:
        return 0

########################################
# prior: orbital period
def prior_period(period):
    from  math import log
    if(period<period_max) and (period>=period_min):
        return 1/(period*log(period_max/period_min))
    else:
        return 0

########################################
# prior: planet mass in Earth mass
def prior_mass(p_mass):
    from  math import log
    if(p_mass<p_mass_max) and (p_mass>=p_mass_min):
        return 1/(p_mass*log(p_mass_max/p_mass_min))
    else:
        return 0

########################################
# prior: unknown errors
def prior_var_uke(var_uk_err):
    from  math import log
    if(var_uk_err<var_uk_err_max) and (var_uk_err>=var_uk_err_min):
        return 1/log(var_uk_err_a+var_uk_err_max/var_uk_err_a)/(var_uk_err+var_uk_err_a)
    else:
        return 0

# NOTE order change for the last 3 varibale
def para_boundary(x):
    from numpy import cos
    if one_planet_fit:
        cos_incl = x[0]
        ecc = x[1]
        anO = x[2]
        po = x[3]
        M0 = x[4]
        period = x[7]
        mp = x[5]
        var_uke = x[6]
        if (cos_incl<cos_incl_max) and (cos_incl>cos_incl_min) and (ecc<ecc_max) and (ecc>ecc_min) and (anO<an_Omega_max) and (anO>an_Omega_min) and (po<p_omega_max) and (po>p_omega_min) and (M0<M0_max) and (M0>M0_min) and (mp<p_mass_max) and (mp>p_mass_min) and (var_uke<var_uk_err_max) and (var_uke>var_uk_err_min) and (period<period_max) and (period>period_min):
            return True
        else:
            return False
    else:
        cos_incl = x[0]
        ecc = x[1]
        anO = x[2]
        po = x[3]
        M0 = x[4]
        period = x[7]
        mp = x[5]
        var_uke = x[6]
        cos_incl2 = x[8]
        ecc2 = x[9]
        anO2 = x[10]
        po2 = x[11]
        M02 = x[12]
        period2 = x[14]
        mp2 = x[13]
        if (cos_incl<cos_incl_max) and (cos_incl>cos_incl_min) and (ecc<ecc_max) and (ecc>ecc_min) and (anO<an_Omega_max) and (anO>an_Omega_min) and (po<p_omega_max) and (po>p_omega_min) and (M0<M0_max) and (M0>M0_min) and (mp<p_mass_max) and (mp>p_mass_min) and (var_uke<var_uk_err_max) and (var_uke>var_uk_err_min) and (period<period_max) and (period>period_min) and (cos_incl2<cos_incl_max) and (cos_incl2>cos_incl_min) and (ecc2<ecc_max) and (ecc2>ecc_min) and (anO2<an_Omega_max) and (anO2>an_Omega_min) and (po2<p_omega_max) and (po2>p_omega_min) and (M02<M0_max) and (M02>M0_min) and (mp2<p_mass_max) and (mp2>p_mass_min) and (period2<period_max) and (period2>period_min):
            return True
        else:
            return False


# NOTE order change for the last 3 varibale
def log_prior(x):
    from numpy import log, cos
    # Gregory 2005
    if one_planet_fit:
        cos_incl = x[0]
        ecc = x[1]
        anO = x[2]
        po = x[3]
        M0 = x[4]
        period = x[7]
        mp = x[5]
        var_uke = x[6]
        l_p = log(prior_cos_i(cos_incl))+log(prior_e(ecc))+log(prior_anO(anO))+log(prior_po(po))+log(prior_M0(M0))+log(prior_period(period))+log(prior_mass(mp))+log(prior_var_uke(var_uke))
    else:
        cos_incl = x[0]
        ecc = x[1]
        anO = x[2]
        po = x[3]
        M0 = x[4]
        period = x[7]
        mp = x[5]
        var_uke = x[6]
        cos_incl2 = x[8]
        ecc2 = x[9]
        anO2 = x[10]
        po2 = x[11]
        M02 = x[12]
        period2 = x[14]
        mp2 = x[13]
        l_p = log(prior_cos_i(cos_incl))+log(prior_e(ecc))+log(prior_anO(anO))+log(prior_po(po))+log(prior_M0(M0))+log(prior_period(period))+log(prior_mass(mp))+log(prior_var_uke(var_uke)) + log(prior_cos_i(cos_incl2))+log(prior_e(ecc2))+log(prior_anO(anO2))+log(prior_po(po2))+log(prior_M0(M02))+log(prior_period(period2))+log(prior_mass(mp2))
    return l_p

def log_likelihood(delta_dx_dy_sig, var_uke_init, delta_dx_dy_obs, delta_dx_dy_mc):
    # Gregory 2005
    import numpy as np
    N = len(delta_dx_dy_sig)
    llkl_ra = 0
    for i in range(N):
        sig_power = (delta_dx_dy_sig[i,0]**2.0+var_uke_init**2.0)
        AC_one = np.log((2*np.pi)**(-0.5)) + np.log(sig_power**(-0.5))
        exp_one = -(delta_dx_dy_mc[i,0]-delta_dx_dy_obs[i,0])**2.0/2/sig_power
        llkl_ra = llkl_ra+AC_one+exp_one
    llkl_dec = 0
    for i in range(N):
        sig_power = (delta_dx_dy_sig[i,1]**2.0+var_uke_init**2.0)
        AC_one = np.log((2*np.pi)**(-0.5)) + np.log(sig_power**(-0.5))
        exp_one = -(delta_dx_dy_mc[i,1]-delta_dx_dy_obs[i,1])**2.0/2/sig_power
        llkl_dec = llkl_dec+AC_one+exp_one
    llkl = llkl_ra + llkl_dec
    return llkl


def init_para():
    import numpy as np
    if one_planet_fit:
        cos_incl_init = np.random.rand()*(cos_incl_max-cos_incl_min)+cos_incl_min
        ecc_init = np.random.rand()*(ecc_max-ecc_min)+ecc_min
        anO_init = np.random.rand()*(an_Omega_max-an_Omega_min)+an_Omega_min
        po_init = np.random.rand()*(p_omega_max-p_omega_min)+p_omega_min
        M0_init = np.random.rand()*(M0_max-M0_min)+M0_min
        mp_init = np.random.rand()*(p_mass_max-p_mass_min)+p_mass_min
        var_uke_init = np.random.rand()*(var_uk_err_max-var_uk_err_min)+var_uk_err_min
        period_init = np.random.rand()*(period_max-period_min)+period_min
        x_init = np.zeros(n_dim, dtype=np.float64)
        x_init[0] = cos_incl_init 
        x_init[1] = ecc_init 
        x_init[2] = anO_init 
        x_init[3] = po_init 
        x_init[4] = M0_init 
        x_init[5] = mp_init 
        x_init[6] = var_uke_init 
        x_init[7] = period_init 
    else:
        cos_incl_init = np.random.rand()*(cos_incl_max-cos_incl_min)+cos_incl_min
        ecc_init = np.random.rand()*(ecc_max-ecc_min)+ecc_min
        anO_init = np.random.rand()*(an_Omega_max-an_Omega_min)+an_Omega_min
        po_init = np.random.rand()*(p_omega_max-p_omega_min)+p_omega_min
        M0_init = np.random.rand()*(M0_max-M0_min)+M0_min
        mp_init = np.random.rand()*(p_mass_max-p_mass_min)+p_mass_min
        var_uke_init = np.random.rand()*(var_uk_err_max-var_uk_err_min)+var_uk_err_min
        period_init = np.random.rand()*(period_max-period_min)+period_min
        #
        cos_incl_init2 = np.random.rand()*(cos_incl_max-cos_incl_min)+cos_incl_min
        ecc_init2 = np.random.rand()*(ecc_max-ecc_min)+ecc_min
        anO_init2 = np.random.rand()*(an_Omega_max-an_Omega_min)+an_Omega_min
        po_init2 = np.random.rand()*(p_omega_max-p_omega_min)+p_omega_min
        M0_init2 = np.random.rand()*(M0_max-M0_min)+M0_min
        mp_init2 = np.random.rand()*(p_mass_max-p_mass_min)+p_mass_min
        period_init2 = np.random.rand()*(period_max-period_min)+period_min
        x_init = np.zeros(n_dim, dtype=np.float64)
        x_init[0] = cos_incl_init 
        x_init[1] = ecc_init 
        x_init[2] = anO_init 
        x_init[3] = po_init 
        x_init[4] = M0_init 
        x_init[5] = mp_init 
        x_init[6] = var_uke_init 
        x_init[7] = period_init 
        x_init[8] = cos_incl_init2 
        x_init[9] = ecc_init2 
        x_init[10] = anO_init2 
        x_init[11] = po_init2 
        x_init[12] = M0_init2
        x_init[13] = mp_init2 
        x_init[14] = period_init2 
    return x_init

def init_para_fix():
    import numpy as np
    cos_incl_init = 0.01
    ecc_init = 0.001
    anO_init = 1.0 
    po_init = 2.2
    M0_init = 0.1
    mp_init = 8.1
    var_uke_init = 0
    period_init = 365
    x_init = np.zeros(n_dim, dtype=np.float64)
    x_init[0] = cos_incl_init 
    x_init[1] = ecc_init 
    x_init[2] = anO_init 
    x_init[3] = po_init 
    x_init[4] = M0_init 
    x_init[5] = mp_init 
    x_init[6] = var_uke_init 
    x_init[7] = period_init 
    return x_init


######################
def gaussian_proposal(x, sigma=0.1):
    """
    Gaussian proposal distribution.
    :param x: Parameter array
    :param sigma:
        Standard deviation of Gaussian distribution. Can be scalar
        or vector of length(x)
    #
    :returns: (new parameters, ratio of proposal densities)
    #
    This proposal is nearly as simple as it gets, mathematically it is:
        q(x∗∣xi)=Normal(xi,σ2),
    that is, a Gaussian centered on the current position xi
    with variance given by a standard deviation parameter σ.
    """
    import numpy as np
    # Draw x_star
    N = len(x)
    dx = np.random.randn(N) 
    if one_planet_fit:
        dx[0] = dx[0]*sigma[0]
        dx[1] = dx[1]*sigma[1]
        dx[2] = dx[2]*sigma[2]
        dx[3] = dx[3]*sigma[3]
        dx[4] = dx[4]*sigma[4]
        dx[5] = dx[5]*sigma[5]
        dx[6] = dx[6]*sigma[6]
        dx[7] = dx[7]*sigma[7]
    else:
        dx[0] = dx[0]*sigma[0]
        dx[1] = dx[1]*sigma[1]
        dx[2] = dx[2]*sigma[2]
        dx[3] = dx[3]*sigma[3]
        dx[4] = dx[4]*sigma[4]
        dx[5] = dx[5]*sigma[5]
        dx[6] = dx[6]*sigma[6]
        dx[7] = dx[7]*sigma[7]
        dx[8] = dx[8]*sigma[8]
        dx[9] = dx[9]*sigma[9]
        dx[10] = dx[10]*sigma[10]
        dx[11] = dx[11]*sigma[11]
        dx[12] = dx[12]*sigma[12]
        dx[13] = dx[13]*sigma[13]
        dx[14] = dx[14]*sigma[14]
    x_star = x + dx
    # proposal ratio factor is 1 since jump is symmetric
    qxx = 1
    return (x_star, qxx)


def gen_gaus_prop_sig():
    #gprop_sigma = [deg_180,0.012,deg_360, deg_360, deg_360, mp_Mearth, vuk_muas p_days, ]
    #(incl, ecc, anO, po, M0, mp, var_uke,period)
    import numpy as np
    gp_sig = np.zeros(n_dim, dtype=np.float64)
    if one_planet_fit:
        gp_sig[0]=(cos_incl_max-cos_incl_min)*init_ratio_gp
        gp_sig[1]=(ecc_max-ecc_min)*init_ratio_gp
        gp_sig[2]=(an_Omega_max-an_Omega_min)*init_ratio_gp
        gp_sig[3]=(p_omega_max-p_omega_min)*init_ratio_gp
        gp_sig[4]=(M0_max-M0_min)*init_ratio_gp
        gp_sig[5]=(p_mass_max-p_mass_min)*init_ratio_gp
        gp_sig[6]=(var_uk_err_max-var_uk_err_min)*init_ratio_gp
        gp_sig[7]=(period_max-period_min)*init_ratio_gp
    else:
        gp_sig[0]=(cos_incl_max-cos_incl_min)*init_ratio_gp
        gp_sig[1]=(ecc_max-ecc_min)*init_ratio_gp
        gp_sig[2]=(an_Omega_max-an_Omega_min)*init_ratio_gp
        gp_sig[3]=(p_omega_max-p_omega_min)*init_ratio_gp
        gp_sig[4]=(M0_max-M0_min)*init_ratio_gp
        gp_sig[5]=(p_mass_max-p_mass_min)*init_ratio_gp
        gp_sig[6]=(var_uk_err_max-var_uk_err_min)*init_ratio_gp
        gp_sig[7]=(period_max-period_min)*init_ratio_gp
        gp_sig[8]=(cos_incl_max-cos_incl_min)*init_ratio_gp
        gp_sig[9]=(ecc_max-ecc_min)*init_ratio_gp
        gp_sig[10]=(an_Omega_max-an_Omega_min)*init_ratio_gp
        gp_sig[11]=(p_omega_max-p_omega_min)*init_ratio_gp
        gp_sig[12]=(M0_max-M0_min)*init_ratio_gp
        gp_sig[13]=(p_mass_max-p_mass_min)*init_ratio_gp
        gp_sig[14]=(period_max-period_min)*init_ratio_gp
    return gp_sig





def mh_sampler(x0, time_con, data, ln_likelhd_func, beta, delta_dx_dy_sig, prior_func, prop_fn, iterations, prop_fn_kwargs={}):
    """Simple metropolis hastings sampler.
    :param x0: Initial array of parameters.
    :param prop_fn: Function to perform jumps.
    :param prop_fn_kwargs: Keyword arguments for proposal function
    :param iterations: Number of iterations to run sampler. Default=100000
    #:returns:
    #    (chain, acceptance, lnprob) tuple of parameter chain , acceptance rate
    #    and log-posterior chain.
    """
    from func_orbit import gen_dx_dy_mc, gen_dx_dy_mc_2p
    import numpy as np
    np.random.seed()
    ##################
    # number of dimensions
    ndim = len(x0)
    # initialize chain, acceptance rate and lnprob
    chain = np.zeros((iterations, ndim), dtype=np.float128)
    lnprob = np.zeros((iterations), dtype=np.float128)
    accept_rate = np.zeros((iterations), dtype=np.float128)
    ##################
    # first samples
    chain[0,:] = x0
    # delta dx dy in mu_as
    if one_planet_fit:
        delta_dx_dy_init = gen_dx_dy_mc(x0, time_con)
    else:
        delta_dx_dy_init = gen_dx_dy_mc_2p(x0, time_con)
    # ln_posterior
    llhd = ln_likelhd_func(delta_dx_dy_sig, x0[6], data, delta_dx_dy_init)*beta
    # prior 
    lprior = prior_func(x0)
    lnprob0 = llhd + lprior
    lnprob[0] = lnprob0
    accept_rate[0] = 0
    # start loop
    # 
    ## this is real calc accpet rate 
    naccept = 0
    #
    ii = 0
    while  ii <  (iterations-1):
        #
        # fake accpet rate, accpet status of each iteration 
        #naccept = 0
        #
        # propose
        xs, factor = prop_fn(x0, **prop_fn_kwargs)
        ####
        #if option 4: Just reject any proposals that fall outside the unit square
        if para_boundary(xs):
            if one_planet_fit:
                delta_dx_dy_one = gen_dx_dy_mc(xs, time_con)
            else:
                delta_dx_dy_one = gen_dx_dy_mc_2p(xs, time_con)
            # ln_posterior
            llhd = ln_likelhd_func(delta_dx_dy_sig, xs[6], data, delta_dx_dy_one)*beta
            # prior 
            lprior = prior_func(xs)
            lnprob_s = llhd + lprior
            # MH
            # NTOE not compute hastings ratio in the 1st_criterion 
            #       to avoid overflow of exp func
            if (lnprob_s - lnprob0) > np.log(1/factor):
                x0 = xs
                lnprob0 = lnprob_s
                #
                ## real accpet rate
                naccept += 1
                # fake rate
                #naccept = 1
                #
            else:
                H = np.exp(lnprob_s - lnprob0) * factor
                # draw random uniform number
                u = np.random.uniform(0, 1)
                # accept/reject step (update acceptance counter)
                if u < H:
                    x0 = xs
                    lnprob0 = lnprob_s
                    #
                    ## real accpet rate
                    naccept += 1
                    # fake rate
                    #naccept = 1
                    #
            ii += 1
            #
            # save parameter in chain
            # for save only jumped chain
            chain[ii,:] = x0
            # for save all chain
            #chain[ii,:] = xs
            #
            lnprob[ii] = lnprob0
            #
            ## real rate
            #accept_rate[ii] = naccept / ii
            # fake rate
            accept_rate[ii] = naccept 
            #
    return chain, accept_rate, lnprob


def run_mcmc(queue, time_con, i, data, sigma_known, n_iter_once, beta_PT_once, sigma_prop):
    # wrapper
    chains_list = queue.get()
    # 
    # set x_init_once to first chain item
    x_init_once = chains_list[i].chain[0,:] 
    #
    # not work when big memory needed! only work if n_iter_once and n_processes are small
    chains_list[i].chain, chains_list[i].accept_rate, chains_list[i].lnprob = mh_sampler(x_init_once, time_con, data, log_likelihood, beta_PT_once, sigma_known, log_prior, gaussian_proposal, n_iter_once, prop_fn_kwargs={'sigma':sigma_prop})
    #
    queue.put(chains_list)

def gen_two_rand(n_PT):
    import numpy as np
    np.random.seed()
    swap_i = int((n_PT-1)*np.random.rand())
    swap_j = int(swap_i+1)
    return swap_i, swap_j




def mh_crit_beta(x_i, x_j, beta_i, beta_j, time_con, ln_likelhd_func, delta_dx_dy_sig, data):
    #:returns:
    #    (chain, acceptance, lnprob) tuple of parameter chain , acceptance rate
    #    and log-posterior chain.
    from func_orbit import gen_dx_dy_mc, gen_dx_dy_mc_2p
    import numpy as np
    np.random.seed()
    ##################
    # delta dx dy in mu_as
    if one_planet_fit:
        delta_dx_dy_xi = gen_dx_dy_mc(x_i, time_con)
        delta_dx_dy_xj = gen_dx_dy_mc(x_j, time_con)
    else:
        delta_dx_dy_xi = gen_dx_dy_mc_2p(x_i, time_con)
        delta_dx_dy_xj = gen_dx_dy_mc_2p(x_j, time_con)
    # ln_posterior
    ll_xi_bi = ln_likelhd_func(delta_dx_dy_sig, x_i[6], data, delta_dx_dy_xi)*beta_i
    ll_xi_bj = ln_likelhd_func(delta_dx_dy_sig, x_i[6], data, delta_dx_dy_xi)*beta_j
    ll_xj_bi = ln_likelhd_func(delta_dx_dy_sig, x_j[6], data, delta_dx_dy_xj)*beta_i
    ll_xj_bj = ln_likelhd_func(delta_dx_dy_sig, x_j[6], data, delta_dx_dy_xj)*beta_j
    # 
    # default False
    accept = False
    #
    if (ll_xj_bi + ll_xi_bj - ll_xi_bi - ll_xj_bj) > 0 :
        accept = True 
    else: 
        H = np.exp(ll_xj_bi + ll_xi_bj - ll_xi_bi - ll_xj_bj)
        # draw random uniform number
        u = np.random.uniform(0, 1)
        #
        if u < H:
            accept = True
    #
    return accept

def gen_init_sigma_prop(n_PT, init_sigma_prop):
    import numpy as np
    n_dim = len(init_sigma_prop)
    sigma_2d_init = np.zeros((n_PT,n_dim), dtype=float)
    for i in range(n_PT):
        # change j_th parameter
        sigma_2d_init[i,:] = init_sigma_prop
    #for i in range(n_dim):
    #    # set n_PT same parameter for each param(n_dim)
    #    sigma_2d_init[:,i] = init_sigma_prop[i]
    #print(sigma_2d_init)
    return sigma_2d_init

def sigma_tune_boundary(scale_min, scale_max):
    import numpy as np
    sigma_bound = np.zeros((n_dim,2), dtype=float)
    if one_planet_fit:
        sigma_bound[0,0]=(cos_incl_max-cos_incl_min)*scale_min
        sigma_bound[0,1]=(cos_incl_max-cos_incl_min)*scale_max
        sigma_bound[1,0]=(ecc_max-ecc_min)*scale_min
        sigma_bound[1,1]=(ecc_max-ecc_min)*scale_max
        sigma_bound[2,0]=(an_Omega_max-an_Omega_min)*scale_min
        sigma_bound[2,1]=(an_Omega_max-an_Omega_min)*scale_max
        sigma_bound[3,0]=(p_omega_max-p_omega_min)*scale_min
        sigma_bound[3,1]=(p_omega_max-p_omega_min)*scale_max
        sigma_bound[4,0]=(M0_max-M0_min)*scale_min
        sigma_bound[4,1]=(M0_max-M0_min)*scale_max
        sigma_bound[5,0]=(p_mass_max-p_mass_min)*scale_min
        sigma_bound[5,1]=(p_mass_max-p_mass_min)*scale_max
        sigma_bound[6,0]=(var_uk_err_max-var_uk_err_min)*scale_min
        sigma_bound[6,1]=(var_uk_err_max-var_uk_err_min)*scale_max
        sigma_bound[7,0]=(period_max-period_min)*scale_min
        sigma_bound[7,1]=(period_max-period_min)*scale_max
    else:
        sigma_bound[0,0]=(cos_incl_max-cos_incl_min)*scale_min
        sigma_bound[0,1]=(cos_incl_max-cos_incl_min)*scale_max
        sigma_bound[1,0]=(ecc_max-ecc_min)*scale_min
        sigma_bound[1,1]=(ecc_max-ecc_min)*scale_max
        sigma_bound[2,0]=(an_Omega_max-an_Omega_min)*scale_min
        sigma_bound[2,1]=(an_Omega_max-an_Omega_min)*scale_max
        sigma_bound[3,0]=(p_omega_max-p_omega_min)*scale_min
        sigma_bound[3,1]=(p_omega_max-p_omega_min)*scale_max
        sigma_bound[4,0]=(M0_max-M0_min)*scale_min
        sigma_bound[4,1]=(M0_max-M0_min)*scale_max
        sigma_bound[5,0]=(p_mass_max-p_mass_min)*scale_min
        sigma_bound[5,1]=(p_mass_max-p_mass_min)*scale_max
        sigma_bound[6,0]=(var_uk_err_max-var_uk_err_min)*scale_min
        sigma_bound[6,1]=(var_uk_err_max-var_uk_err_min)*scale_max
        sigma_bound[7,0]=(period_max-period_min)*scale_min
        sigma_bound[7,1]=(period_max-period_min)*scale_max
        #
        sigma_bound[8,0]=(cos_incl_max-cos_incl_min)*scale_min
        sigma_bound[8,1]=(cos_incl_max-cos_incl_min)*scale_max
        sigma_bound[9,0]=(ecc_max-ecc_min)*scale_min
        sigma_bound[9,1]=(ecc_max-ecc_min)*scale_max
        sigma_bound[10,0]=(an_Omega_max-an_Omega_min)*scale_min
        sigma_bound[10,1]=(an_Omega_max-an_Omega_min)*scale_max
        sigma_bound[11,0]=(p_omega_max-p_omega_min)*scale_min
        sigma_bound[11,1]=(p_omega_max-p_omega_min)*scale_max
        sigma_bound[12,0]=(M0_max-M0_min)*scale_min
        sigma_bound[12,1]=(M0_max-M0_min)*scale_max
        sigma_bound[13,0]=(p_mass_max-p_mass_min)*scale_min
        sigma_bound[13,1]=(p_mass_max-p_mass_min)*scale_max
        sigma_bound[14,0]=(period_max-period_min)*scale_min
        sigma_bound[14,1]=(period_max-period_min)*scale_max
    return sigma_bound

def sigma_tune_2d_alltune(sigma_prop_slice, tune_arr, sigma_bound):
    # n: n_th parameter to change
    import numpy as np
    n_PT = len(tune_arr)
    n_dim = len(sigma_prop_slice)
    #
    sigma_2d_alltune = np.zeros((n_PT,n_dim), dtype=float)
    #
    for i in range(n_dim):
        # set n_PT same parameter for each param(n_dim)
        sigma_2d_alltune[:,i] = sigma_prop_slice[i]
    for i in range(n_dim):
        # change j_th parameter: n_PT arr for each para(n_dim)
        sigma_2d_alltune[:,i] = tune_arr*sigma_prop_slice[i]
        # out of boundary
        if sigma_2d_alltune[0,i] < sigma_bound[i,0]:
            sigma_2d_alltune[:,i] = np.logspace(sigma_bound[i,0], sigma_bound[i,0]*scale_bound, n_PT)
            sigma_2d_alltune[:,i] = np.logspace(np.log10(sigma_bound[i,0]), np.log10(sigma_bound[i,0]*scale_bound), n_PT)
        if sigma_2d_alltune[-1,i] > sigma_bound[i,1]:
            sigma_2d_alltune[:,i] = np.logspace(np.log10(sigma_bound[i,1]/scale_bound), np.log10(sigma_bound[i,1]), n_PT)
    return  sigma_2d_alltune

def sigma_tune_2d_1tune(sigma_prop_alltune, n, current_sigma_prop_iPT):
    # n: n_th parameter to change
    import numpy as np
    n_PT = sigma_prop_alltune.shape[0]
    n_dim = len(current_sigma_prop_iPT)
    #
    sigma_2d_1tune = np.zeros((n_PT,n_dim), dtype=float)
    #
    for i in range(n_dim):
        # n_PT same parameter
        sigma_2d_1tune[:,i] = current_sigma_prop_iPT[i]
    # change n_th parameter to sigma_prop_alltune
    sigma_2d_1tune[:,n] = sigma_prop_alltune[:,n]
    return  sigma_2d_1tune

def racist(sigma_prop_alltune, ar_2d_alltune):
    import numpy as np
    n_PT, n_dim = ar_2d_alltune.shape[:]
    d_ar_2d = abs(ar_2d_alltune - target_ar)
    sigma_prop_new = np.zeros(n_dim, dtype=float)
    for i in range(n_dim):
        iPT_d_min = d_ar_2d[:,i].argmin()
        #print(iPT_d_min)
        sigma_prop_new[i] = sigma_prop_alltune[iPT_d_min,i]
    #print("d_ar_2d\n")
    #print(d_ar_2d)
    #print("sigma_prop_alltune\n")
    #print(sigma_prop_alltune)
    #print("sigma_prop_new\n")
    #print(sigma_prop_new)
    return sigma_prop_new

def judge_new_sigma(ar_2d_alltune):
    import numpy as np
    n_PT, n_dim = ar_2d_alltune.shape[:]
    #
    print("ar_2d_alltune:\n")
    print(ar_2d_alltune)
    ######### case 1
    # find largest std for n_dim para
    std_arr = np.zeros(n_dim, dtype=float)
    for i in range(n_dim):
        ar_dim = ar_2d_alltune[:,i]
        std_arr[i] = np.std(ar_dim)
    i_change1 = np.argmax(std_arr)
    #print(i_change, "i_change")
    #
    # find closest ar compared to target_ar
    ar_i_ch = ar_2d_alltune[:,i_change1]
    #print(ar_i_ch, "ar_i_ch")
    d_good_ar =  abs(ar_i_ch - target_ar)
    d_min_std = d_good_ar.min()
    i_PT_use1 = np.argmin(d_good_ar)
    #
    ######### case 2 global closest ar
    err_2d_all =  abs(ar_2d_alltune-target_ar)
    i_PT_use2, i_change2 = np.unravel_index(err_2d_all.argmin(), err_2d_all.shape)
    d_min_2d =  err_2d_all.min()
    #
    ######### judgement
    # d_good_ar good
    if (d_min_std < ar_diff_crit):
        # d_min_2d good
        if (d_min_2d < ar_diff_crit):
            np.random.seed()
            num_rand = np.random.rand()
            #print(num_rand)
            if num_rand > 0.5:
                i_change = i_change1
                i_PT_use = i_PT_use1
            else:
                i_change = i_change2
                i_PT_use = i_PT_use2
        else:
            i_change = i_change1
            i_PT_use = i_PT_use1
    # d_good_ar bad, either too low or too high
    else:
        # if d_min_2d good, use case 2
        if (d_min_2d < ar_diff_crit):
            i_change = i_change2
            i_PT_use = i_PT_use2
        # if d_min_2d bad
        else:
            np.random.seed()
            num_rand = np.random.rand()
            if num_rand > 0.5:
                # try , reset to default 10% of all parameters
                i_change = -1
                i_PT_use = 0
            else:
                # e.g. set largest ar in each parameter
                i_change = -1
                i_PT_use = -1
    #
    f=open('judge_new_sigma.dat','a+')
    str_list = "%s%s%s%s%s%s%s%s" % (str(d_min_std)," ",str(d_min_2d)," ",str(i_change)," ",str(i_PT_use),"\n")
    f.write(str_list)
    f.close()
    return  i_change, i_PT_use

def find_bad_chain(ar_arr, ar_min_crit, ar_max_crit):
    import numpy as np
    # now work in case empty array
    #ind_sm = np.where((ar_arr < ar_min_crit))[0]
    #ind_bg = np.where((ar_arr > ar_max_crit))[0]
    #ind = np.concatenate(ind_sm, ind_bg)
    #_, i = np.unique(z, return_index=True)
    #ind[np.sort(i)]
    bad_ind = []
    for i in range(len(ar_arr)):
        if (ar_arr[i] < ar_min_crit) or (ar_arr[i] > ar_max_crit):
            bad_ind.append(i)
    return bad_ind



###  reduce turning times 
def ratio_true(ratio):
    import numpy as np
    randnum = np.random.uniform(0, 1)
    if randnum < ratio:
        return True
    else:
        return False



##########################
## parallel tempering
# philosiphy: n_PT long chains_list with many iterations: simulatanously:
#        simulatanously short(n_iter_swap_av, dn_swap) chains by Queue:
#        in which shuffle last parameter-set and set as x_init* inital para-set
def mcmc_PT(time_con, delta_dx_dy_1refstar, delta_dx_dy_sig):
    import multiprocessing as mp
    import numpy as np
    import class_mcmc
    from inputbayes import n_iter, n_dim, beta_PT, n_iter_swap_av, dn_swap
    # 
    # make print of array elegant
    np.set_printoptions(precision=2)
    #
    # i_tmp, star point of current sub_chain (chains_tmp); 
    i_tmp = int(0)
    # i_nt, star point of next sub_chain (next chains_tmp); 
    i_nt = int(0)
    # number of beta and number chains
    n_PT = len(beta_PT)
    #
    # chains_list : the whole chain
    chains_list = {x: class_mcmc.Chains(n_iter,n_dim) for x in range(n_PT)}
    #
    # first x_init for first chains_tmp
    x_init_tmp = {x: init_para() for x in range(n_PT)}
    # test case: fix value, e.g., para of the smltorbit, to see mcmc behaviour
    #x_init_tmp = {x: init_para_fix() for x in range(n_PT)}
    #
    # set init sigma_prop
    init_sigma_prop = gen_gaus_prop_sig()
    # sigma_prop in 2d. init to same value for each chain 
    sigma_prop = gen_init_sigma_prop(n_PT, init_sigma_prop)
    #
    # to calc real ar
    ar_arr = np.zeros(n_PT)
    ar_arr_tmp = np.zeros(n_PT)
    #
    # mp stuff
    que = mp.Queue()
    #
    # cald sigma_bound for tune limit
    sigma_bound = sigma_tune_boundary(scale_min, scale_max)
    #
    #
    while i_nt < (n_iter):
        #
        print("###################################################################")
        print("###################################################################")
        print("###################################################################")
        print("i_tmp(current starting iteration):",i_tmp)
        #
        # reset random number generator
        np.random.seed()
        # int-length of sub_chain, chain_tmp
        dn_swap_once = np.random.uniform(-dn_swap, dn_swap)
        n_iter_once =  round(n_iter_swap_av+dn_swap_once) 
        # new i_nt 
        i_nt = i_tmp + n_iter_once
        # new i_nt for the final tail
        if i_nt > n_iter:
            n_iter_once = int(n_iter-i_tmp)
            i_nt = i_tmp + n_iter_once
        #
        ####################################
        # parallel stuff
        # 
        # create new chain_tmp
        chains_tmp = {x: class_mcmc.Chains(n_iter_once,n_dim) for x in range(n_PT)}
        # set 1st chians_tmp item to the last item of previous chains_tmp
        for i in range(n_PT):
            chains_tmp[i].chain[0,:]  = x_init_tmp[i] 
        #
        #  put chains_tmp to Queue
        que.put(chains_tmp)
        #
        # Setup a list of processes that we want to run
        processes = [mp.Process(target=run_mcmc, args=(que, time_con, i, delta_dx_dy_1refstar, delta_dx_dy_sig, n_iter_once,  beta_PT[i], sigma_prop[i,:])) for i in range(n_PT)]
        #
        # Run processes
        for p in processes:
            p.start()
        # blocks the calling thread until the thread whose join() method is called is terminated.
        for p in processes:
            p.join()
#
        chains_tmp = que.get(chains_tmp)
        #
        ####################################
        # get back the paralell runs and put in the main long chains_list
        for i in range(n_PT):
            chains_list[i].chain[i_tmp:i_nt,:] = chains_tmp[i].chain
            chains_list[i].lnprob[i_tmp:i_nt] = chains_tmp[i].lnprob
            # wrapper accpet_rate 
            nb_array = np.arange(start=i_tmp+1, stop=i_nt+1, step=1)
            if i_tmp > 0:
                ar_old = chains_list[i].accept_rate[i_tmp-1] 
                chains_list[i].accept_rate[i_tmp:i_nt] = (ar_old*i_tmp+chains_tmp[i].accept_rate)/nb_array
                ar_arr_tmp[i] = chains_tmp[i].accept_rate[-1]/n_iter_once
            else:
                ar_old = 0
                chains_list[i].accept_rate[i_tmp:i_nt] = (ar_old*i_tmp+chains_tmp[i].accept_rate)/nb_array
                ar_arr_tmp[i] = chains_tmp[i].accept_rate[-1]/n_iter_once
        #
        ####################################
        # set new x_init from the last iterm of previous chain_tmp 
        for i in range(n_PT):
            x_init_tmp[i] = chains_tmp[i].chain[-1,:] 
            #print("i xinit:", x_init_tmp[i])
        # 
        ###################################
        # add accept rate check and automatic proposal change
        if (i_tmp > (test_ar_start)) and ratio_true(tune_ar):
            #
            # find bad chains
            for i in range(n_PT):
                ar_arr[i] = chains_list[i].accept_rate[i_nt-1] 
            print("ar_array_cumul:", ar_arr)
            print("ar_array_tmp:", ar_arr_tmp)
            #
            # bad chains list
            i_PT_bad_list = find_bad_chain(ar_arr_tmp, ar_min_crit, ar_max_crit)
            #
            # if i_PT_bad_list not empty
            if i_PT_bad_list:
                #
                # create a scale array to tune para for all n_PT
                tune_arr = np.logspace(para_tune_min, para_tune_max, n_PT)
                #
                ################################################
                ################################################
                # tune for the sigma_prop_i_PT for each i_PT
                print("##################################################")
                print("ERR. bad chains:", i_PT_bad_list)
                f=open('bad_ar_chains.dat','a+')
                i_PT_bad_list_str = "%s%s" % (str(i_PT_bad_list),  "\n")
                f.write(i_PT_bad_list_str)
                f.close()
                # 
                # tuned sigma_prop_new
                sigma_prop_new = copy.deepcopy(sigma_prop)
                print("sigma_prop_before", sigma_prop)
                #
                for i in range(len(i_PT_bad_list)):
                    #
                    i_PT_totune = i_PT_bad_list[i]
                    #
                    i_PT_sigma_prop = sigma_prop[i,:]
                    print("i_PT_totune",i_PT_totune)
                    #
                    # create a scaled sigma matrix of n_PT*n_dim dimensition
                    sigma_prop_alltune = sigma_tune_2d_alltune(i_PT_sigma_prop, tune_arr, sigma_bound)
                    print("sigma_prop_alltune:",sigma_prop_alltune)
                    #
                    ################
                    # vary 1 para at a time
                    # create ar matrix  
                    ar_2d_alltune = np.zeros((n_PT,n_dim), dtype=float)
                    for j in range(n_dim):
                        #print(j, "_th para")
                        #
                        # n_PT simultaneous chains corespding to n_PT changes for each parameter
                        # create temp chain_seg 
                        chains_tmp = {x: class_mcmc.Chains(n_iter_tunesig,n_dim) for x in range(n_PT)}
                        #
                        # set sigma_prop  for each parameter
                        sigma_prop_1tune = sigma_tune_2d_1tune(sigma_prop_alltune, j, i_PT_sigma_prop)
                        #print("sigma_prop_1tune",sigma_prop_1tune)
                        #
                        # set init parameter set
                        for k in range(n_PT):
                            chains_tmp[k].chain[0,:]  = x_init_tmp[i_PT_totune] 
                        #print("bad i, bad xinit:", i_PT_totune, x_init_tmp[i_PT_totune])
                        # beta for param in tune
                        beta_PT_totune = beta_PT[i_PT_totune]
                        #
                        #############
                        # parallel chains
                        que.put(chains_tmp)
                        #
                        processes = [mp.Process(target=run_mcmc, args=(que, time_con, k, delta_dx_dy_1refstar, delta_dx_dy_sig, n_iter_tunesig, beta_PT_totune, sigma_prop_1tune[k,:])) for k in range(n_PT)]
                        # Run processes
                        for p in processes:
                            p.start()
                        # blocks the calling thread until the thread whose join() method is called is terminated.
                        for p in processes:
                            p.join()
                        chains_tmp = que.get(chains_tmp)
                        #
                        for k in range(n_PT):
                            ar_2d_alltune[k,j] = chains_tmp[k].accept_rate[-1]/n_iter_tunesig 
                        #
                    ##########
                    ## still i_PT loop
                    ## set new i_PT_sigma_prop
                    i_change, i_PT_use = judge_new_sigma(ar_2d_alltune)
                    #print("i_PT_use, i_change_para",i_PT_use, i_change)
                    #print("ar_2d_alltune",ar_2d_alltune)
                    print("i_PT, i_para, tuned to ar_tuned", i,i_change, ar_2d_alltune[i_PT_use, i_change])
                    # find good i_change and i_PT_use
                    if i_change >= 0:
                        sigma_prop_new[i_PT_totune,i_change] = sigma_prop_alltune[i_PT_use,i_change]
                    # bad i_change and i_PT_use, set to initial ratio
                    else: 
                        if i_PT_use >= 0:
                            sigma_prop_new[i_PT_totune,:] = gen_gaus_prop_sig()
                        else:
                            sigma_prop_new[i_PT_totune,:] = racist(sigma_prop_alltune, ar_2d_alltune)
                    #
                #############
                ## done tune all bad PT
                #print("loop_before_i", sigma_prop)
                sigma_prop = copy.deepcopy(sigma_prop_new)
                #print("last sigma_prop", sigma_prop)
            #
            # above crit ok
            else:
                print("OK",sigma_prop)
        #
        # first head iteration not check ok
        else:
            print("OK",sigma_prop)
        #
        ###############################
        #  swap_or_not for all 
        for i_swap in range(N_swap):
            # sequential swap_ test
            #swap_i = int(swap_i)
            #swap_j = int(swap_i + 1)
            # generante two random integer
            swap_i, swap_j = gen_two_rand(n_PT)
            # calc the swap probability
            do_swap = mh_crit_beta(x_init_tmp[swap_i], x_init_tmp[swap_j], beta_PT[swap_i], beta_PT[swap_j], time_con, log_likelihood, delta_dx_dy_sig, delta_dx_dy_1refstar) 
            if (do_swap):
                print("DO  swap,i,j: T,", swap_i, swap_j)
                # swap the initial parameter-set of two short chains
                x_init_tmp[swap_i], x_init_tmp[swap_j] = x_init_tmp[swap_j], x_init_tmp[swap_i] 
            else:
                print("NOT swap,i,j: F,", swap_i, swap_j)
        #
        # save sigma_prop
        for i in range(n_PT):
            fname = "%s%s%s" % ("sigma_prop_", i, ".dat")
            f=open(fname,'a+')
            sig_prop_i_str = "%s%s" % (str(sigma_prop[i,:]),  "\n")
            f.write(sig_prop_i_str) 
            f.close()
        #
        #
        ###############################
        i_tmp = i_nt
        #
    return chains_list

