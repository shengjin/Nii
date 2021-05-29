from constant import *
from numpy import sin,cos,pi,sqrt,arctan2,radians,array #,tan,arcsin,arccos,arctan,degrees,median,mean,array,log
#from numpy import sin,cos,tan,arcsin,arccos,arctan,sqrt,pi,radians,degrees,arctan2,median,mean,array,log



def calc_osi_etc(mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega):
    ##############
    # alpha, e*sin(omega), e*cos(omega), i, Omega, M0, period
    # calc osi_alpha
    osi_alpha = 3.0*(mp_Mearth)*(ms_Msun**-1.0)*(a_AU)*(d_pc**-1) # mu as
    #osi_alpha = 5.0*(mplanet_JupiterMass)*(msun_SolarMass**-1.0)*(a_5AU)*(d_pc**-1) # mas
    # calc orbital period
    period_second = 2*PI*((a_AU*AU2cm)**3.0/(GG*Msun*ms_Msun))**0.5 #second
    period = period_second/3600.0/24.0 #/year
    return osi_alpha, period

def newton(EE,e,M):
    # overflow encountered in double_scalars
    import numpy as np
    eps=1
    eps = np.float64(eps)
    while(abs(eps)>1e-8):
        E1=EE-(EE-e*sin(EE)-M)/(1-e*cos(EE))
        eps=E1-EE
        EE=E1
    return EE

def newton_solver(M, e, EE0=None):
    import numpy as np
    tolerance=1e-9
    max_iter=100
    M = np.array(M)
    if EE0 is None:
        EE = np.copy(M)
    else:
        EE = np.copy(EE0)
    EE -= (EE - (e * sin(EE)) - M) / (1.0 - (e * cos(EE)))
    diff = (EE - (e * sin(EE)) - M) / (1.0 - (e * cos(EE)))
    abs_diff = abs(diff)
    ind = np.where(abs_diff > tolerance)
    niter = 0
    while ((ind[0].size > 0) and (niter <= max_iter)):
        EE[ind] -= diff[ind]
        if niter == (max_iter//2):
            EE[ind] = pi
        diff[ind] = (EE[ind] - (e * sin(EE[ind])) - M[ind]) / (1.0 - (e * cos(EE[ind])))
        abs_diff[ind] = abs(diff[ind])
        ind = np.where(abs_diff > tolerance)
        niter += 1
    return EE


def func_as(time, ecc, osi, inc, OmegaO, M0, omega, per):
    # alpha, e*sin(omega), e*cos(omega), i, Omega, M0, period
    omega=radians(omega)
    inc=radians(inc)
    OmegaO=radians(OmegaO)
    M0=radians(M0)
    coso=cos(omega)
    sino=sin(omega)
    cosOg=cos(OmegaO)
    sinOg=sin(OmegaO)
    cosi=cos(inc)
    A=osi*(coso*cosOg-sino*sinOg*cosi)
    B=osi*(coso*sinOg+sino*cosOg*cosi)
    F=osi*(-sino*cosOg-coso*sinOg*cosi)
    G=osi*(-sino*sinOg+coso*cosOg*cosi)
    E_as= newton_solver(((2*pi)/per*time-M0), ecc) 
    #E_as=array(list(map(lambda x:newton(x,ecc,x),(2*pi)/per*time-M0)))
    #print("A,B,F,G",A,B,F,G)
    X=cos(E_as)-ecc
    Y=sqrt(1-ecc**2)*sin(E_as)
    ra=(B*X+G*Y)
    dec=(A*X+F*Y)
    return array([ra,dec])
    

def cal_t_radec(time_con,mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega,i_orbit, ascend_node_Omega, M0): 
    # calc osi_alpha, esino, ecoso, period 
    osi_alpha, period = calc_osi_etc(mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega)
    as_cal_con=func_as(time_con, e_orbit, osi_alpha, i_orbit, ascend_node_Omega, M0, periapsis_omega, period)
    return as_cal_con
    





##############################


def period_2_au(period_days, ms_Msun):
    AU2cm = 1.4959787e13   # Astronomical unit [cm]
    GG    = 6.67259e-8     # Gravitational constant
    Msun  = 1.9891e33      # Solar mass [g]
    PI    = 3.14159265358
    period_second = period_days*24.0*3600.0
    #period_second = 2*PI*((a_AU*AU2cm)**3.0/(GG*Msun*ms_Msun))**0.5 #second
    # todo M/m
    a_AU = ( (period_second/2.0/PI)**2.0*(GG*Msun*ms_Msun) )**(1.0/3.0) / AU2cm 
    #a_AU = ( (period_second/2.0/PI)**2.0*(GG*Msun*(ms_Msun+mp_Msun)) )**(1.0/3.0) / AU2cm 
    return a_AU

def au_2_period(a_AU, ms_Msun):
    AU2cm = 1.4959787e13   # Astronomical unit [cm]
    GG    = 6.67259e-8     # Gravitational constant
    Msun  = 1.9891e33      # Solar mass [g]
    PI    = 3.14159265358
    period_second = 2*PI*((a_AU*AU2cm)**3.0/(GG*Msun*ms_Msun))**0.5 #second
    period = period_second/3600.0/24.0 #/year
    return period

def mp_osi_au(a_AU, osi, ms_Msun, d_pc):
    #osi = 3.0*(mp_Mearth)*(ms_Msun**-1.0)*(a_AU)*(d_pc**-1) # mu as
    mp_Mearth = osi / 3.0 / (ms_Msun**-1.0) / (a_AU) / (d_pc**-1) 
    return mp_Mearth

#mp_Mearth_got =  mp_osi_au(a_AU_got, osi_max, ms_Msun, d_pc)





def gen_dx_dy_mc(x, time_con):
    from inputsimobs import ms_Msun,d_pc
    from func_orbit import period_2_au, cal_t_radec
    from func_ref import create_delta_dxdy_pureAsMu
    import numpy as np
    au = period_2_au(x[7], ms_Msun)
    as_cal_con = cal_t_radec(time_con,x[5],ms_Msun,au,d_pc,x[1],x[3],x[0],x[2],x[4])
    as_mu = np.transpose(as_cal_con)
    delta_dx_dy = create_delta_dxdy_pureAsMu(time_con, as_mu)
    return delta_dx_dy 

def gen_dx_dy_mc_2p(x, time_con):
    from inputsimobs import ms_Msun,d_pc
    from func_orbit import period_2_au, cal_t_radec
    from func_ref import create_delta_dxdy_pureAsMu
    import numpy as np
    au = period_2_au(x[7], ms_Msun)
    as_cal_con = cal_t_radec(time_con,x[5],ms_Msun,au,d_pc,x[1],x[3],x[0],x[2],x[4])
    au2 = period_2_au(x[14], ms_Msun)
    as_cal_con2 = cal_t_radec(time_con,x[13],ms_Msun,au,d_pc,x[9],x[11],x[8],x[10],x[12])
    as_mu1 = np.transpose(as_cal_con)
    as_mu2 = np.transpose(as_cal_con2)
    as_mu = as_mu1 + as_mu2
    delta_dx_dy = create_delta_dxdy_pureAsMu(time_con, as_mu)
    return delta_dx_dy 

