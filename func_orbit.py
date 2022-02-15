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

    


def func_as(time, ecc, osi, cos_inc, OmegaO, M0, omega, per):
    # alpha, e*sin(omega), e*cos(omega), i, Omega, M0, period
    omega=radians(omega)

    OmegaO=radians(OmegaO)
    M0=radians(M0)
    coso=cos(omega)
    sino=sin(omega)
    cosOg=cos(OmegaO)
    sinOg=sin(OmegaO)
    cosi=cos_inc
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
    

def cal_t_radec(time_con,mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega,acos_i_orbit, ascend_node_Omega, M0): 
    # calc osi_alpha, esino, ecoso, period 
    osi_alpha, period = calc_osi_etc(mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega)
    as_cal_con=func_as(time_con, e_orbit, osi_alpha, acos_i_orbit, ascend_node_Omega, M0, periapsis_omega, period)
    return as_cal_con
    
def cal_t_radec_fromAEI(M0, time_con,aei_file,mp_Mearth,ms_Msun,d_pc):
    #
    from inputsimobs import ms_Msun,d_pc
    from math import cos,radians
    import numpy as np
    #
    #NOTE
    M_j = M0 #t_seg_para[0,6]
    #
    # Nrows of aei file: 0day,1a,2e,3i,4periomega,5nodeOmega,6M,7mass
    Ndays_aei=aei_file.shape[0]
    # Nrows of time_con
    Ndays_timecon=time_con.shape[0]
    #
    # each time_con segment
    para_start = aei_file[0,0:7]
    # 
    # calc mp_Mearth
    # NOTE: this mass is generated from the aei file, may not = value in inputsimobs
    # divided by earth and moon(1/81earth)
    #print(aei_file[0,7])
    mp_Mearth = aei_file[0,7]/3.0034146856628466e-06 # divide by earth/solar
    #
    as_cal_con=np.zeros((len(time_con),2), dtype=np.float64)
    for i in range(len(time_con)-1):
        # interpol end time paras
        para_end = t_end_para(time_con[i+1], aei_file)
        #
        #print(para_start[0])
        #print(para_end[0])
        # derive all seg time para sequence
        between_min = int(time_con[i])+1
        between_max = int(time_con[i+1])
        t_seg_para = t_seg_create(para_start, para_end, between_min, between_max, aei_file)
        para_start = para_end
        #
        ###########################
        ##### per time step not work right now
        # should divide per time grid to sub_grid
        ## cal as_cal_con_varipara
        ###########################
        #
        # run t_seg_para as a whole list
        t_j = t_seg_para[:,0]
        a_j = t_seg_para[0,1]
        e_j = t_seg_para[0,2]
        i_j = t_seg_para[0,3]
        cos_i_j = cos(radians(i_j))
        perio_j = t_seg_para[0,4]
        nodeO_j = t_seg_para[0,5]
        #
        # tail, not work at large time gap
        #N_tail_arr = -3
        #t_seg_para_tail = t_seg_para[N_tail_arr:,:]
        #t_j = t_seg_para_tail[:,0]
        #print(t_seg_para_tail)
        #a_j = t_seg_para_tail[0,1]
        #e_j = t_seg_para_tail[0,2]
        #i_j = t_seg_para_tail[0,3]
        #cos_i_j = cos(radians(i_j))
        #perio_j = t_seg_para_tail[0,4]
        #nodeO_j = t_seg_para_tail[0,5]
        #
        osi_alpha, period = calc_osi_etc(mp_Mearth,ms_Msun,a_j,d_pc,e_j,perio_j)
        #print("osi,peri",osi_alpha,period)
        as_cal_con_seg=func_as(t_j, e_j, osi_alpha, cos_i_j, nodeO_j, M_j, perio_j, period)
        #print(as_cal_con_seg.shape)
        as_cal_con[i,:]=as_cal_con_seg[:,0]
        as_cal_con[i+1,:]=as_cal_con_seg[:,-1]
        #print(as_cal_con[i,:])
        #print(as_cal_con[i+1,:])
        #return as_cal_con
    as_mu = as_cal_con
    #np.savetxt("as_mu_2.dat", np.transpose([as_mu[0,:],as_mu[1,:]]))
    return as_mu


def t_seg_create(para_start, para_end, between_min, between_max, aei_file):
    import numpy as np
    n_t_seg = between_max-between_min+1+2
    #print(n_t_seg)
    n_para_add1 = aei_file.shape[1]-1
    t_seg_para = np.zeros((n_t_seg, n_para_add1), dtype=np.float64)
    t_seg_para[0,:] = para_start
    t_seg_para[-1,:] = para_end
    for i in range(n_t_seg-2):
        t_seg_para[i+1,:] = aei_file[between_min+i,0:n_para_add1]
    return t_seg_para


def t_end_para(t_end,aei_file):
    import numpy as np
    #print(t_end)
    #print(aei_file[:,0])
    end_err = aei_file[:,0] - t_end
    #print(end_err)
    #print(abs(end_err).argmin())
    j = abs(end_err).argmin()
    #print(aei_file[j-1:j+2,0])
    para_end = np.zeros(aei_file.shape[1]-1)
    para_end[0] = t_end
    for l in range(len(para_end)-1):
        para_end[l+1] = lagrange(t_end, aei_file[j-1:j+2,0], aei_file[j-1:j+2,l+1])
    return para_end



# Implementing Lagrange Interpolation
def lagrange(xp, x_arr, y_arr):
    #print(x_arr)
    #print(y_arr)
    if len(x_arr) != len(y_arr):
        quit()
    else:
        yp = 0
        for i in range(len(x_arr)):
            p = 1
            for j in range(len(y_arr)):
                if i != j:
                    p = p * (xp - x_arr[j])/(x_arr[i] - x_arr[j])
            yp = yp + p * y_arr[i] 
    return yp





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
    as_cal_con2 = cal_t_radec(time_con,x[13],ms_Msun,au2,d_pc,x[9],x[11],x[8],x[10],x[12])
    as_mu1 = np.transpose(as_cal_con)
    as_mu2 = np.transpose(as_cal_con2)
    as_mu = as_mu1 + as_mu2
    delta_dx_dy = create_delta_dxdy_pureAsMu(time_con, as_mu)
    return delta_dx_dy 

