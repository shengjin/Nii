import numpy as np
from inputsimobs import *

##########################
# observation time sequrence in days
# drop save random part
time_con0 = np.linspace(t0,t1,int(N_time_init))
#
seg_rand_min = 0.5
seg_rand_max = 1.5
fill_rand_min = 0.2
fill_rand_max = 1.8

def create_time_seg():
    N_miss_seg = np.random.randint(drop_seg_min, drop_seg_max+1) 
    N_miss_point = round(len(time_con0)*drop_ratio)
    N_miss_point_left = N_miss_point
    Len_N_miss_seg = np.zeros(N_miss_seg)
    for i in range(N_miss_seg-1):
        N_of_i_miss_seg = np.random.randint(round(N_miss_point_left/(N_miss_seg-i)*seg_rand_min), round(N_miss_point_left/(N_miss_seg-i)*seg_rand_max))
        Len_N_miss_seg[i] = N_of_i_miss_seg 
        N_miss_point_left = N_miss_point_left - N_of_i_miss_seg
    Len_N_miss_seg[-1] = N_miss_point_left
    #print(Len_N_miss_seg)
    #print(Len_N_miss_seg.sum())
    
    N_fill_point = len(time_con0)-int(Len_N_miss_seg.sum())
    N_fill_seg = N_miss_seg + 1
    N_fill_point_left = N_fill_point
    Len_N_fill_seg = np.zeros(N_fill_seg)
    for i in range(N_fill_seg-1):
        N_of_i_fill_seg = np.random.randint(round(N_fill_point_left/(N_miss_seg-i)*fill_rand_min), round(N_fill_point_left/(N_fill_seg-i)*fill_rand_max))
        Len_N_fill_seg[i] = N_of_i_fill_seg 
        N_fill_point_left = N_fill_point_left - N_of_i_fill_seg
    Len_N_fill_seg[-1] = N_fill_point_left
    #print(Len_N_fill_seg)
    #print(Len_N_fill_seg.sum())
    
    #print(len(time_con0))
    t_transit = np.zeros(len(time_con0))
    # fill N_miss to 1
    N_miss_start = 0
    N_miss_end = 0
    #print(N_miss_start,N_miss_end)
    for i in range(N_miss_seg):
        N_miss_start = N_miss_start + Len_N_fill_seg[i]
        N_miss_end = N_miss_start + Len_N_miss_seg[i]
        #print(N_miss_start,N_miss_end)
        t_transit[int(N_miss_start):int(N_miss_end)] = 1
        N_miss_start = N_miss_end
    
    #np.savetxt("t_transit.dat", np.transpose([t_transit]))
    retain_ind = np.where(t_transit < 1.0)
    return retain_ind
    
    
