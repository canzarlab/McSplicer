
# coding: utf-8

# In[1]:

"""
test_EM.py: We used this file to test and debug EM algorithm 
using a single gtf and a single cnt files
"""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


# In[2]:



import os
import math
from global_vars import *
from subpath_info import get_binary_and_startend_list_for_all_subpaths
from compute_f00_f01_f10_f11 import compute_f00_f01_f10_f11
from populate_PZ import populate_PZ
from populate_ws import populate_ws
from expected_prefix_suffix_length import *
from cases_code_new import *


# In[3]:

cnt_file = DIR + gtf_filename.split('.')[0].upper()+'/' + gtf_filename.split('.')[0]+'_0.cnt'
subpath_list, subpath_freq_list = parse_cnt_file(DIR, cnt_file)
binary_subpath_list, start_end_nodes_subpath_list = get_binary_and_startend_list_for_all_subpaths(subpath_list)
p_arr = [0.2] * len(start_sites_dict)
q_arr = [0.2] * (len(end_sites_dict)-1)


# In[4]:

all_start_end_list = []
for i in range(len(start_sites_dict)):
    all_start_end_list.append('s'+str(i))
for i in range(len(end_sites_dict)):
    all_start_end_list.append('e'+str(i))
    
length_list = get_length_list()


no_subpaths = len(binary_subpath_list)
Fs_list = [get_index_in_loc_list(x[0]) for x in start_end_nodes_subpath_list]


# In[5]:

prev_likelihood = 0.
curr_likelihood = 0.
num_iter = 1000
step = 0
count = 0
max_count = 3
epsilon = 1e-5
offset_val = 0.0001 # to avoid division by zero error


# In[6]:

while step < num_iter:
    prev_likelihood = curr_likelihood
    curr_likelihood = 0.
    
    PS = [0. for _ in range(no_subpaths)]
    A = [0. for _ in range(len(start_sites_dict))]
    B = [0. for _ in range(len(start_sites_dict))]
    C = [0. for _ in range(len(end_sites_dict)-1)]
    D = [0. for _ in range(len(end_sites_dict)-1)]
    E = 0.
    F = 0.
    numerator_e = 0.
    numerator_f = 0.
    
    
    f00,f01,f10,f11 = compute_f00_f01_f10_f11(p_arr, q_arr)
    PZ = populate_PZ(q_arr,p_arr)
    w = populate_ws(q_arr,p_arr,binary_subpath_list, start_end_nodes_subpath_list,subpath_list)
    ls_in,ls_out = get_expected_suffix_length(length_list,PZ,f00,f01,f10,f11)
    lp_in,lp_out = get_expected_prefix_length(length_list,PZ,f00,f01,f10,f11)
    
    for i in range(no_subpaths):
        start_end_list = start_end_nodes_subpath_list[i]
        freq = subpath_freq_list[i]
        PS[i] = (1-p_arr[0])*f01[0][Fs_list[i]]*w[i] + (p_arr[0])*f11[0][Fs_list[i]]*w[i]
        
        #print "start_end_list",start_end_list
        #print "PS[i]",PS[i]
        #print "Fs_list[i]",Fs_list[i]
        #print "i = ",i
        #print "PS[i] = ", PS[i]
        #print "f01[1][Fs_list[i]] =",f01[1][Fs_list[i]]
        #print "f11[1][Fs_list[i]] = ",f11[1][Fs_list[i]]
        #print "w[i] = ",w[i]
        
        if Fs_list[i]==0:
            numerator_e+=freq*1.
        else:
            numerator_e+=freq*(p_arr[0]*f11[0][Fs_list[i]]*w[i])/PS[i]
            numerator_f+=freq*((1-p_arr[0])*f01[0][Fs_list[i]]*w[i])/PS[i]
        for start_end_value in all_start_end_list:#all_start_end_list:
            if start_end_value == 's0' or start_end_value == ('e'+str(len(end_sites_dict)-1)):
                continue
            elif start_end_value[0] == 's':
                start_value = start_end_value
                m = int(start_value[1:])
                sm = start_sites_list[m] # location val
                Ism = get_index_in_loc_list(start_value)-1
                A[m] += freq * builder_start_01(binary_subpath_list[i], start_end_list, Fs_list[i], sm, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]                
                B[m] += freq * builder_start_00(binary_subpath_list[i], start_end_list, Fs_list[i], sm, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]
            elif start_end_value[0] == 'e':
                end_value = start_end_value
                m = int(end_value[1:])
                em = end_sites_list[m] # location val
                Iem = get_index_in_loc_list(end_value)-1
                
                C[m] += freq * builder_end_10(binary_subpath_list[i], start_end_list, Fs_list[i], em, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]
                D[m] += freq * builder_end_11(binary_subpath_list[i], start_end_list, Fs_list[i], em, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]
                                
            else:
                print "Unknown error"
        



    E = numerator_e/ls_in[0] 
    F = numerator_f/ls_out[0]
    p_arr[0] = (E+offset_val)/(E+F+2*offset_val)

    for m in range(1,len(start_sites_dict)):
        Ism = get_index_in_loc_list('s'+str(m))-1
        A[m] /= (lp_out[Ism] + ls_in[Ism+1])
        B[m] /= (lp_out[Ism] + ls_out[Ism+1])
        p_arr[m] = (A[m]+offset_val)/(A[m]+B[m]+2*offset_val)


    for m in range(len(end_sites_dict)-1):
        Iem = get_index_in_loc_list('e'+str(m))-1
        C[m] /= (lp_in[Iem] + ls_out[Iem+1])
        D[m] /= (lp_in[Iem] + ls_in[Iem+1])
        q_arr[m] = (C[m]+offset_val)/(C[m]+D[m]+2*offset_val)
    
    for i in range(no_subpaths):
        #if freq * PS[i] > 0 :
        curr_likelihood += math.log(freq * PS[i]) 
    
    print "step", step
    print "p_arr", p_arr
    print "q_arr", q_arr
    print "Prev likelihood", prev_likelihood
    print "Curr likelihood", curr_likelihood
    
    
    
    if curr_likelihood-prev_likelihood <= epsilon:
        count += 1
    else:
        count = 0
        
    step += 1
    
    if count == max_count:
    	print "Breaking after %d iterations"%step
        break
    
    #if curr_likelihood-prev_likelihood <= epsilon:
    #    break
    


# In[7]:

print start_sites_dict
print p_arr
print end_sites_dict
print q_arr


# In[8]:

#print start_end_nodes_subpath_list
#print binary_subpath_list
#print Fs_list

