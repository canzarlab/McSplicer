
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""populate_PZ.py: See section 3.3 in Alternative Splicing document."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"



import numpy as np
#from global_vars import *



  
strand_dir = ''
loc_index_dict = []
start_sites_dict = {}
end_sites_dict = {}
subexon_ids_dict = {}
subexon_trans_ids_dict = {}
end_sites_list = []
start_sites_list = []
loc_list = []
    

def set_global_vars_in_compute_fs(_strand_dir,_loc_index_dict,_start_sites_dict,_end_sites_dict,_subexon_ids_dict,_subexon_trans_ids_dict,_end_sites_list,_start_sites_list,_loc_list):
    
    global strand_dir
    global loc_index_dict
    global start_sites_dict
    global end_sites_dict
    global subexon_ids_dict
    global subexon_trans_ids_dict
    global end_sites_list
    global start_sites_list
    global loc_list
    
    strand_dir = _strand_dir
    loc_index_dict = _loc_index_dict
    start_sites_dict = _start_sites_dict 
    end_sites_dict = _end_sites_dict
    subexon_ids_dict = _subexon_ids_dict
    subexon_trans_ids_dict = _subexon_trans_ids_dict
    end_sites_list = _end_sites_list
    start_sites_list = _start_sites_list 
    loc_list = _loc_list 
    




def compute_f00_f01(p_arr, q_arr):
    """ Compute the following quantities for all valid i and j
        f00(i,j) := P(Z_j=0|Z_i=0) 
        f01(i,j) := P(Z_j=1|Z_i=0) 
        returns two matrices, f00 and f01
        ** We assume that loc list is sorted either ascending or descending based on strand dir
    """

    # Ignoring the first and last sites.
    M = len(loc_list) - 1   # The gene has M segments, see 1.1 notation in RNAsplicing.pdf
    
    f00 = [[0. for i in range(M)] for j in range(M)]
    f01 = [[0. for i in range(M)] for j in range(M)]
    
    for i in range(0,M):
        for j in range(i,M):
            location = loc_list[j]
            #print "i= "+str(i) + " j= "+str(j)
            #print "location = "+str(location)
            # Base case.
            if i == j-1:
                if location in start_sites_dict:
                    f00[i][j] = 1-p_arr[start_sites_dict[location]]
                    f01[i][j] = p_arr[start_sites_dict[location]]
                # segements j-1 and j seperated by end site
                elif location in end_sites_dict:
                    f00[i][j] = 1.
                    f01[i][j] = 0.
                else:
                    print "Error1: \nfile -> compute_f00_f01_f10_f11.py \nMethod -> compute_f00_f01"
            elif i < j-1:
                # segements j-1 and j seperated by start site
                if location in start_sites_dict:
                    f00[i][j] = f00[i][j-1]*(1-p_arr[start_sites_dict[location]])
                    f01[i][j] = f01[i][j-1] + f00[i][j-1]*p_arr[start_sites_dict[location]]
                # segements j-1 and j seperated by end site
                elif location in end_sites_dict:
                    f00[i][j] = f00[i][j-1] + f01[i][j-1]*q_arr[end_sites_dict[location]]
                    f01[i][j] = f01[i][j-1]*(1-q_arr[end_sites_dict[location]])
                else:
                    print "Error2: \nfile -> compute_f00_f01_f10_f11.py \nMethod -> compute_f00_f01"
            elif i == j:
                f00[i][j] = 1.
                f01[i][j] = 0. 
                    
     
    return f00,f01


# In[ ]:

def compute_f10_f11(p_arr, q_arr):
    """ Compute the following quantities for all valid i and j
        f10(i,j) := P(Z_j=0|Z_i=1) 
        f11(i,j) := P(Z_j=1|Z_i=1) 
        returns two matrices, f00 and f01
    """
    M = len(loc_list) - 1
    f10 = [[0. for i in range(M)] for j in range(M)]
    f11 = [[0. for i in range(M)] for j in range(M)]
    
    for i in range(0,M):
        for j in range(i,M):
            location = loc_list[j]
            #print "i= "+str(i) + " j= "+str(j)
            #print "location = "+str(location)
            # Base case.
            if i==j-1:
                if location in start_sites_dict:
                    f10[i][j] = 0.
                    f11[i][j] = 1.
                # segements j-1 and j seperated by end site
                elif location in end_sites_dict:
                    f10[i][j] = q_arr[end_sites_dict[location]]
                    f11[i][j] = 1-q_arr[end_sites_dict[location]]
                else:
                    print "Error3: \nfile -> compute_f00_f01_f10_f11.py \nMethod -> compute_f10_f11"
            elif i < j-1:
                # segements j-1 and j seperated by start site
                if location in start_sites_dict:
                    f10[i][j] = f10[i][j-1]*(1-p_arr[start_sites_dict[location]])
                    f11[i][j] = f11[i][j-1] + f10[i][j-1]*p_arr[start_sites_dict[location]]
                # segements j-1 and j seperated by end site
                elif location in end_sites_dict:
                    f10[i][j] = f10[i][j-1] + f11[i][j-1]*q_arr[end_sites_dict[location]]
                    f11[i][j] = f11[i][j-1]*(1-q_arr[end_sites_dict[location]])
                else:
                    print "Error4: \nfile -> compute_f00_f01_f10_f11.py \nMethod -> compute_f10_f11"
            elif i == j:
                f10[i][j] = 0.
                f11[i][j] = 1.
    return f10,f11


# In[ ]:

def compute_f00_f01_f10_f11(p_arr, q_arr):
    
    f00,f01 = compute_f00_f01(p_arr, q_arr)
    f10,f11 = compute_f10_f11(p_arr, q_arr)
    
    return f00,f01,f10,f11

