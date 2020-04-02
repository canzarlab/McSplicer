
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""populate_PZ.py: See section 3.1 in Alternative Splicing document."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


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
    

def set_global_vars_in_populate_PZ(_strand_dir,_loc_index_dict,_start_sites_dict,_end_sites_dict,_subexon_ids_dict,_subexon_trans_ids_dict,_end_sites_list,_start_sites_list,_loc_list):
    
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
    


def populate_PZ(q_arr,p_arr):
    
    """ See section 3.1 in Alternative Splicing document
    Compute P(Z_i) for all possible i
    Z = (Z_1,...,Z_M): a sequence of hidden R.Vs, where Z_i is an idicator of whether a 
    segment X_i is apart of an isoform, where M = M_s + M_e -1 (total num of end and start sites) 
    Parameters:
        loc_list:  all locations sorted based on strand direction
        q_arr: q1,q2,q3,..
        p_arr: p0,p1,p2,... Note that p[0] = pi
        start_sites_dict -> key: start location, index (helps to determine s1, s2, ...)
        end_sites_dict ->   key: end location, index (helps to determine e1, e2, ...)
     Returns:
     an array of P(Z_i = 1) where i in [0,M], P(Z_i=0) is the complement of this
    """
    M_s = len(start_sites_dict)
    M_e = len(end_sites_dict)
    M = M_s + M_e - 1     # The gene has M segments, see 1.1 notation in RNAsplicing.pdf    
    PZ = [0. for i in range(M)] 
    PZ[0] = p_arr[0]      # defined as pi in the doc
    
    
    
    i = 1
    for location in loc_list[1:-1]: # We skip the very first start site and the very last end site, see Figure 1 in RNASPlicing.pdf 
        if location in start_sites_dict.keys():
            PZ[i] = PZ[i-1]*1 + (1-PZ[i-1])*(p_arr[start_sites_dict[location]]) # See equation 29
        elif location in end_sites_dict.keys():
            PZ[i] = PZ[i-1]*(1-q_arr[end_sites_dict[location]]) + (1-PZ[i-1])*0 #TODO:- remove 0
        else:
            print "Error: populate_PZ"
        i += 1
        
    return PZ

