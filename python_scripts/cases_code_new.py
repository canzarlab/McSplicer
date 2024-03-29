
# coding: utf-8

# In[ ]:

#!/usr/bin/env python

"""cases_builders.py: See section 3.5 in Alternative Splicing document.
    make sure to call compute_f00_f01_f10_f11(p_arr, q_arr) before calling the functions here
"""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"



# In[ ]:

#from global_vars import *

err_case4 = 'Case 4 occurs'
err_case5 = 'Case 5 occurs'
err_unknown = 'Unknown case'


strand_dir = ''
loc_index_dict = []
start_sites_dict = {}
end_sites_dict = {}
subexon_ids_dict = {}
subexon_trans_ids_dict = {}
end_sites_list = []
start_sites_list = []
loc_list = []
    

def set_global_vars_in_cases_builders(_strand_dir,_loc_index_dict,_start_sites_dict,_end_sites_dict,_subexon_ids_dict,_subexon_trans_ids_dict,_end_sites_list,_start_sites_list,_loc_list):
    
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






# In[8]:

def get_index_in_loc_list(index_str):
    '''
        index_str has a format of s5 or e10
    '''
    if index_str[0:1] == 's':
        return loc_index_dict[start_sites_list[int(index_str[1:])]]
    else:
        return loc_index_dict[end_sites_list[int(index_str[1:])]]


# In[33]:

def builder_start_01(binary_list,start_end_list, F_s, sm, w_s,p_arr, q_arr, f00, f01, f10, f11):
    '''
        Section 3.5.1 - Alternative Splicing doc
        binary_list: the subpath binary list
        start_end_list: list of start and end indeces
        w_s: the conditional probability of subpath s on X_F(s) is a part of an isoform
        F_s:  first segment index of a subpath s
        sm: potential exon start site (location value)
        BIG NOTE: we assume sm belongs to location list, we don't consider the case 
        when sm is not an official start site, do we need to change that??
    '''
  
    #f00,f01,f10,f11 = compute_f00_f01_f10_f11(p_arr, q_arr)
    
    
    L_str = start_end_list[-1]    
    L_s = get_index_in_loc_list(L_str)-1             # last segment index of the subpath
    sm_index = start_sites_dict[sm]                  # index according to start site locations
    Ism = loc_index_dict[sm]-1                       # index according to all location
    ## |_____|_____|_____|
    ## Ism represents the index of the segment before sm
    
    
    Ism_str = 's'+str(sm_index) # the string index of sm according to start sites
        
    return_val = 0.
    
    #print start_end_list
    #print binary_list
    #print 'F_s =',F_s
    #print 'L_s =',L_s
    #print 'sm_index =',sm_index
    #print 'Ism =',Ism
    #print 'Ism_str =',Ism_str
    
    
    # Case1: exon start site sm appears left side of a subpath s
    if Ism < F_s:
        #print 'Case1: Ism < F_s'
        # equ 125
        try:
            return_val = ((1-p_arr[0])*f00[0][Ism] + p_arr[0]*f10[0][Ism])* p_arr[sm_index]*f11[Ism+1][F_s]*w_s
        except IndexError:
            print('Warning: Index out of range in case 125 (cases code new)')
            return_val = 0.            
            
    
    # Case2: exon start site sm appears right side of a subpath s
    elif L_s < Ism:
        #print 'Case2: L_s < Ism'
        # equ 131
        try:
            return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])* p_arr[sm_index]*f10[L_s][Ism]*w_s
        except IndexError:
            print('Warning: Index out of range in case 131 (cases code new)')
            return_val = 0.
        
    elif F_s <= Ism and Ism < L_s: 
        i = start_end_list.index(Ism_str)-1
        # check if (Z_Ism = 0, Z_Ism+1 = 1) is subset of s
        if binary_list[i]  == 0 and binary_list[i+1] == 1:
            #print  'Case3: F_s <= Ism and Ism < L_s'
            # equ 137
            try:
                return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*w_s
                
            except IndexError:                
                print('Warning: Index out of range in case 137 (cases code new)')
                return_val = 0.
        #else:
            #print "builder_start_01: "+err_case4
    #elif Ism == L_s:
        #print "builder_start_01: "+err_case5
    #else:
        #print "builder_start_01: "+err_unknown  

    return return_val


# In[35]:

def builder_start_00(binary_list,start_end_list, F_s, sm, w_s,p_arr, q_arr, f00, f01, f10, f11):
    '''
        Section 3.5.2 - Alternative Splicing doc
        binary_list: the subpath binary list
        start_end_list: list of start and end indeces
        w_s: the conditional probability of subpath s on X_F(s) is a part of an isoform
        F_s:  first segment index of a subpath s
        sm: potential exon start site (location value)
        BIG NOTE: we assume sm belongs to location list, we don't consider the case 
        when sm is not an official start site, do we need to change that??
    '''
  
    #f00,f01,f10,f11 = compute_f00_f01_f10_f11(p_arr, q_arr)
    
    L_str = start_end_list[-1]    
    L_s = get_index_in_loc_list(L_str)-1             # last segment index of the subpath
    sm_index = start_sites_dict[sm]                  # index according to start site locations
    Ism = loc_index_dict[sm]-1                       # index according to all location
    
    
    Ism_str = 's'+str(sm_index) # the string index of sm according to start sites
       
    return_val = 0.
    '''
    print start_end_list
    print binary_list
    print 'F_s =',F_s
    print 'L_s =',L_s
    print 'sm_index =',sm_index
    print 'Ism =',Ism
    print 'Ism_str =',Ism_str
    '''
    
    # Case1: exon start site sm appears left side of a subpath s
    if Ism < F_s-1:
        #print 'Case1: Ism < F_s'
        # equ 145
        try:
            return_val =  ((1-p_arr[0])*f00[0][Ism] + p_arr[0]*f10[0][Ism])*(1-p_arr[sm_index])*f01[Ism+1][F_s]*w_s
        except IndexError:                
            print('Warning: Index out of range in case 145 (cases code new)')
            return_val = 0.
            
        #print 'return_val=',return_val
    # Case2: exon start site sm appears right side of a subpath s
    elif L_s < Ism:
        #print 'Case2: L_s < Ism'
        # equ 151
        try:
            return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*(1-p_arr[sm_index])*f10[L_s][Ism]*w_s 
        except IndexError:                
            print('Warning: Index out of range in case 151 (cases code new)')
            return_val = 0.
            
        
    elif F_s <= Ism and Ism < L_s:
        i = start_end_list.index(Ism_str) - 1
        #print "i is",i
        # check if (Z_Ism = 0, Z_Ism+1 = 1) is subset of s
        if binary_list[i]  == 0 and binary_list[i+1] == 0:
            #print  'Case3: F_s <= Ism and Ism < L_s'
            # equ 157
            try:
                return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*w_s
            except IndexError:
                print('Warning: Index out of range in case 157 (cases code new)')
                return_val = 0.
            
            
        #else:
            #print "builder_start_00: "+err_case4
    #elif Ism == L_s:
        #print "builder_start_00: "+err_case5
    #else:
        #print "builder_start_00: "+err_unknown 

    return return_val


# In[36]:

def builder_end_11(binary_list,start_end_list, F_s, em, w_s,p_arr, q_arr, f00, f01, f10, f11):
    '''
        Section 3.5.3 - Alternative Splicing doc
        binary_list: the subpath binary list
        start_end_list: list of start and end indeces
        w_s: the conditional probability of subpath s on X_F(s) is a part of an isoform
        F_s:  first segment index of a subpath s
        em: potential exon end site (location value)
        BIG NOTE: we assume em belongs to location list, we don't consider the case 
        when em is not an official start site, do we need to change that??
    '''
    
    #f00,f01,f10,f11 = compute_f00_f01_f10_f11(p_arr, q_arr)
    
    L_str = start_end_list[-1]    
    L_s = get_index_in_loc_list(L_str)-1                # last segment index of the subpath
    em_index = end_sites_dict[em]                       # index according to end site locations
    Iem = loc_index_dict[em]-1                          # index according to all location
    
    
    Iem_str = 'e'+str(em_index) # the string index of em according to end sites
        
    return_val = 0.
    
    #print start_end_list
    #print binary_list
    #print 'F_s =',F_s
    #print 'L_s =',L_s
    #print 'em_index =',em_index
    #print 'Iem =',Iem
    #print 'Iem_str =',Iem_str
    
    
    # Case1: exon end site em appears left side of a subpath s
    if Iem < F_s:
        #print 'Case1: Iem < F_s'
        # equ 163
        try:
            return_val =  ((1-p_arr[0])*f01[0][Iem] + p_arr[0]*f11[0][Iem])*(1-q_arr[em_index])*f11[Iem+1][F_s]*w_s
            
        except IndexError:
            print('Warning: Index out of range in case 163 (cases code new)')
            return_val = 0.
        
        
    
    # Case2: exon end site em appears right side of a subpath s
    elif L_s <= Iem:
        #print 'Case2: L_s <= Iem'
        # equ 167
        try:
            return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*(1-q_arr[em_index])*f11[L_s][Iem]*w_s
            
        except IndexError:
            print('Warning: Index out of range in case 167 (cases code new)')
            return_val = 0.
            
        
    elif F_s <= Iem and Iem < L_s:
        i = start_end_list.index(Iem_str)-1
        # check if (Z_Iem = 1, Z_Iem+1 = 1) is subset of s
        if binary_list[i]  == 1 and binary_list[i+1] == 1:
            #print  'Case3: F_s <= Iem and Iem < L_s'
            # equ 172
            try:
                return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*w_s
                
            except IndexError:
                print('Warning: Index out of range in case 167 (cases code new)')
                return_val = 0.
            
            
        #else:
            #print "builder_end_11: "+err_case4
    #elif Iem == L_s:
        #print "builder_end_11: "+err_case5
    #else:
        #print "builder_end_11: "+err_unknown
    
    return return_val


# In[37]:

def builder_end_10(binary_list,start_end_list, F_s, em, w_s,p_arr, q_arr, f00, f01, f10, f11):
    '''
        Section 3.5.4 - Alternative Splicing doc
        binary_list: the subpath binary list
        start_end_list: list of start and end indeces
        w_s: the conditional probability of subpath s on X_F(s) is a part of an isoform
        F_s:  first segment index of a subpath s
        em: potential exon end site (location value)
        BIG NOTE: we assume em belongs to location list, we don't consider the case 
        when em is not an official start site, do we need to change that??
    '''
    
    #f00,f01,f10,f11 = compute_f00_f01_f10_f11(p_arr, q_arr)
    
    L_str = start_end_list[-1]    
    L_s = get_index_in_loc_list(L_str)-1                # last segment index of the subpath
    em_index = end_sites_dict[em]                       # index according to end site locations
    Iem = loc_index_dict[em]-1                          # index according to all location
    
    
    Iem_str = 'e'+str(em_index) # the string index of em according to end sites
        
    return_val = 0.
    
    #print start_end_list
    #print binary_list
    #print 'F_s =',F_s
    #print 'L_s =',L_s
    #print 'em_index =',em_index
    #print 'Iem =',Iem
    #print 'Iem_str =',Iem_str
    
    
    # Case1: exon end site em appears left side of a subpath s
    if Iem < F_s-1:
        #print 'Case1: Iem < F_s-1'
        # equ 177
        try:
            return_val =  ((1-p_arr[0])*f01[0][Iem] + p_arr[0]*f11[0][Iem])*q_arr[em_index]*f01[Iem+1][F_s]*w_s
        except IndexError:
            print('Warning: Index out of range in case 177 (cases code new)')
            return_val = 0.
    
    # Case2: exon end site em appears right side of a subpath s
    elif L_s <= Iem:
        #print 'Case2: L_s <= Iem'
        # equ 181
        try:
            return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*q_arr[em_index]*f11[L_s][Iem]*w_s
        except IndexError:
            print('Warning: Index out of range in case 181 (cases code new)')
            return_val = 0.   
           
        
    elif F_s <= Iem and Iem < L_s:
        i = start_end_list.index(Iem_str)-1
        # check if (Z_Iem = 1, Z_Iem+1 = 0) is subset of s
        if binary_list[i]  == 1 and binary_list[i+1] == 0:
            #print  'Case3: F_s <= Iem and Iem < L_s'
            # equ 185
            try:
                return_val =  ((1-p_arr[0])*f01[0][F_s] + p_arr[0]*f11[0][F_s])*w_s
                
            except IndexError:
                print('Warning: Index out of range in case 185 (cases code new)')
                return_val = 0.  

                
        #else:
            #print "builder_end_10: "+err_case4
    #elif Iem == F_s-1:
        #print "builder_end_10: "+err_case5
    #else:
        #print "builder_end_10: "+err_unknown

    return return_val


# In[ ]:



