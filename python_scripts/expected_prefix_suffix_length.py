
#!/usr/bin/env python

"""expected_prefix_suffix_length.py: See section 3.2 in Alternative Splicing document.
"""

__authors__      = "Yash Kumar Sonthalia, Israa Alqassem"
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
    

def set_global_vars_in_prefix_suffix(_strand_dir,_loc_index_dict,_start_sites_dict,_end_sites_dict,_subexon_ids_dict,_subexon_trans_ids_dict,_end_sites_list,_start_sites_list,_loc_list):
    
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


def get_length_list(read_length):
	'''
		This function computes the segment lengths according to the number of bases
			x_____y
			|_____|
			
		We have 4 cases here:
		
		1) if x and y are both start sites OR both are end sites, then:
			length = y - x
		2) if x is a start site y is an end site, then:
			length = y - x + 1
		3) if x is an end site and y is a start site, then:
			length = y - x - 1
			
	'''
	
	length_list = []
	for i in range(len(loc_list)-1):
		first = 's'
		second = 's'
		if loc_list[i] in end_sites_dict:
			first = 'e'
		if loc_list[i+1] in end_sites_dict:
			second = 'e'
		if first=='s' and second=='e':
			length_list.append(loc_list[i+1]-loc_list[i]+1)
		elif first=='e' and second=='s':
			length_list.append(loc_list[i+1]-loc_list[i]-1)
		else:
			length_list.append(loc_list[i+1]-loc_list[i]) ## s and s, e and e
	length_list = [abs(x) for x in length_list]
	
    # Subtract the read length from the last segment(s)
	read_length -= 1
	for j in range(len(length_list)-1, -1 ,-1):
		length_list[j] = length_list[j] - read_length 
		if length_list[j] >= 0: 
			break
		else:
			read_length = abs(length_list[j]) # avoid negative segment lengths
			length_list[j] = 0
    
    
	return length_list



#TO-DO : check division by error for both prefix and suffix length func's
def get_expected_prefix_length(length_list,PZ,f00,f01,f10,f11):
	'''
		This function computes the expected prefix lengths as described in 3.2.1
		Two types of prefix length for the i-th segment:
		lp_in[i]: the expected length of a subpath Z_1:i given that X_i is a part of an isoform; Z_i=1
		lp_out[i]: the expected length of a subpath Z_1:i given that X_i is not a part of an isoform Z_i=0
		
		This function returns two lists lp_in,lp_out for all possible i's
		
	'''
	lp_in = [0. for _ in range(len(length_list))]
	lp_out = [0. for _ in range(len(length_list))]
	lp_in[0] = length_list[0] #first index
	
	for i in range(1,len(length_list)):
		if PZ[i] == 0 or 1-PZ[i] == 0:
			lp_in[i] = 0
			lp_out[i] = 0
			print('get_expected_prefix_length: Division by zero for i =',i)
			continue
		if PZ[i] > 0 :
			lp_in[i] = length_list[i] + (lp_in[i-1]*PZ[i-1]*f11[i-1][i])/PZ[i] + (lp_out[i-1]*(1-PZ[i-1])*f01[i-1][i])/PZ[i]
		if PZ[i] != 1 : 
			lp_out[i] = (lp_in[i-1]*PZ[i-1]*f10[i-1][i])/(1-PZ[i]) + (lp_out[i-1]*(1-PZ[i-1])*f00[i-1][i])/(1-PZ[i])
	return lp_in,lp_out



def get_expected_suffix_length(length_list,PZ,f00,f01,f10,f11):
	'''
		This function computes the expected suffix lengths as described in 3.2.2
		Two types of prefix length for the i-th segment:
		ls_in[i]: the expected length of a subpath Z_i:M given that X_i is a part of an isoform; Z_i=1
		ls_out[i]: the expected length of a subpath Z_i:M given that X_i is not a part of an isoform Z_i=0
		
		This function returns two lists ls_in,ls_out for all possible i's
		
	'''	
	ls_in = [0. for _ in range(len(length_list))]
	ls_out = [0. for _ in range(len(length_list))]
	ls_in[-1] = length_list[-1] #last index
	for i in range(len(length_list)-2,-1,-1):
		ls_in[i] = length_list[i] + ls_in[i+1]*f11[i][i+1] + ls_out[i+1]*f10[i][i+1]
		ls_out[i] = ls_in[i+1]*f01[i][i+1] + ls_out[i+1]*f00[i][i+1]
	return ls_in,ls_out
