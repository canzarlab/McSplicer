
# coding: utf-8

# In[1]:



#!/usr/bin/env python

"""populate_pe_ws.py: See section 3.5 
    w(sl,sr): The probability of a of read coming from a fragment starting with sl and ending with sr
"""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


# In[2]:

import argparse
import time
import os,sys,math,operator
import itertools
import scipy.stats
from pandas import *
import numpy as np
from subpath_info import *
from parse_cnt_file import parse_cnt_file
from parse_gtf_file import *
from new_gtf_genome_wide_parser import *


# In[3]:

strand_dir = ''
loc_index_dict = []
start_sites_dict = {}
end_sites_dict = {}
subexon_ids_dict = {}
subexon_trans_ids_dict = {}
end_sites_list = []
start_sites_list = []
loc_list = []
    

def set_global_vars_in_ws_pe(_strand_dir,_loc_index_dict,_start_sites_dict,_end_sites_dict,_subexon_ids_dict,_subexon_trans_ids_dict,_end_sites_list,_start_sites_list,_loc_list):
    
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


# In[4]:

def get_index_in_loc_list(index_str):
    '''
        index_str has a format of s5 or e10
    '''
    if index_str[0:1] == 's':
        return loc_index_dict[start_sites_list[int(index_str[1:])]]
    else:
        return loc_index_dict[end_sites_list[int(index_str[1:])]]


# In[5]:

def filter_subpath_list_by_subexon_ids(subpath_l_list,subpath_r_list, subpath_freq_list, subexon_ids):
    '''
        For paired-end read data
        get only the subpaths that belong to a sepcific gene (filtered by it's node ids)
    '''
    indices = []
    for i in range(len(subpath_l_list)):
        if set(subpath_l_list[i]) < set(subexon_ids) and set(subpath_r_list[i]) < set(subexon_ids):
            indices.append(i)
    subset_subpath_list_l = [subpath_l_list[idx] for idx in indices]
    subset_subpath_list_r = [subpath_r_list[idx] for idx in indices]
    subset_subpath_freq_list = [subpath_freq_list[idx] for idx in indices]
    return subset_subpath_list_l,subset_subpath_list_r,subset_subpath_freq_list


# In[6]:

def get_segment_length():
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
			length = y - x + 1
			
        This function will return a dictionary of segment id and it's length
        segment ids span from 0 to len(loc_list)-1
        
        |__0__|__1__|__2__| 
         
         The indices of segments start from 0
         
        The definition of a segment according to our model is all the base pairs between subexon boundries 
        s-s, s-e, e-s, e-e
	'''
	segment_length_dict = {}
	
	for i in range(len(loc_list)-1):
		first = 's'
		second = 's'
		if loc_list[i] in end_sites_dict:
			first = 'e'
		if loc_list[i+1] in end_sites_dict:
			second = 'e'
		if first=='s' and second=='e':
			segment_length_dict[i] = abs(loc_list[i+1]-loc_list[i]+1)
		elif first=='e' and second=='s':
			segment_length_dict[i] = abs(loc_list[i+1]-loc_list[i]-1)
		else:
			segment_length_dict[i] = abs(loc_list[i+1]-loc_list[i]) ## s and s, e and e
            
            
	return segment_length_dict


# In[59]:

def get_fragment_length(binary_read,s_e_read_seq):
    ''' For paired-end reads data
        if for example we have the following binary and start end lists:
            binary read = [1, 0, 0, 0, 1]
            s_e_read_seq = ['s8', 'e6', 's9', 'e7', 's10', 'e8']
        the length of the fragment would be the length of the first and last segments only (isoform is formed from these two segments)
    '''
    segment_length_dict = get_segment_length() # dict of segment id and its length
    frag_len = 0
    frag_len_list = []
    for i in range(len(binary_read)):
        #print binary_read[i]
        if binary_read[i] == 1:
            element = s_e_read_seq[i]
            if element[0:1] == 'e':
                
                genome_loc = end_sites_list[int(element[1:])]
                
                loc_idx = loc_index_dict[end_sites_list[int(element[1:])]]
                #print element,loc_idx,segment_length_dict[loc_idx]
                
                #print loc_idx,genome_loc,segment_length_dict[loc_idx]
                frag_len_list.append(segment_length_dict[loc_idx])
                frag_len += segment_length_dict[loc_idx]
            elif element[0:1] == 's':
                
                genome_loc = start_sites_list[int(element[1:])]
                
                loc_idx = loc_index_dict[start_sites_list[int(element[1:])]]
                
                #print loc_idx,genome_loc,segment_length_dict[loc_idx]
                
                frag_len_list.append(segment_length_dict[loc_idx])
                
                frag_len += segment_length_dict[loc_idx]
    return frag_len,frag_len_list


# In[36]:

def compute_wp_for_subpath(p_arr,q_arr,binary_node_list,start_end_node_list):
    '''
        This function is the same as populate ws function for single end reads
        w(s): denotes the conditional probability of S_n=s conditional on X_F(s) is a part of an isoform
    '''
    
    probability = 1.
    # If only one segment then probability is just 1.
    if len(binary_node_list) > 1:
        trimmed_start_end_node_list = start_end_node_list[1:-1]
        for i in range(len(binary_node_list)-1):
            segment_1 = binary_node_list[i]
            segment_2 = binary_node_list[i+1]
            type_site = trimmed_start_end_node_list[i][0]
            site_index = int(trimmed_start_end_node_list[i][1:])
            #print '>>>site_index',site_index
            if type_site == 's':
                if segment_1==0 and segment_2==0:
                    probability*=(1-p_arr[site_index])
                elif segment_1==0 and segment_2==1:
                    probability*=p_arr[site_index]
                elif segment_1==1 and segment_2==1:
                    probability*=1
                else:
                    probability=0.
                    #print "Error3: \nfile -> compute_ws.py \nPair of 1-0 segments is seperated by start site"    

            if type_site == 'e':
                if segment_1==1 and segment_2==1:
                    probability*=(1-q_arr[site_index])
                elif segment_1==1 and segment_2==0:
                    probability*=q_arr[site_index]
                elif segment_1==0 and segment_2==0:
                    probability*=1
                else:
                    probability=0.
                    #print "Error4: \nfile -> compute_ws.py \nPair of 0-1 segments seperated by end site"
    return probability


# In[9]:

def w_prime(p_arr,q_arr,subpath_l,subpath_r,idx1,val1,val2):
    '''
        Compute w_prime (wp') as in equation 123
    '''
    frag_len_dist = scipy.stats.norm(250, 50) # fragment length distribution

    if (subpath_l == subpath_r) or (subpath_l[-1] >= subpath_r[0]-1):
        print 'Error:\nNo gap between the left and right reads'
        return -1
    else:

        _, start_end_s_tmp = get_subpath_info(get_subexon_inbetween(subpath_l,subpath_r))
        binary_s_l, start_end_s_l =  get_subpath_info(subpath_l)
        binary_s_r, start_end_s_r = get_subpath_info(subpath_r)
        start_end_list = start_end_s_l + start_end_s_tmp + start_end_s_r
        start_end_list = list(unique(start_end_list))

        #print '>>start_end_s_tmp',start_end_list

        ## get the index of idx1 element in start_end_location
        for k in range(len(start_end_list)):
            print start_end_list[k],get_index_in_loc_list(start_end_list[k])
            if idx1 == get_index_in_loc_list(start_end_list[k]):
                break

        #print start_end_list[k],k

        # 2^N tuples :- number of segments between the left and right reads
        n = len(start_end_list)-len(start_end_s_l)-len(start_end_s_r)+1
        tuple_list = map(list, itertools.product([0, 1], repeat=n))

        # Get all the unique binary lists that meet the criteria (idx1=val,idx2=val2)
        binary_subpath_list = []
        for _tuple in tuple_list:
            binary_list = binary_s_l + _tuple + binary_s_r
            binary_list[k-1] = val1
            binary_list[k] = val2
            binary_subpath_list.append(binary_list)

        binary_subpath_list = [list(x) for x in set(tuple(x) for x in binary_subpath_list)]
        #print binary_subpath_list

        sum_p = 0.
        for binary_list in binary_subpath_list:
            subpath_length,_ = get_fragment_length(binary_list,start_end_list)
            frag_length_pdf = frag_len_dist.pdf(subpath_length)
            prob_v = compute_wp_for_subpath(p_arr,q_arr,binary_list,start_end_list)
            sum_p += (prob_v*frag_length_pdf)

    return sum_p


# In[10]:

def get_subexon_inbetween(s_l,s_r):
    '''
        s_l: the left read as a list of subexon ids
        s_r: the right read as a list of subexon ids
        
        Returns a list of subexon ids (node ids) between the left and right reads
    
    '''
 
    if strand_dir == '+':
        #subexon_ids = sorted(subexon_ids_dict.keys())
        end_l = max(s_l)
        start_r = min(s_r)
        
    elif strand_dir == '-':
        #subexon_ids = sorted(subexon_ids_dict.keys(), reverse=True)
        end_l = min(s_l)
        start_r = max(s_r)
        

    subexons_in_between = []
    for subexon_id in subexon_ids_dict.keys():
        
        if strand_dir == '+':
            if subexon_id > end_l and subexon_id < start_r:
                subexons_in_between.append(subexon_id)
                
        elif strand_dir == '-':
            if subexon_id < end_l and subexon_id > start_r:
                subexons_in_between.append(subexon_id)
            
    return subexons_in_between


# In[12]:

def populate_pe_ws(q_arr,p_arr,subpath_l_list,subpath_r_list,fragment_mean,fragment_sd):


    if len(subpath_l_list) != len(subpath_r_list):
        print 'Error: the number of left paired end reads doesn\'t match the number of right reads! '
        return []
    

    frag_len_dist = scipy.stats.norm(fragment_mean, fragment_sd) # fragment length distribution

    w = []
    
    frag_length_list = []

    for i in range(len(subpath_r_list)):
        
        if (strand_dir == '+'):
            subpath_l = sorted(subpath_l_list[i])
            subpath_r = sorted(subpath_r_list[i])
        else:
            subpath_l = sorted(subpath_l_list[i],reverse=True)
            subpath_r = sorted(subpath_r_list[i],reverse=True)
            
            
        # Get the min and max possible fragment length for a specific subpath  
        binary_list_l,start_end_list_l = get_subpath_info(subpath_l)
        binary_list_r,start_end_list_r = get_subpath_info(subpath_r)
        _,left_frag_len_list = get_fragment_length(binary_list_l,start_end_list_l)
        _,right_frag_len_list = get_fragment_length(binary_list_r,start_end_list_r)
        min_frag_length, max_frag_length = get_min_max_frag_length(subpath_l,subpath_r,left_frag_len_list,right_frag_len_list,read_length=100)

    
        # Case1: Left and right reads are obtained from the same exon
        if (subpath_l == subpath_r):

            subpath = subpath_l # either left or right (The same)
            binary_list,start_end_list = get_subpath_info(subpath)
            prob_v = compute_wp_for_subpath(p_arr,p_arr,binary_list,start_end_list)
            
            #frag_length_pdf = frag_len_dist.pdf(subpath_length)
            
            #w.append(prob_v*frag_length_pdf)
            
            #print 'subpath=',subpath
            #print 'subpath_length=',subpath_length,'prob_v=',prob_v,'frag_length_pdf=',frag_length_pdf

        # Case2: Left and right reads intersects or right read is a continuation of the left read
        elif (strand_dir == '+' and subpath_l[-1] >= subpath_r[0]-1) or            (strand_dir == '-' and subpath_l[-1] <= subpath_r[0]-1):

            subpath = list(unique(subpath_l+subpath_r))
            binary_list,start_end_list = get_subpath_info(subpath)
            
            prob_v = compute_wp_for_subpath(p_arr,p_arr,binary_list,start_end_list)
            
            #frag_length_pdf = frag_len_dist.pdf(subpath_length)
            
            #w.append(prob_v*frag_length_pdf)
            
            
            #print 'subpath=',subpath
            #print 'subpath_length=',subpath_length,'prob_v=',prob_v,'frag_length_pdf=',frag_length_pdf

        else:
            # Case3: There's a gab between left read and right read

            # Get the left and right reads and the path between them
            _, start_end_s_tmp = get_subpath_info(get_subexon_inbetween(subpath_l,subpath_r))
            
            start_end_list = start_end_list_l + start_end_s_tmp + start_end_list_r
            start_end_list = list(unique(start_end_list))

            sum_p = 0. # sum of probabilities of each possible subpath between left and right reads

            # 2^N tuples :- number of segments between the left and right reads
            n = len(start_end_list)-len(start_end_s_l)-len(start_end_s_r)+1
            tuple_list = map(list, itertools.product([0, 1], repeat=n))

            ## TODO-; check list length
            for _tuple in tuple_list:

                binary_list = binary_list_l + _tuple + binary_list_r
                subpath_length,frag_len_list = get_fragment_length(_tuple,start_end_s_tmp)
                
                #frag_length_pdf = frag_len_dist.pdf(subpath_length)
                
                prob_v = compute_wp_for_subpath(p_arr,q_arr,binary_list,start_end_list)
                
                #sum_p += (prob_v*frag_length_pdf)

                #print 'subpath=',subpath_l,subpath_r
                #print 'subpath_length=',subpath_length,'prob_v=',prob_v,'frag_length_pdf=',frag_length_pdf

            
            #print 'sum_p=',sum_p
            w.append(sum_p)
    return w,frag_length_list


# In[51]:

#if __name__ == "__main__":


DIR = '../data/'
gtf_file = 'Drosophila_melanogaster.BDGP6.89.chr.transid.refined.gtf'
gene_id = 'FBgn0264695' # gene_name = Mhc
paired_end = True
cnt_file_pe = 'out_paired_end.cnt'

all_gene_dict = get_all_genes_dict(DIR+gtf_file)
subpath_l_list,subpath_r_list, subpath_freq_list = parse_cnt_file(DIR+cnt_file_pe,paired_end)

strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict = get_gene_data(gene_id,all_gene_dict[gene_id])

if strand_dir == '+' :
    end_sites_list = sorted(end_sites_dict.keys())
    start_sites_list = sorted(start_sites_dict.keys())
    loc_list = sorted(loc_index_dict.keys())
elif strand_dir == '-' :
    end_sites_list = sorted(end_sites_dict.keys(), reverse=True)
    start_sites_list = sorted(start_sites_dict.keys(), reverse=True)
    loc_list = sorted(loc_index_dict.keys(), reverse=True)
else:
    print 'Error: unknown strand direction',strand_dir



######################  Set global variables in this file and in subpath info 
set_global_vars_in_ws_pe(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)
set_global_vars_in_subpath(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)



######################  Extract only the reads that belong to this gene
subpath_l_list,subpath_r_list, subpath_freq_list = filter_subpath_list_by_subexon_ids(subpath_l_list,subpath_r_list, subpath_freq_list, subexon_ids_dict.keys())

