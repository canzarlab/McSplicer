
# coding: utf-8

# In[14]:



#!/usr/bin/env python

"""subpath_info.py: This file contains functions 
used to compute subpath binary and start_end list."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


# In[15]:

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
    

def set_global_vars_in_subpath(_strand_dir,_loc_index_dict,_start_sites_dict,_end_sites_dict,_subexon_ids_dict,_subexon_trans_ids_dict,_end_sites_list,_start_sites_list,_loc_list):
    
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



# In[16]:

def find_nearest_start_loc(element):
    """ Find the start location that's located exactly before the passed element (location)
    """ 
    if strand_dir == '+' :
        for i in range(1,len(loc_list)):
            if loc_list[i] >= element:
                return loc_list[i-1]
    else:
        for i in range(1,len(loc_list)):
            if loc_list[i] <= element:
                return loc_list[i-1]
    
    return loc_list[-1]


# In[17]:

def find_nearest_end_loc(element):
    """ Find the end location that's located exactly after the passed element (location)
    """ 
    if strand_dir == '+' :
        for i in range(1,len(loc_list)):
            if loc_list[i] >= element:
                return loc_list[i]
    else:
        for i in range(1,len(loc_list)):
            if loc_list[i] <= element:
                return loc_list[i]
            
    return loc_list[0]


# In[18]:

def compute_binary_and_startend_list(first_most_start,last_most_end,subpath_loc_list):
    """ 
    We assume that subpath_loc_list is sorted
    
    For a subpath this function returns two lists, the binary list and start and end sites list 
    e.g., [1, 0, 0, 0, 1, 1] 
          ['s1', 'e1', 's2', 'e2', 's3', 's4', 'e3']
    """ 
    
    num_of_segments = loc_index_dict[last_most_end] - loc_index_dict[first_most_start]
    binary_node_list = [0]*(num_of_segments)


    # loop through either start or end sites
    # i -> iterator over all locations starting from subpath_start to subpath_end
    # j -> iteratot over all locations within subpath (these locations aren't necessary 'real' ones)
    i = loc_index_dict[first_most_start] 
    j = 0
    count = 0
    while i < loc_index_dict[last_most_end] and j < len(subpath_loc_list):
        
        # cond: if this subexon coordinates lies within this segment
        # (+):  i___x___i+1 , x is the subexon lies within a segment
        # (-):  i+1___x___i
        if strand_dir=='+':
            cond1 = subpath_loc_list[j] >= loc_list[i]  and subpath_loc_list[j] < loc_list[i+1] 
            ## Fix for subexon that starts and ends at the same position
            cond2 =  subpath_loc_list[j] >= loc_list[i] and subpath_loc_list[j] == subpath_loc_list[j+1] and subpath_loc_list[j] <= loc_list[i+1]
            
            cond = cond1 or cond2
            
        else:
            cond1 = subpath_loc_list[j] <= loc_list[i]  and  subpath_loc_list[j] > loc_list[i+1]
            cond2 = subpath_loc_list[j] <= loc_list[i] and subpath_loc_list[j] == subpath_loc_list[j+1] and subpath_loc_list[j] >= loc_list[i+1]
            
            cond = cond1 #or cond2 
            
        if cond: 
            # set that segment to 1 (if it's not already 1)
            if binary_node_list[count]!=1:
                binary_node_list[count] = 1
                count +=1
                i+=1
            j+=2 # move to the next subexon coordinates (start and end)
        else:
             # If no other subexon lies in here then set this segment to 0
            if binary_node_list[count]!=1:
                binary_node_list[count] = 0
                count +=1
            i+=1
            
    
    start_end_node_list = []
    try:
        for i in range(loc_index_dict[first_most_start],loc_index_dict[last_most_end]+1):
            if i < len(loc_list):
                if loc_list[i] in start_sites_dict:
                    start_end_node_list.append('s'+str(start_sites_dict[loc_list[i]]))
                elif loc_list[i] in end_sites_dict:
                    start_end_node_list.append('e'+str(end_sites_dict[loc_list[i]]))
                else:
                    print("Error2: \nfile -> subpath_info.py \nUndefined Location-> ",loc_list[i])

    except IndexError as error:
        print('subpath_info.py:')
        print(error)
        #print 'i=%d, j=%d'%(i,j)
            


    return binary_node_list,start_end_node_list


# In[19]:

def get_subpath_info(subpath):
    '''
        subpath: contains a list of node ids in that subpath (an entry in cnt file)
        Similar to compute_binary_and_startend_list this function returns two lists, 
        the binary list and start and end sites list, but we only need to pass the subpath as a 
        parameter
            e.g., [1, 0, 0, 0, 1, 1] 
                  ['s1', 'e1', 's2', 'e2', 's3', 's4', 'e3']
    '''
    
    subpath_loc_list = []          # A list of start and end sites of all subexons within this read span
    subpath_nodes_start_sites = [] # A list of start locations of nodes within this read span
    subpath_nodes_end_sites = []   # A list of end locations of nodes within this read span


    for node in subpath:
        subpath_nodes_start_sites.append(subexon_ids_dict[node][0])
        subpath_nodes_end_sites.append(subexon_ids_dict[node][1])
        subpath_loc_list.extend(subexon_ids_dict[node])


    # Get the boundaries of this read span, i.e.,
    # the shortest subpath from which this read span (s) is derived
    
    if strand_dir == '+':
        first_most_start = min(subpath_nodes_start_sites)
        last_most_end = max(subpath_nodes_end_sites)
        subpath_loc_list.sort()
    else:
        first_most_start = max(subpath_nodes_start_sites)
        last_most_end = min(subpath_nodes_end_sites)
        subpath_loc_list.sort(reverse=True)

    # In caes that the subexon's start/end coordinate is not a "real" start/end site,
    # find the nearest start and end based on strand direction
    # last most end can be a start site (Stefan confirmed), that'w why we pass loc_list
    if first_most_start not in start_sites_dict:
        first_most_start = find_nearest_start_loc(first_most_start)

    if last_most_end not in end_sites_dict:
        last_most_end = find_nearest_end_loc(last_most_end)


    if last_most_end not in loc_index_dict:
        print("Error1: \nfile -> subpath_info.py \nSubpath doesn't have an end!\nBreaking the loop -> ",subpath)
        return [],[]
    
    return compute_binary_and_startend_list(first_most_start,last_most_end,subpath_loc_list)

# In[20]:

def get_binary_and_startend_list_for_all_subpaths(subpath_list):
    
    binary_subpath_list = []
    start_end_nodes_subpath_list = []
    
    for subpath in subpath_list:
        binary_node_list,start_end_node_list = get_subpath_info(subpath)
        binary_subpath_list.append(binary_node_list)
        start_end_nodes_subpath_list.append(start_end_node_list)
        
    return binary_subpath_list, start_end_nodes_subpath_list



def filter_subpath_list_by_subexon_ids(subpath_list, subpath_freq_list, subexon_ids):
    indices = []
    for i in range(len(subpath_list)):
        if set(subpath_list[i]) <= set(subexon_ids):
            indices.append(i)
    subset_subpath_list = [subpath_list[idx] for idx in indices]
    subset_subpath_freq_list = [subpath_freq_list[idx] for idx in indices]
    return subset_subpath_list,subset_subpath_freq_list

