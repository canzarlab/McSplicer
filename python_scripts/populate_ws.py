
# coding: utf-8

# In[ ]:



#!/usr/bin/env python

"""compute_ws.py: See section 3.4 in Alternative Splicing pdf."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"



    

def populate_ws(q_arr,p_arr,binary_subpath_list,start_end_nodes_subpath_list,subpath_list):
    """ 
    w(s): denotes the conditional probability of S_n=s conditional on X_F(s) is 
    a part of an isoform
    """ 
    
    w = []
    
    for subpath_index in range(0,len(subpath_list)): 
       
        subpath = subpath_list[subpath_index]
        
        binary_node_list = binary_subpath_list[subpath_index]
        
        
        start_end_node_list = start_end_nodes_subpath_list[subpath_index]
        
        probability = 1.
        # If only one segment then probability is just 1.
        if len(binary_node_list) > 1:
            trimmed_start_end_node_list = start_end_node_list[1:-1]
            
            
            for i in range(len(binary_node_list)-1):
                segment_1 = binary_node_list[i]
                segment_2 = binary_node_list[i+1]
                
     
                if i >= len(trimmed_start_end_node_list):
                    break
                
                type_site = trimmed_start_end_node_list[i][0]
                site_index = int(trimmed_start_end_node_list[i][1:])

                if type_site == 's':
                    if segment_1==0 and segment_2==0:
                        probability*=(1-p_arr[site_index])
                    elif segment_1==0 and segment_2==1:
                        probability*=p_arr[site_index]
                    elif segment_1==1 and segment_2==1:
                        probability*=1
                    else:
                        probability = 1e-10
                        #print 'compute_ws.py  -> subpath',subpath,i,i+1
                        print("Error3: \nfile -> compute_ws.py \nPair of 1-0 segments is seperated by start site")    

                if type_site == 'e':
                    if segment_1==1 and segment_2==1:
                        probability*=(1-q_arr[site_index])
                    elif segment_1==1 and segment_2==0:
                        probability*=q_arr[site_index]
                    elif segment_1==0 and segment_2==0:
                        probability*=1
                    else:
                        probability = 1e-10
                        print("Error4: \nfile -> compute_ws.py \nPair of 0-1 segments seperated by end site")
        
        w.append(probability)

    return w

