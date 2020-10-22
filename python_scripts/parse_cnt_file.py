#!/usr/bin/env python

"""parse_cnt_file.py: Parses cnt file."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


import os
from collections import defaultdict

def parse_cnt_file(filename, paired_end = False):
    
    """ Parses .cnt file
    Returns two lists:
        subpath_list ->   list of read spans, each read span contains a list node ids in that read
        subpath_freq_list -> the frequency of each read
    """
    
    with open(filename) as f:
        cnt_data = f.readlines()


    subpath_list = []
    subpath_freq_list = []
    
    subpath_dict = defaultdict(int)
    
    for line in cnt_data:
        
            
        if line.startswith('#'):
            continue
        
        if line.find('^') == -1:
            # single end read
            entry = line.split()
            subpath_dict[entry[0]] += float(entry[1])

        else:
            # Left and right suppaths
            entry = line.split()
            s_l,s_r = entry[0].split('^')             
            read_count = float(entry[1])
            subpath_dict[s_l] += read_count
            subpath_dict[s_r] += read_count
            

            
    for key in subpath_dict:
        subpath = [int(i) for i in key.split('-')]
        subpath_list.append(subpath)
        subpath_freq_list.append(subpath_dict[key])
            
            
        
    return subpath_list, subpath_freq_list

