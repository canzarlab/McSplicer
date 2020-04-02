#!/usr/bin/env python

"""parse_cnt_file.py: Parses cnt file."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


import os
from pandas import *

def parse_cnt_file(filename, paired_end = False):
    
    """ Parses .cnt file
    Returns two lists:
        subpath_list ->   list of read spans, each read span contains a list node ids in that read
        subpath_freq_list -> the frequency of each read
    """
    
    if not paired_end:
        
        with open(filename) as f:
            cnt_data = f.readlines()

        subpath_list = []
        subpath_freq_list = []

        for line in cnt_data:
            if line.startswith('#'):
                continue
            entry = line.split()
            subpath_list.append([int(i) for i in entry[0].split('-')])
            subpath_freq_list.append(float(entry[1]))
        
        return subpath_list, subpath_freq_list
    
    else:
        cnt_data = read_csv(filename, delim_whitespace=True, header=None)
        subpath_l_list = []
        subpath_r_list = []
        subpath_freq_list = []
        for _,entry in cnt_data.iterrows():

            s_l = entry[0].split('^')[0] # Left suppath
            s_r = entry[0].split('^')[1] # Right ssubpath
            subpath_l_list.append([int(i) for i in s_l.split('-')])
            subpath_r_list.append([int(i) for i in s_r.split('-')])
            subpath_freq_list.append(float(entry[1]))

        return subpath_l_list,subpath_r_list, subpath_freq_list

