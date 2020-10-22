
# coding: utf-8


#!/usr/bin/env python

"""global_vars.py: This file contains all the global variables."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


import os
import sys
import numpy as np
from parse_gtf_file import *
from parse_cnt_file import *



'''
This file contains the following variables:
    strand_dir: forward (+) or backward (-)
    loc_index_dict: a dict of locations (keys) and their indices (values)
    start_sites_dict: start locations (keys) and their indices (values)
    end_sites_dict: end locations (keys) and their indices (values)
    subexon_ids_dict: subexon_ids (keys) and their coordinates [start, end] as values
    end_sites_list: sorted list of end locations
    start_sites_list: sorted list of start locations
    loc_list: sorted list of all locations
'''



DIR = '../data/'


gtf_file_list = ['altacc.refined.gtf','altdon.refined.gtf','alttes.refined.gtf','alttss.refined.gtf','exskip.refined.gtf']
#gtf_filename = gtf_file_list[0]

gtf_file_strn_mlck = 'Strn_Mlck/Drosophila_melanogaster.BDGP6.89.chr.transid.refined_Strn_Mlck_gene.gtf'
gtf_file_sls = 'sls/Drosophila_melanogaster.BDGP6.89.chr.transid.refined_sls_gene.gtf'
gtf_file_wupA = 'wupA/Drosophila_melanogaster.BDGP6.89.chr.transid.refined.wupA_gene.gtf'
gtf_file_up = 'up/Drosophila_melanogaster.BDGP6.89.chr.transid.refined_up_gene.gtf'
gtf_file_Mhc = 'Mhc/Drosophila_melanogaster.BDGP6.89.chr.transid.refined_Mhc_gene.gtf'
gtf_file_pfk = 'pfk/Drosophila_melanogaster.BDGP6.89.chr.transid.refined_pfk_gene.gtf'

gtf_filename = gtf_file_list[4]


strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict,subexon_trans_ids_dict = parse_gtf_file_general(DIR+gtf_filename)

if strand_dir == '+' :
    end_sites_list = sorted(end_sites_dict.keys())
    start_sites_list = sorted(start_sites_dict.keys())
    loc_list = sorted(loc_index_dict.keys())
else:
    end_sites_list = sorted(list(end_sites_dict.keys()), reverse=True)
    start_sites_list = sorted(list(start_sites_dict.keys()), reverse=True)
    loc_list = sorted(list(loc_index_dict.keys()), reverse=True)
