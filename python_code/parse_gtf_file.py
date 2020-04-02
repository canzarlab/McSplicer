
# coding: utf-8

# In[1]:

#!/usr/bin/env python

"""parse_gtf_file.py: Parses gtf file."""

__authos__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"

import os
import csv
from pandas import *
import random
import math


# In[2]:

def create_location_dicts(start_sites, end_sites, strand_dir):
    """ This func expects sorted start_sites and end_sites based on strand direction
    Returns 3 dicts:
    loc_index_dict ->   Key: location, val: index
    start_sites_dict -> key: start location, index (helps to determine s1, s2, ...)
    end_sites_dict ->   key: end location, index (helps to determine e1, e2, ...)
    """
    
    loc_index_dict = {}       
    start_sites_dict = {}     
    end_sites_dict = {}       
    
    location_list = []      # List of all start and end locations
    location_list.extend(start_sites)
    location_list.extend(end_sites)

    if strand_dir == '+':
        location_list.sort()
    else:
        location_list.sort(reverse=True)

    index = 0
    for location in location_list:
        loc_index_dict[location] = index
        index+=1

    index = 0
    for location in start_sites:
        start_sites_dict[location] = index
        index+=1

    index = 0
    for location in end_sites:
        end_sites_dict[location] = index
        index+=1
        
    return loc_index_dict, start_sites_dict, end_sites_dict




def parse_gtf_file_general(filename):
    """ Parses .gtf file
    Returns three dict:
        loc_index_dict ->   Key: location, val: index
        start_sites_dict -> key: start location, val: index (helps to determine s1, s2, ...)
        end_sites_dict ->   key: end location, val: index (helps to determine e1, e2, ...)
        subexon_ids_dict -> key: subexon_id or node_id, val: [subexon_start_siye, subexon_end_site]
    """
    subexon_data = read_csv(filename, sep=';',header=None,names='abcdefghi')
    subexon_start_sites = []
    subexon_end_sites = []
    splice_ends = []
    splice_directions = []
    subexon_ids = []

    subexon_ids_dict = {}
    subexon_trans_ids_dict = {}

    for _,row in subexon_data.iterrows():
        #print row
        entry0 = row[0].split('\t') # 1st entry in the row contains subexon coordinates, dir and type of row.
        if  entry0[2] != 'subexon':
            continue
    
        #entry3 = row[3]   # 3rd entry to determine whether a subexon is R, L, or B for Right, Left, or Both respectively 
        #entry4 = row[4]   # 4th entry to determine if the subexon belongs to forward (+) or backward (-) strand
        #entry5 = row[5]   # 5th entry in the row has the node id, i.e., subexon id

        strand_dir = entry0[6]
        splice_directions.append(strand_dir)
        
        start_site = int(entry0[3])
        end_site = int(entry0[4])
        if strand_dir == '+':
            subexon_start_sites.append(start_site)
            subexon_end_sites.append(end_site)
        else:
            subexon_start_sites.append(end_site)
            subexon_end_sites.append(start_site)

        for col in row[1:]:
            if isinstance(col, float) and math.isnan(col):
                continue
            #print col
            tag = col.split(' ')[1]
            value = col.split(' ')[2]
            if tag=='SpliceEnd':
                splice_ends.append(value[1:2]) # R, L, or B
            elif tag=='NodeId':
                subexon_id = int(value)
                subexon_ids.append(subexon_id)
                if strand_dir == '+':
                    subexon_ids_dict[subexon_id] = [start_site,end_site]
                else:
                    subexon_ids_dict[subexon_id] = [end_site,start_site]
            elif tag=='transcript_id':
                trans_id = value[1:-1]
                if trans_id in subexon_trans_ids_dict:
                    subexon_trans_ids_dict[trans_id].append(subexon_id)
                else:
                    subexon_trans_ids_dict[trans_id] = [subexon_id]

    #print subexon_start_sites
    #print subexon_end_sites
    #print splice_ends
    #print splice_directions
    """
    Forward strand (+):
       potential start sites ->  Left of L 
                                 Left of B 
       potential end sites ->    Right of R
                                 Right of B
                                 
            s_____s____ e_____e
            |__L__|__B__|__R__|
 
 
    Backward strand (+):
        potential start site -> right of L
                                right of B
        potential end site ->  Left of R
                               Left of B

            e_____e_____s_____s
            |__R__|__B__|__L__|

    """ 
    start_sites = []   # start sites s1, s2, s3, etc.
    end_sites = []     # end sites e1, e2, e3, etc.

    strand_dir = splice_directions[0] # strand direction is the same in the whole file
    
    if strand_dir == '+':
        forward_strand = True
    elif strand_dir == '-':
        forward_strand = False 
    else:
        print 'Error: Splice direction must be + or -. Undefined splice direction  -> ' + splice_ends[i]

    for i in range(len(subexon_start_sites)):

        if forward_strand == False:      
            if splice_ends[i]=='R':
                end_sites.append(subexon_end_sites[i])

            elif splice_ends[i]=='L':
                 start_sites.append(subexon_start_sites[i])

            elif splice_ends[i]=='B':
                start_sites.append(subexon_start_sites[i])
                end_sites.append(subexon_end_sites[i])
            
            elif splice_ends[i] !='-': # dash means internal exon, just ignore it, otherwise show error
                print 'Error: Splice end value must be L, R, or B. Undefined splice end -> ' + splice_ends[i]

        elif forward_strand == True:   
            if splice_ends[i]=='R':
                end_sites.append(subexon_end_sites[i])

            elif splice_ends[i]=='L':
                 start_sites.append(subexon_start_sites[i])

            elif splice_ends[i]=='B':
                start_sites.append(subexon_start_sites[i])
                end_sites.append(subexon_end_sites[i])
                
            elif splice_ends[i] !='-': # dash means internal exon, just ignore it, otherwise show error
                print 'Error: Splice end value must be L, R, or B. Undefined splice end -> ' + splice_ends[i]

    if forward_strand == False:
        end_sites.sort(reverse=True)
        start_sites.sort(reverse=True)
    else:
        end_sites.sort()
        start_sites.sort()

    #m_s = len(start_sites)
    #m_e = len(end_sites)
    #print end_sites
    #print m_s
    #print m_e
    loc_index_dict, start_sites_dict, end_sites_dict = create_location_dicts(list(unique(start_sites)), list(unique(end_sites)), strand_dir)
    
    return strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict, subexon_trans_ids_dict



def parse_gtf_file_per_gene(filename, gene_id):
    """ Parses .gtf file
    Returns three dict:
        loc_index_dict ->   Key: location, val: index
        start_sites_dict -> key: start location, val: index (helps to determine s1, s2, ...)
        end_sites_dict ->   key: end location, val: index (helps to determine e1, e2, ...)
        subexon_ids_dict -> key: subexon_id or node_id, val: [subexon_start_siye, subexon_end_site]
    """
    subexon_data = read_csv(filename, sep=';',header=None,names='abcdefghi',low_memory=False) ## I added low_memory option to avoid getting an annoying warning!
    subexon_start_sites = []
    subexon_end_sites = []
    splice_ends = []
    splice_directions = []
    subexon_ids = []

    subexon_ids_dict = {}
    subexon_trans_ids_dict = {}
    
    strand_dir = ''
    gene_found = False
    loc_index_dict = {}
    start_sites_dict = {}
    end_sites_dict = {}

    for _,row in subexon_data.iterrows():
        entry0 = row[0].split('\t') # 1st entry in the row contains subexon coordinates, dir and type of row.
        if  entry0[2] != 'subexon':
            continue
        
         
        for idx in range(len(row)):
            item = row[idx].split(' ')
            if item[1] == 'gene_id':
                this_gene_id =  item[2][1:len(item[2])-1] #remove double qouts
                break
                

        if this_gene_id != gene_id:
            continue
        else:
            gene_found = True
    
        #entry3 = row[3]   # 3rd entry to determine whether a subexon is R, L, or B for Right, Left, or Both respectively 
        #entry4 = row[4]   # 4th entry to determine if the subexon belongs to forward (+) or backward (-) strand
        #entry5 = row[5]   # 5th entry in the row has the node id, i.e., subexon id

        strand_dir = entry0[6]
        splice_directions.append(strand_dir)
        
        start_site = int(entry0[3])
        end_site = int(entry0[4])
        if strand_dir == '+':
            subexon_start_sites.append(start_site)
            subexon_end_sites.append(end_site)
        else:
            subexon_start_sites.append(end_site)
            subexon_end_sites.append(start_site)

        # Parsing a list of tag-value pairs
        for col in row[1:]:
            if isinstance(col, float) and math.isnan(col):
                continue
            #print col
            tag = col.split(' ')[1]
            value = col.split(' ')[2]
            
            #print tag,value
            
            if tag=='SpliceEnd':
                splice_ends.append(value[1:2]) # R, L, or B
            elif tag=='NodeId':
                subexon_id = int(value)
                subexon_ids.append(subexon_id)
                if strand_dir == '+':
                    subexon_ids_dict[subexon_id] = [start_site,end_site]
                else:
                    subexon_ids_dict[subexon_id] = [end_site,start_site]
            elif tag=='transcript_id':
                trans_id = value[1:-1]
                if trans_id in subexon_trans_ids_dict:
                    subexon_trans_ids_dict[trans_id].append(subexon_id)
                else:
                    subexon_trans_ids_dict[trans_id] = [subexon_id]
            #elif tag=='gene_id':
             #   print 'gene_id',value

    #print subexon_start_sites
    #print subexon_end_sites
    #print splice_ends
    #print splice_directions
    """
    Forward strand (+):
       potential start sites ->  Left of L 
                                 Left of B 
       potential end sites ->    Right of R
                                 Right of B
                                 
            s_____s____ e_____e
            |__L__|__B__|__R__|
 
 
    Backward strand (+):
        potential start site -> right of L
                                right of B
        potential end site ->  Left of R
                               Left of B

            e_____e_____s_____s
            |__R__|__B__|__L__|

    """ 
    
    if not gene_found: # if gene not found
        print 'Error: No values returned.\nThe gene ID', gene_id, 'was not found in the gtf file,', filename
        
    else:
    
        start_sites = []   # start sites s1, s2, s3, etc.
        end_sites = []     # end sites e1, e2, e3, etc.

        strand_dir = splice_directions[0] # strand direction is the same in the whole file

        if strand_dir == '+':
            forward_strand = True
        elif strand_dir == '-':
            forward_strand = False 
        else:
            print 'Error: Splice direction must be + or -. Undefined splice direction  -> ' + splice_ends[i]

        for i in range(len(subexon_start_sites)):

            if forward_strand == False:      
                if splice_ends[i]=='R':
                    end_sites.append(subexon_end_sites[i])

                elif splice_ends[i]=='L':
                     start_sites.append(subexon_start_sites[i])

                elif splice_ends[i]=='B':
                    start_sites.append(subexon_start_sites[i])
                    end_sites.append(subexon_end_sites[i])

                elif splice_ends[i] != '-': # - means internal node, just ignore it, otherwise show an error
                    print 'Error: Splice end value must be L, R, or B. Undefined splice end -> ' + splice_ends[i]

            elif forward_strand == True:   
                if splice_ends[i]=='R':
                    end_sites.append(subexon_end_sites[i])

                elif splice_ends[i]=='L':
                     start_sites.append(subexon_start_sites[i])

                elif splice_ends[i]=='B':
                    start_sites.append(subexon_start_sites[i])
                    end_sites.append(subexon_end_sites[i])
                
                elif splice_ends[i] != '-':
                    print 'Error: Splice end value must be L, R, or B. Undefined splice end -> ' + splice_ends[i]

        if forward_strand == False:
            end_sites.sort(reverse=True)
            start_sites.sort(reverse=True)
        else:
            end_sites.sort()
            start_sites.sort()

        #m_s = len(start_sites)
        #m_e = len(end_sites)
        #print end_sites
        #print m_s
        #print m_e
        loc_index_dict, start_sites_dict, end_sites_dict = create_location_dicts(list(unique(start_sites)), list(unique(end_sites)), strand_dir)
    
    return strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict, subexon_trans_ids_dict


def get_list_of_genes_from_gtf_file(gtf_file):
    subexon_data = read_csv(gtf_file, sep=';',header=None,names='abcdefghi')
    gene_id_list = []
    for row,_ in subexon_data.iterrows():
        entry0 = row[0].split('\t') # 1st entry in the row contains subexon coordinates, dir and type of row.
        if  entry0[2] != 'subexon':
            continue


        for idx in range(len(row)):
            item = row[idx].split(' ')
            if item[1] == 'gene_id':
                gene_id =  item[2][1:len(item[2])-1] #remove double qouts
                gene_id_list.append(gene_id)
                break
    
    return unique(gene_id_list)


