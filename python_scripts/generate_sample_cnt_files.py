
# coding: utf-8
#!/usr/bin/env python

"""generate_sample_cnt_files.py: generate sample cnt files from the following gtf files.
    1. alttss.refined.gtf
    2. altacc.refined.gtf
    3. altdon.refined.gtf
    4. alttes.refined.gtf
    5. exskip.refined.gtf
"""
__authors__      = "Yash Kumar Sonthalia, Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


import os
import random
from parse_gtf_file import *
from collections import Counter


def sample_reads_from_transcript(read_length,transcript,no_samples,subexon_ids_dict,trans_length):
    '''
    Parameters:
        read_length: length of read (default 75)
        transcript: a transcript consists of a list of subexon ids forming that transcript
        no_samples: number of read samples to generate from this transcript
        subexon_ids_dict: a dict of subexon_ids and their corresponding start and end sites
        trans_length: transcript length
        
    Returns:
        read_freq_dict: read span over subexon(s) and the number of times we sample a read from the same subexon(s)
    
    '''
    
    read_freq_dict = {}
    for i in range(no_samples):
        read_start = random.randint(0, trans_length) # each read start at a random location and spans over a single or multiple of subexons
        #read_start = 15
        #print "read_start: ",read_start
        read_left = read_length   # default read length is 75
        read_found = False
        read_span = []            # a list of subexon ids that the current read spreads over
        
        for subexon_id in transcript:

            if read_left <= 0:
                break # stop because we alread read a span of length read_length
                
            subexon_coordinates = subexon_ids_dict[subexon_id]
            subexon_length = subexon_coordinates[1] - subexon_coordinates[0]
            
            # The start of this read exists within this subexon
            if  read_start < subexon_length or read_found:
                read_found = True
                read_left -= (subexon_length - read_start)
                read_start = 0                 # start the read from the begining of the next subexon if read_left > 0
                read_span.append(str(subexon_id))  #append the current subexon to the read span
                
            elif read_start > subexon_length:
                read_start -= subexon_length
                
        #before_read_span = read_span
        
        read_span = '-'.join(read_span)
        if read_span == '':
            continue
            # sth went wrong while reading read counts, and the read is not assigned to any exon/subexon, skip.
            #print '>>>>',read_span
            #print '>>>>',before_read_span
        
        if read_span not in read_freq_dict:
            read_freq_dict[read_span] = 1.
        else:
            read_freq_dict[read_span] += 1.
            
    return read_freq_dict


'''
 This function generates cnt files from a single gtf file
''' 
def create_cnt_files(gtf_file,total_num_files,total_samples,in_DIR,out_DIR):
    
    for idx in range(0,total_num_files):

        out_file = out_DIR+gtf_file.split('.')[0]+'_'+str(idx)+'.cnt'

        strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict,trans_id_subexons_dict = parse_gtf_file_general(in_DIR+gtf_file)

        trans_subexon_list = list(trans_id_subexons_dict.values()) # list of subexon ids for each transcript
        trans_subexon_list.sort()

        # Initial vars
        #random.seed(0) #to get the same results everytime we run the code
        random.seed(idx) # to get different results

        no_trans = len(trans_subexon_list)
        read_length = 75
        a_prop = [1./no_trans] * no_trans
        N_list = [0.] * no_trans
        

        lengths_trans = []
        for subexon_list in trans_subexon_list:
            length = -read_length+1   # we read up to transcript_length - read_length + 1
            # the length of the transcript is the total length of its subexons
            for subexon in subexon_list:
                subexon_start_end_locs = subexon_ids_dict[subexon]
                length+=subexon_start_end_locs[1] - subexon_start_end_locs[0] 
            lengths_trans.append(length)


        denum = 0.
        for i in range(no_trans):
            N_list[i] = (lengths_trans[i]*a_prop[i])
            denum+= N_list[i]

        N_list = [i/denum for i in N_list]

        #print N_list
        #print subexon_ids_dict
        #print "trans_subexon_list>> ",trans_subexon_list
        #print lengths_trans

        read_freq_list = []
        for i in range(no_trans):
            num_samples = int(total_samples * N_list[i]) # total number of reads from that transcript
            #print trans_subexon_list[i]
            _read_dict = sample_reads_from_transcript(75,trans_subexon_list[i],num_samples,subexon_ids_dict,lengths_trans[i])
            read_freq_list.append(_read_dict)

        read_freq_dict = dict(sum((Counter(dict(x)) for x in read_freq_list),Counter()))
        #print 'read_freq_dict>>>',read_freq_dict

        ## write to cnt file
        file = open(out_file,'w')  
        for key,val in read_freq_dict.items():    
            file.write(key+' '+str(val)+'\n') 
        file.close() 


def main():
    
    in_DIR = '../data/'
    gtf_file_list = ['alttss.refined.gtf','altacc.refined.gtf','altdon.refined.gtf','alttes.refined.gtf','exskip.refined.gtf']  
    total_num_files = 100   # get 100 files for each case
    total_samples = 1000    # total number of reads
    
    for gtf_file in gtf_file_list:
        
        print('Generating count files for',gtf_file.split('.')[0].upper())
        
        out_DIR = in_DIR+gtf_file.split('.')[0].upper()+'/'
        
        if not os.path.isdir(out_DIR):
            os.system('mkdir '+out_DIR)
        
        create_cnt_files(gtf_file,total_num_files,total_samples,in_DIR,out_DIR)
        
        print('Output is written to',out_DIR)


main()

