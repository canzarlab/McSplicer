'''
 This code analyzes the junction reads and outputs green, black, red counts as follows:
 if a junction read spans over an end site, it can be:
    ------> a. red: if a junction read spans over the subexon of the current end site and any later non-contiguous subexon
    ------> b. black: if a junction read spans over the subexon of the current end site and the next contiguous neighbouring subexon
 if a junction read spans over a start site, it can be:
    ------> c. green: if a junction read  spans over the subexon of the current start site and a previous subexon 
    ------> but NOT the one that is exactly preceding to the subexon in question
'''



__author__      = "Israa Alqassem"
__copyright__   = "Copyright 2018, McSplicer"


import os,sys,csv
import argparse
from new_gtf_genome_wide_parser import *
from parse_cnt_file import parse_cnt_file
from subpath_info import *


def get_ss_usage_from_read_counts(gene_datalist):

    strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict = get_gene_data(gene_id,gene_datalist)

    reverse_val = False
    if strand_dir == '-':
        reverse_val = True

    end_sites_list = sorted(list(end_sites_dict.keys()), reverse=reverse_val)
    start_sites_list = sorted(list(start_sites_dict.keys()), reverse=reverse_val)
    loc_list = sorted(list(loc_index_dict.keys()), reverse=reverse_val)


    #subexon_ids_dict
    loc_subexon_dict = {}
    for subexon_id in subexon_ids_dict:
        start,end = subexon_ids_dict[subexon_id]
        #print subexon_id,start,end
        loc_subexon_dict[start] = subexon_id
        loc_subexon_dict[end] = subexon_id


    output = []
    output.append(['Location','Index','Green count','Red count','Black_count'])

    for site in loc_list:

        subexon = loc_subexon_dict[site]
        this_start, this_end = subexon_ids_dict[subexon]

        if site in end_sites_dict:
            key = 'e'+str(end_sites_dict[site])
            idx = -1 # the last elemet, the path ends at the end of this subexon
        elif site in start_sites_dict:
            key = 's'+str(start_sites_dict[site])
            idx = 0 # the first elemet, the path starts at the start of this subexon
        else:
            continue

        black_count = 0.
        red_count = 0.
        green_count = 0.

        if key == 's0':
            output.append([site,key,green_count,red_count,black_count])
            continue

        #print '**************',key,site,subexon,idx,'**************'


        for i in range(len(subpath_list)):        
            if len(subpath_list[i]) > 1 and subexon in subpath_list[i]:

                # Index of this subexon in the list
                idx_subexon = subpath_list[i].index(subexon)

                #print subpath_list[i],subexon,idx_subexon

                if key[0:1]=='e':

                    next_subexon_idx = idx_subexon + 1 # id in the list


                    if len(subpath_list[i]) > next_subexon_idx:
                        next_subexon = subpath_list[i][next_subexon_idx]
                        next_start, next_end = subexon_ids_dict[next_subexon]
                        #print 'exons boundaries:', this_start,this_end, next_start, next_end
                        
                        if strand_dir == '+' :
                            if (this_end == next_start or this_end == next_start-1):
                                black_count += subpath_freq_list[i]
                            elif next_start > this_end+1:
                                red_count += subpath_freq_list[i]
                                
                        elif strand_dir == '-':
                            if (this_end == next_start or this_end == next_start+1):
                                black_count += subpath_freq_list[i]
                            elif next_start < this_end-1:
                                red_count += subpath_freq_list[i]
                            

                    elif  subpath_list[i][idx_subexon] == subpath_list[i][-1] and key[0:1]=='e':
                        #print 'read ends at current subexon'
                        continue


                elif key[0:1]=='s' and idx_subexon != 0:
                    prev_subexon_idx = idx_subexon - 1       # id in the list
                    prev_subexon = subpath_list[i][prev_subexon_idx]
                    prev_start, prev_end = subexon_ids_dict[prev_subexon]
                    #print 'exons boundaries:', prev_start, prev_end, this_start,this_end
                    if strand_dir == '+' :
                        if this_start != prev_end  and this_start != prev_end+1:
                            green_count += subpath_freq_list[i]

                    elif strand_dir == '-':
                        if this_start != prev_end  and this_start != prev_end-1:
                            green_count += subpath_freq_list[i]


        output.append([site,key,green_count,red_count,black_count])               
        #print "red_count=",red_count
        #print "black_count=",black_count
        #print "green_count=",green_count
        

    return output

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf_file', type=str, help='Input refined gtf file name.', required=True)
    parser.add_argument('--cnt_file', type=str, help='Input count file name.', required=True)
    parser.add_argument('--gene_id', type=str, help='Input gene ID.', required=True)
    parser.add_argument('--out_dir', type=str, required=True,help='Destination for output file(s).')
    args = parser.parse_args()
    
    
    return args



if __name__ == "__main__":

    params = vars(get_args())

    refined_gtf = params['gtf_file']
    cnt_file = params['cnt_file']
    gene_id = params['gene_id']
    out_DIR = params['out_dir'] + '/'
    
    #refined_gtf = '../data/alttes_b.refined.gtf' #'../data/exskip.refined.gtf'
    #cnt_file = '../data/ALTTES_b/alttes_b_0.cnt'   #'../data/EXSKIP/exskip_0.cnt' 
    #gene_id = 'A'
    #out_DIR = '/home/isra/Desktop/splice_site_usage/'
    

    out_filename = out_DIR+gene_id+'.csv'

    all_gene_dict = get_all_genes_dict(refined_gtf)
    all_subpath_list, all_subpath_freq_list = parse_cnt_file(cnt_file)


    gene_datalist = all_gene_dict[gene_id]
    subexon_ids = [row[0] for row in gene_datalist]
    subpath_list, subpath_freq_list = filter_subpath_list_by_subexon_ids(all_subpath_list, all_subpath_freq_list, subexon_ids)

    output = get_ss_usage_from_read_counts(gene_datalist)

    f = open(out_filename, 'w')
    writer = csv.writer(f)


    for row in output:
         writer.writerow(row)
    f.close()

    print('\nOutput is written to the folder:\n\t',out_filename)

