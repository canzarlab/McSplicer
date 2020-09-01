
# coding: utf-8

# <h3>Modified version</h3>
# 
# <ul>
#     <li>Added start and stop nodes to resolve the identifiability issue</li>
#     <li>End sites are computed from junction reads, not in EM</li>
# </ul>

# In[1]:




import argparse
import os,sys,math,operator
import os.path
from os import listdir
from collections import defaultdict
from pandas import *
import numpy as np
from subpath_info import *
from compute_f00_f01_f10_f11 import *
from populate_PZ import *
from populate_ws import *
from expected_prefix_suffix_length import *
from cases_code_new import *
from parse_cnt_file import parse_cnt_file
from new_gtf_genome_wide_parser import *
from add_start_end_read_counts import *
np.random.seed(42)


# In[2]:


def get_gene_tx_dict(abundance_file):
    
    gene_tx_dict = {}
    header_line = True
    
    with open(abundance_file,'rb') as in_file:
        for line in in_file:
            if header_line:
                header_line = False
                continue
            line = line.split() 

            tx_id = line[2][1:-1]
            gene_id = line[1][1:-1]
            counts = float(line[12])

            if gene_id not in gene_tx_dict:
                gene_tx_dict[gene_id] = []
            gene_tx_dict[gene_id].append([tx_id,counts])
            
    return gene_tx_dict


# In[3]:


def filter_out_unexpressed_transcripts(gene_datalist,gene_tx_list):
    
    unexpressed_tx_count = 0
    total_tx_count = 0
    gene_data = []
    

    subexon_trans_ids_dict = defaultdict(list)
    for row in gene_datalist:
        subexon_id = row[0]
        trans_id = row[5]
        subexon_trans_ids_dict[trans_id].append(subexon_id)

    unexpressed_tx = [tx for tx,cnt in gene_tx_list if cnt <=10e-6]
    unexpressed_tx_count += len(unexpressed_tx)
    total_tx_count += len(gene_tx_list)


    for row in gene_datalist:
        tx_id = row[-1]
        if tx_id not in unexpressed_tx:
            gene_data.append(row)
            
            
    return gene_data


# In[4]:


#python2 python_code/bootstrap.py --gtf_file /data/israa/RNASeqMixology/GRCh37_latest_genomic.refined.gtf --cnt_file /data/israa/RNASeqMixology/new_cnt_hg19/R1-050_C5N4YANXX_TTAGGCAT.cnt --gene_id gene17406 --out_dir /data/israa/RNASeqMixology/output2/ --no_steps 10 --read_length 100
def bootstrap(subpath_freq_list):

    #print '(Before) Freq list',subpath_freq_list
    
    total = sum(subpath_freq_list)
    #print 'total=',total
    probabilities = [float(x)/total for x in subpath_freq_list]
    # Generate random samples drawn from the same probability distributions
    samples = np.random.choice(len(probabilities), int(total), p=probabilities)
    subpath_freq_list = [0.] * len(probabilities)
    for x in samples:
        subpath_freq_list[int(x)] += 1

    #print '(After) Freq list',subpath_freq_list
    return subpath_freq_list


# In[5]:


def get_q_vector_from_junction_reads(strand_dir,subexon_ids_dict,loc_list,end_sites_dict,subpath_list,subpath_freq_list):

    
    loc_subexon_dict = {}
    for subexon_id in subexon_ids_dict:
        start,end = subexon_ids_dict[subexon_id]
        loc_subexon_dict[start] = subexon_id
        loc_subexon_dict[end] = subexon_id
      

    end_sites_arr = []
    #end_sites_arr.append(['Location','Index','Probabilty','Red count','Black_count'])
    
    for site in loc_list:

        subexon = loc_subexon_dict[site]
        this_start, this_end = subexon_ids_dict[subexon]

        if site in end_sites_dict:
            key = 'e'+str(end_sites_dict[site])
            idx = -1 # the last elemet, the path ends at the end of this subexon
        else:
            continue

        black_count = 0.
        red_count = 0.

        for i in range(len(subpath_list)):        
            if len(subpath_list[i]) > 1 and subexon in subpath_list[i]:
                # Index of this subexon in the list
                idx_subexon = subpath_list[i].index(subexon)
                
                if key[0:1]=='e':
                    next_subexon_idx = idx_subexon + 1 # id in the list


                    if len(subpath_list[i]) > next_subexon_idx:
                        next_subexon = subpath_list[i][next_subexon_idx]
                        next_start, next_end = subexon_ids_dict[next_subexon]
                        
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


        if key[0:1]!='s':
            prob_val = 0
            if red_count+black_count > 0:
                prob_val = red_count/(red_count+black_count)

            prob_val = 1e-6 if prob_val == 0. else prob_val
            prob_val = 1.0-1e-6 if prob_val == 1. else prob_val
            end_sites_arr.append([site,key,prob_val,red_count,black_count])               
    
    end_sites_arr[-1][2] = 1. # The prob of last end site is 1
    
    return end_sites_arr


# In[6]:


def run_EM_bootstrap( cnt_file, gene_id, no_steps, read_length,gene_datalist,disable_end_sites_computation,add_start_stop_nodes):
    
    ######################  Step 1: Read gtf file
    
    strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict = get_gene_data(gene_id,gene_datalist)
    
    
    subexon_trans_ids_dict = {}
    for row in gene_datalist:
        subexon_id = row[0]
        trans_id = row[5]
        if trans_id not in subexon_trans_ids_dict:
            subexon_trans_ids_dict[trans_id] = []
        subexon_trans_ids_dict[trans_id].append(subexon_id)
    
    if strand_dir == '':
        return {},{},[],[]
    
    if strand_dir == '+' :
        end_sites_list = sorted(end_sites_dict.keys())
        start_sites_list = sorted(start_sites_dict.keys())
        loc_list = sorted(loc_index_dict.keys())
    else:
        end_sites_list = sorted(end_sites_dict.keys(), reverse=True)
        start_sites_list = sorted(start_sites_dict.keys(), reverse=True)
        loc_list = sorted(loc_index_dict.keys(), reverse=True)
        
    
    ######################  Step 2: Read cnt file and get all subpaths list that belong to this gene
    p_arr_list = []
    q_arr_list = []
    subpath_list, subpath_freq_list = parse_cnt_file(cnt_file)
    subpath_list, subpath_freq_list = filter_subpath_list_by_subexon_ids(subpath_list, subpath_freq_list, subexon_ids_dict.keys())
   

#    print 'subpath_list:',subpath_list
#    print 'subpath_freq_list:',subpath_freq_list

 
    if len(subpath_list) < 1:
        print 'No counts are found for gene', gene_id, '. Continuing with the next gene.'
        return {},{},[],[]
    
            
    if add_start_stop_nodes:
        
        ###################### Step 2.a: add fake synthetic counts to the start and end nodes
        subpath_list, subpath_freq_list = preprocess_read_counts(subexon_ids_dict,subexon_trans_ids_dict,read_length,subpath_list,subpath_freq_list)

        ###################### Step 2.b: Add starting and stopping nodes
        subexon_ids_dict,start_sites_dict,end_sites_dict,loc_index_dict,start_sites_list,end_sites_list,loc_list = add_starting_and_stopping_nodes(strand_dir,subexon_ids_dict,subexon_trans_ids_dict,start_sites_list,end_sites_list,loc_list,start_sites_dict,end_sites_dict,loc_index_dict)

        
    
    orig_subpath_freq_list = subpath_freq_list
    
    for curr_step in range(no_steps+1): # Total bootstrap steps
        
        ######################  Step 3: Set global varaibales in other files
        set_global_vars_in_cases_builders(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)
        set_global_vars_in_subpath(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)
        set_global_vars_in_prefix_suffix(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)
        set_global_vars_in_compute_fs(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)
        set_global_vars_in_populate_PZ(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,subexon_trans_ids_dict,end_sites_list,start_sites_list,loc_list)
        
        
        ###################### Step 4: Initialize variables for EM 
        binary_subpath_list, start_end_nodes_subpath_list = get_binary_and_startend_list_for_all_subpaths(subpath_list)

        all_start_end_list = []
        for i in range(len(start_sites_dict)):
            all_start_end_list.append('s'+str(i))
        for i in range(len(end_sites_dict)):
            all_start_end_list.append('e'+str(i))

        length_list = get_length_list(read_length)

        no_subpaths = len(binary_subpath_list)
	#print '>>>>>>>>>>>start_end_nodes_subpath_list:',start_end_nodes_subpath_list
        Fs_list = [get_index_in_loc_list(x[0]) for x in start_end_nodes_subpath_list if len(x) > 0 ]

        #print 'Fs_list=',Fs_list

        p_arr = [0.2] * len(start_sites_dict)
        
        ###################### End site computation from junction reads
        q_arr = [0.2] * (len(end_sites_dict)-1)
        
        if disable_end_sites_computation:
            end_sites_arr = get_q_vector_from_junction_reads(strand_dir,subexon_ids_dict,loc_list,end_sites_dict,subpath_list,subpath_freq_list)
            q_arr = [row[2] for row in end_sites_arr]
            #print '(before) q_arr =',q_arr
        
        
        ###################### 
        
        prev_likelihood = 0.
        curr_likelihood = 0.
        num_iter = 1000
        step = 0
        count = 0
        max_count = 30
        epsilon = 1e-5
        offset_val = 0.0001 # to avoid division by zero error
        
        ###################### Step 5: Run EM
        while step < num_iter:
            
            
            prev_likelihood = curr_likelihood
            curr_likelihood = 0.

            PS = [0. for _ in range(no_subpaths)]
            A = [0. for _ in range(len(start_sites_dict))]
            B = [0. for _ in range(len(start_sites_dict))]
            C = [0. for _ in range(len(end_sites_dict)-1)]
            D = [0. for _ in range(len(end_sites_dict)-1)]
            E = 0.
            F = 0.
            numerator_e = 0.
            numerator_f = 0.


            f00,f01,f10,f11 = compute_f00_f01_f10_f11(p_arr, q_arr)
            PZ = populate_PZ(q_arr,p_arr)
            w = populate_ws(q_arr,p_arr,binary_subpath_list, start_end_nodes_subpath_list,subpath_list)
            ls_in,ls_out = get_expected_suffix_length(length_list,PZ,f00,f01,f10,f11)
            lp_in,lp_out = get_expected_prefix_length(length_list,PZ,f00,f01,f10,f11)

            for i in range(no_subpaths):
                start_end_list = start_end_nodes_subpath_list[i]
                freq = subpath_freq_list[i]
                PS[i] = (1-p_arr[0])*f01[0][Fs_list[i]]*w[i] + (p_arr[0])*f11[0][Fs_list[i]]*w[i]


                if Fs_list[i]==0:
                    numerator_e+=freq*1.
                else:
                    numerator_e+=freq*(p_arr[0]*f11[0][Fs_list[i]]*w[i])/PS[i]
                    numerator_f+=freq*((1-p_arr[0])*f01[0][Fs_list[i]]*w[i])/PS[i]
                for start_end_value in all_start_end_list:#all_start_end_list:
                    if start_end_value == 's0' or start_end_value == ('e'+str(len(end_sites_dict)-1)):
                        continue
                    elif start_end_value[0] == 's':
                        start_value = start_end_value
                        m = int(start_value[1:])
                        sm = start_sites_list[m] # location val
                        Ism = get_index_in_loc_list(start_value)-1
                        A[m] += freq * builder_start_01(binary_subpath_list[i], start_end_list, Fs_list[i], sm, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]                
                        B[m] += freq * builder_start_00(binary_subpath_list[i], start_end_list, Fs_list[i], sm, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]
                    elif start_end_value[0] == 'e' and not disable_end_sites_computation:
                        end_value = start_end_value
                        m = int(end_value[1:])
                        em = end_sites_list[m] # location val
                        Iem = get_index_in_loc_list(end_value)-1

                        C[m] += freq * builder_end_10(binary_subpath_list[i], start_end_list, Fs_list[i], em, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]
                        D[m] += freq * builder_end_11(binary_subpath_list[i], start_end_list, Fs_list[i], em, w[i],p_arr, q_arr, f00, f01, f10, f11)/PS[i]

                    #else:
                        #print "Unknown error"

            
            
            ## BIG NOTE: I ADDED OFFSET HERE, CHECK IF THAT'S RIGHT!!!
            E = numerator_e/(ls_in[0]+offset_val)
            F = numerator_f/(ls_out[0]+offset_val)
            p_arr[0] = (E+offset_val)/(E+F+2*offset_val)

            for m in range(1,len(start_sites_dict)):
                Ism = get_index_in_loc_list('s'+str(m))-1
                
                if Ism+1 >= len(ls_in) or (lp_out[Ism] + ls_in[Ism+1]) == 0:
                    #print "Division by zero", 'Ism=',Ism,'m=',m
                    A[m] = 0.
                else:
                    A[m] /= (lp_out[Ism] + ls_in[Ism+1])
            
                if Ism+1 >= len(ls_in) or (lp_out[Ism] + ls_out[Ism+1]) == 0:
                    #print "Division by zero", 'Ism=',Ism,'m=',m
                    B[m] = 0.
                else:
                    B[m] /= (lp_out[Ism] + ls_out[Ism+1])
                
                p_arr[m] = (A[m]+offset_val)/(A[m]+B[m]+2*offset_val)


            if not disable_end_sites_computation:
                for m in range(len(end_sites_dict)-1):
                    Iem = get_index_in_loc_list('e'+str(m))-1
                    C[m] /= (lp_in[Iem] + ls_out[Iem+1])
                    D[m] /= (lp_in[Iem] + ls_in[Iem+1])
                    q_arr[m] = (C[m]+offset_val)/(C[m]+D[m]+2*offset_val)

            for i in range(no_subpaths):
                if subpath_freq_list[i] * PS[i] > 0 :
			curr_likelihood += math.log(subpath_freq_list[i] * PS[i]) 


            if curr_likelihood-prev_likelihood <= epsilon:
                count += 1
            else:
                count = 0

            step += 1

            if count == max_count:
                #print "Breaking after %d iterations"%step
                break


        #print '\n\nStart Site probabilities = '
        #print str(['{0:0.2f}'.format(i) for i in p_arr]).replace("'", "")

        #print '\nEnd Site probabilities = '
        #print str(['{0:0.2f}'.format(i) for i in q_arr]).replace("'", "")

        #print '\n'
       #print '(after) q_arr =',q_arr
        p_arr_list.append(p_arr)
        q_arr_list.append(q_arr)        
        subpath_freq_list = bootstrap(orig_subpath_freq_list)

    return start_sites_dict, end_sites_dict, p_arr_list, q_arr_list


# In[7]:


def write_output(gene_id, start_sites_dict, end_sites_dict, p_arr_list, q_arr_list, out_DIR, no_steps,out_file,add_start_stop_nodes,chr_id,strand):
    sorted_ss_dict = sorted(start_sites_dict.items(), key=operator.itemgetter(1))
    sorted_es_dict = sorted(end_sites_dict.items(), key=operator.itemgetter(1))

    if  out_file == '':
        out_filename = out_DIR + gene_id+'.csv'
    else:
        out_filename = out_DIR + gene_id+'_'+out_file+'.csv'
        

    f = open(out_filename, 'w')
    writer = csv.writer(f)
    writer.writerow(['step', 'index', 'strand' ,'chr','splice site locus', 'usage estimate'])
    
    ss_count = len(sorted_ss_dict)
    es_count = len(sorted_es_dict)-1
    
    start_idx = 0
    if add_start_stop_nodes:
        ss_count -= 1 # to skip the start site of the last fake node
        es_count -= 1 # to skip the end site of the last fake node
        start_idx = 1 # to skip the first start and end sites of the first fake node
        

    for curr_step in xrange(no_steps+1):
        for i in range(start_idx,ss_count):
            writer.writerow([curr_step, 's'+str(sorted_ss_dict[i][1]-start_idx),strand,chr_id, sorted_ss_dict[i][0], p_arr_list[curr_step][i]])

        for i in range(start_idx,es_count):
            writer.writerow([curr_step, 'e'+str(sorted_es_dict[i][1]-start_idx),strand,chr_id ,sorted_es_dict[i][0], q_arr_list[curr_step][i]])


    f.close()
    return out_filename


# In[8]:


def filter_subpath_list_by_subexon_ids(subpath_list, subpath_freq_list, subexon_ids):
    indices = []
    for i in range(len(subpath_list)):
        if set(subpath_list[i]) <= set(subexon_ids):
            indices.append(i)
    subset_subpath_list = [subpath_list[idx] for idx in indices]
    subset_subpath_freq_list = [subpath_freq_list[idx] for idx in indices]
    return subset_subpath_list,subset_subpath_freq_list


# In[9]:


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', type=str, help='Input refined gtf file.', required=True)
    parser.add_argument('--count_file', type=str, help='Input read count file (generated from bam file).', required=True)
    parser.add_argument('--out_dir', type=str, required=True,help='Output directory prefix.')
    parser.add_argument('--read_len', type=str, help='Input read length.', required=True)
    parser.add_argument('--gene_id', type=str, help='Input gene ID, use this parameter when running McSplicer on a single gene. Default, run on all genes provided in the gtf annotation file.', default = '' ,required=False)
    #parser.add_argument('--gene_list', type=str, help='Input file with gene IDs, where each gene ID is written in a separate line, e.g., gene1\\ngene2\\ngene3, use this parameter when running McSplicer on multiple genes.', default = '',required=False)
    parser.add_argument('--bootstraps', type=int, default="1",help='Number of bootstraps')
    #parser.add_argument('--add_start_end_fake_nodes', type=int,default=0,help='Add a reference points at the start and end of each gene.')
    #parser.add_argument('--use_junction_reads_for_end_site_estimates',type=int,default=0,help='Compute end site probabilities via junction reads data.')
    parser.add_argument('--prefix', type=str, default="", help='Output file prefix.')
    #args = parser.parse_args()
    args, unknown = parser.parse_known_args()

    return args


# In[10]:



def read_gene_list(gene_ids_file):

    gene_list = []

    with open(gene_ids_file,'rb') as in_file:
        for line in in_file:
            if line != '':
#                line = line[:-2] #remove \n
                gene_id = line.split()[0]
                gene_list.append(gene_id)

    return gene_list






if __name__ == "__main__":

    params = vars(get_args())
    gtf_file = params['gtf']
    cnt_file = params['count_file']
    #input_gene_ids_file = params['gene_list']
    gene_id = params['gene_id']
    out_dir = params['out_dir']+'/'
    read_length = int(params['read_len'])
    no_steps = params['bootstraps']
    
    #add_start_stop_nodes = bool(int(params['add_start_end_fake_nodes']))
    #disable_end_sites_computation= bool(int(params['use_junction_reads_for_end_site_estimates']))
    out_file = params['prefix']
    disable_end_sites_computation = True
    add_start_stop_nodes = True

   

  

    gene_list = []

    #if input_gene_ids_file != '':
	#gene_list = read_gene_list(input_gene_ids_file)

    all_gene_dict = get_all_genes_dict(gtf_file) # all gene data dictionary
    
    if gene_id != '':
	gene_list.append(gene_id)
    else:
	gene_list = [gene_id for gene_id in all_gene_dict]

    
    
    for gene_id in gene_list:
        
        
        if not os.path.exists(out_dir):
            #continue
            os.system('mkdir -p %s'%out_dir)
        
        if gene_id not in all_gene_dict:
            print 'Gene %s not found in the refined gtf!'%gene_id
            continue


        print 'Running McSplicer on gene id',gene_id,'...'


        gene_datalist = all_gene_dict[gene_id]


    
        if not gene_datalist:
            print 'No entry found for the gene %s in the provided refined gtf file %s'%(gene_id,gtf_file)
            continue

        #try:
        start_sites_dict, end_sites_dict, p_arr_list, q_arr_list = run_EM_bootstrap(cnt_file, gene_id, no_steps, read_length,gene_datalist,disable_end_sites_computation,add_start_stop_nodes)
        #except IndexError:
         #    print '>>>>>>>>Index out of range error, skipped gene %s'%gene_id             

        if len(start_sites_dict) > 0:
	    chr_id = gene_datalist[0][-1]
            strand = gene_datalist[0][1]
            out_file = write_output(gene_id, start_sites_dict, end_sites_dict, p_arr_list, q_arr_list, out_dir, no_steps,out_file,add_start_stop_nodes,chr_id,strand)

            print '\tOutput is written to the folder',out_file
        else:
            print '\tNo counts are found for gene', gene_id

        #break
        

