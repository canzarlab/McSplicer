import os,csv,math
import numpy as np
from os import listdir
from os.path import isfile, join
from collections import defaultdict



def read_bootstrap_output(filename):
    
    first_line = True
    loc_prob_dict = {}
    loc_idx_dict = {}
    p_hat_dict = {}

    with open(filename, "rb") as f:
        reader = csv.reader(f, delimiter=",")
        for row in reader:
            if first_line:
                first_line = False
                continue
            
            step = int(row[0])
            idx = row[1]
            loc = int(row[2])
            prob = float(row[3])
            
            key = str(loc) #str(loc)+'_'+idx

            if step == 0:
                p_hat_dict[key] = prob

            if key not in loc_prob_dict:
                loc_prob_dict[key] = []

            loc_prob_dict[key].append(prob)
            
    return loc_prob_dict,p_hat_dict




def read_bootstrap_output_modified(mcsplicer_out_dir,gene_id):
    '''
        Modified for reading the output for shortened gene list
    '''
    
    file_list = [f for f in listdir(mcsplicer_out_dir) if isfile(join(mcsplicer_out_dir, f))]
    gene_output_files = [item for item in file_list if item.startswith(gene_id)]
    gene_output_files

    loc_prob_dict = defaultdict(list)
    loc_idx_dict = {}
    p_hat_dict = {}

    for filename in gene_output_files:

        #print '>> file:',filename

        filename = mcsplicer_out_dir+filename

        first_line = True

        with open(filename, "rb") as f:
            reader = csv.reader(f, delimiter=",")
            #gene_id = filename.split('_')[0]

            for row in reader:
                if first_line:
                    first_line = False
                    continue

                step = int(row[0])
                idx = row[1]
                loc = int(row[2])
                prob = float(row[3])

                key = str(loc) #str(loc)+'_'+idx
                #key = gene_id+'_'+key 

                if step == 0:
                    p_hat_dict[key] = prob

    #             if key not in loc_prob_dict:
    #                 loc_prob_dict[key] = []

                loc_prob_dict[key].append(prob)



    return loc_prob_dict,p_hat_dict



# def get_McSplicer_mean_and_std_estimates_for_splice_site(mcsplicer_out_dir,gene_id, splice_site):
    
#     mcsplicer_filename = mcsplicer_out_dir+gene_id+'.csv'
    
#     if not os.path.exists(mcsplicer_filename):
#         return -1,-1
    
#     loc_prob_dict,p_hat_dict = read_bootstrap_output(mcsplicer_filename)
#     try:
#         #return np.median(loc_prob_dict[splice_site]),np.std(loc_prob_dict[splice_site])
#         return np.mean(loc_prob_dict[splice_site]),np.std(loc_prob_dict[splice_site])
#         #return loc_prob_dict[splice_site][0],np.std(loc_prob_dict[splice_site])
    
#     except KeyError as e:
#         return -1,-1



def get_McSplicer_mean_and_std_estimates_for_splice_site(mcsplicer_out_dir,gene_id, splice_site, modified_version=True):
    

    if modified_version:
        loc_prob_dict,p_hat_dict = read_bootstrap_output_modified(mcsplicer_out_dir,gene_id)
        
    else:
        
        mcsplicer_filename = mcsplicer_out_dir+gene_id+'.csv'
    
        if not os.path.exists(mcsplicer_filename):
            return -1,-1
    
        loc_prob_dict,p_hat_dict = read_bootstrap_output(mcsplicer_filename)
        
    try:
        #return np.median(loc_prob_dict[splice_site]),np.std(loc_prob_dict[splice_site])
        return np.mean(loc_prob_dict[splice_site]),np.std(loc_prob_dict[splice_site])
        #return loc_prob_dict[splice_site][0],np.std(loc_prob_dict[splice_site])
    
    except KeyError as e:
        return -1,-1






def get_ground_truth_probabilities(filename):
    
    if not os.path.exists(filename):
        return {}

    first_line = True

    loc_prob_dict = {}

    with open(filename, "rb") as f:
        reader = csv.reader(f, delimiter=",")
        for row in reader:
            if first_line:
                first_line = False
                continue

            loc = row[0]
            idx = row[1]
            p = row[2] # sample1, the row has 19 other values for the other samples, we don't look into them now 
            #print loc,idx,p

            #key = str(loc)+'_'+idx
            key = str(loc)
            loc_prob_dict[key]=p
    return loc_prob_dict
