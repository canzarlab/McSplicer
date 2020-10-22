
# coding: utf-8

# In[1]:


from new_gtf_genome_wide_parser import *
from subpath_info import *


# In[2]:


def get_ground_truth_probability_vector(tx_dict,start_end_full_path,loc_list):
    '''
        tx_dict: 
            key -> tx_id
            val -> [abundance, subexon_list, binary_path, start_end_path]
        
        start_end_full_path: a list of start and end sites with their indecies, e.g., s0, e0, s1, e2, etc.
        loc_list: a sorted gene sites list which represents subexons' start and end sites.
        
    This function returns a list of gene location, whether its start/end site, and the probability of that site
        
    '''

    ground_truth_prob = []

    for site in start_end_full_path[:-1]: # excluding the last end site
        denum = 0.
        num = 0.
        

        for tx in tx_dict:
            item = tx_dict[tx]
            abundance = item[0]
            binary_path = item[2]
            start_end_path = item[3]
            segment = -1

            #print 'binary_path:',binary_path
            #print 'start_end_path:',start_end_path
            
            if site in start_end_path:
                segment = start_end_path.index(site)
                #print 'segment-->',start_end_path.index(site)
                
            

            if site[0] == 's0' and site in start_end_path:
                denum += abundance
                num += abundance

            elif site[0] == 's0' and site not in start_end_path:
                denum += abundance

            elif site[0] == 's' and site in start_end_path:

                if len(binary_path) > segment:

                    #print segment,':',binary_path[segment-1],binary_path[segment]
                    if site == start_end_path[0]: # in case of altss, the tx doesn't start at s0
                        denum += abundance
                        num += abundance
                    elif binary_path[segment-1] == 0 and binary_path[segment] == 1:
                        denum += abundance
                        num += abundance
                    elif binary_path[segment-1] == 0 and binary_path[segment] == 0:
                        denum += abundance
                        
                    ##########################################    
                    ## just added this for testing
                    #elif binary_path[segment-1] == 1 and binary_path[segment] == 1:
                    #    print '>>>>>>>>>>>>>here1'
                    #    denum += abundance
                    #    num += abundance
                    ##########################################
                    
            elif site[0] == 's' and site not in start_end_path:
                    num += abundance
                    denum += abundance

            elif site[0] == 'e' and site in start_end_path:
                if len(binary_path) > segment:
                    
                    #print segment,':',binary_path[segment-1],binary_path[segment]
                    if binary_path[segment-1] == 1 and binary_path[segment] == 0:
                        denum += abundance
                        num += abundance
                    elif binary_path[segment-1] == 1 and binary_path[segment] == 1:
                        denum += abundance
                    ########################################## 
                    ## just added this for testing    
                    #elif binary_path[segment-1] == 1 and binary_path[segment] == 1:
                    #    print '>>>>>>>>>>>>>here2'
                    #    num += abundance
                    #    denum += abundance  
                    ##########################################
                else: # If end site comes at the end of this transcript,assume the next segment is zero
                    
                    if binary_path[segment-1]==1:
                        denum += abundance
                        num += abundance
                    

            elif site[0] == 'e' and site not in start_end_path:
                denum += abundance
                num += abundance      

        #print site,'=',num/denum        
        ground_truth_prob.append(num/denum)       

    ground_truth_prob.append(1.) # for the last end site
    
    loc_probability = []
    for i in range(len(loc_list)):
        loc_probability.append([loc_list[i],start_end_full_path[i],ground_truth_prob[i]])
    
    return loc_probability


# In[3]:


in_DIR = '/home/israa/PyScript/McSplicer/data/'
gtf_list = ['altss_b.refined.gtf','alttes_d.refined.gtf','altacc_b.refined.gtf','exskip_b.refined.gtf','altdon_b.refined.gtf']

refined_gtf = in_DIR+gtf_list[1]

all_gene_dict = get_all_genes_dict(refined_gtf) # all gene data dictionary
gene_id = 'A'
gene_datalist = all_gene_dict[gene_id]
#gene_datalist

tx_node_dict = {} # map between tx and its subexon ids
for row in gene_datalist:#'ENST00000371588'
    tx_id = row[5]
    if tx_id not in tx_node_dict:
        tx_node_dict[tx_id] = []
        
    tx_node_dict[tx_id].append(row[0])
#tx_node_dict

# tx_list
#unique([row[5] for row in gene_datalist])


# In[5]:


strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict = get_gene_data(gene_id,gene_datalist)

if strand_dir == '+' :
    end_sites_list = sorted(end_sites_dict.keys())
    start_sites_list = sorted(start_sites_dict.keys())
    loc_list = sorted(loc_index_dict.keys())
else:
    end_sites_list = sorted(list(end_sites_dict.keys()), reverse=True)
    start_sites_list = sorted(list(start_sites_dict.keys()), reverse=True)
    loc_list = sorted(list(loc_index_dict.keys()), reverse=True)
    
set_global_vars_in_subpath(strand_dir,loc_index_dict,start_sites_dict,end_sites_dict,subexon_ids_dict,{},end_sites_list,start_sites_list,loc_list)

print('loc_index_dict:\t\t',loc_index_dict)
print('start_sites_dict:\t',start_sites_dict)
print('end_sites_dict:\t\t',end_sites_dict)
print('subexon_ids_dict:\t',subexon_ids_dict)
print('loc_list:\t',loc_list)


# In[6]:


binary_full_path, start_end_full_path = get_subpath_info(list(subexon_ids_dict.keys()))
print(binary_full_path, start_end_full_path)
print(tx_node_dict)


# In[13]:


tx1_id = list(tx_node_dict.keys())[0]
tx1 = tx_node_dict[tx1_id]
binary_path1, start_end_path1 = get_subpath_info(tx1)
#print binary_path1, start_end_path1

tx2_id = list(tx_node_dict.keys())[1]
tx2 = tx_node_dict[tx2_id]
binary_path2, start_end_path2 = get_subpath_info(tx2)
#print binary_path2, start_end_path2

tx_dict = {tx1_id:[0.95,tx_node_dict[tx1_id],binary_path1,start_end_path1],
    tx2_id:[0.05,tx_node_dict[tx2_id],binary_path2,start_end_path2]}


# In[14]:


get_ground_truth_probability_vector(tx_dict,start_end_full_path,loc_list)

#print loc_list
#print binary_full_path
#print start_end_full_path
#print ground_truth_prob
#for i in range(len(loc_list)):
    #print loc_list[i],start_end_full_path[i],ground_truth_prob[i]


# In[16]:


if 's0' in start_end_path1:
    print(start_end_path1.index('s0'))
else:
    print(-1)

