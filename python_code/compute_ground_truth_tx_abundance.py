
# coding: utf-8

# In[1]:


from new_gtf_genome_wide_parser import *
from subpath_info import *
import re
import os.path


# In[2]:


## Get a list of genes with differentially expressed transcripts
def gene_list(DE_genes_file):

    DE_gene_list = []

    with open(DE_genes_file, "rb") as in_file:
        for line in in_file:
            DE_gene_list.append(line[0:-1])
            
    return DE_gene_list


# In[3]:


def stringSplitByNumbers(x):
    r = re.compile('(\d+)')
    l = r.split(x)
    return [int(y) if y.isdigit() else y for y in l]


# In[4]:


def get_gene_tx_dict(abundance_file):
    
    gene_tx_dict = {}
    header_line = True
    
    with open(abundance_file,'rb') as in_file:
        for line in in_file:
            if header_line:
                header_line = False
                continue
            line = line.split() 

            tx_id = line[0]
            gene_id = line[1]
            tpm = float(line[5])

            if gene_id not in gene_tx_dict:
                gene_tx_dict[gene_id] = []
            gene_tx_dict[gene_id].append([tx_id,tpm])
            
    return gene_tx_dict


# In[5]:


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
                
                #if site == 's7':
                    #print 'segment=',segment
                    #print 'len(binary_path)=',len(binary_path) 
                    #print binary_path[segment-1] ,binary_path[segment] 

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

            elif site[0] == 's' and site not in start_end_path:
                    denum += abundance

            elif site[0] == 'e' and site in start_end_path:
                if len(binary_path) > segment:
                    
                    #print segment,':',binary_path[segment-1],binary_path[segment]
                    if binary_path[segment-1] == 1 and binary_path[segment] == 0:
                        denum += abundance
                        num += abundance
                    elif binary_path[segment-1] == 1 and binary_path[segment] == 1:
                        denum += abundance
                else: # If end site comes at the end of this transcript,assume the next segment is zero
                    
                    if binary_path[segment-1]==1:
                        denum += abundance
                        num += abundance
                    

            elif site[0] == 'e' and site not in start_end_path:
                denum += abundance
                num += abundance      

        
        if num==0. and denum==0.:
            #print 'site=',site,'num=',num,'denum=',denum
            ground_truth_prob.append(0.)
        else:
            ground_truth_prob.append(num/denum)       

    ground_truth_prob.append(1.) # for the last end site
    
    loc_probability = []
    for i in range(len(loc_list)):
        loc_probability.append([loc_list[i],start_end_full_path[i],ground_truth_prob[i]])
    
    return loc_probability


# In[6]:


def get_gene_transcripts_info(gene_id,abundance_file):
    
    gene_tx_dict = get_gene_tx_dict(abundance_file)
    
    if gene_id not in gene_tx_dict:
        return {},[],[],'not found'
    
    gene_tx_dict[gene_id].sort(key=lambda x: x[1],reverse=True)
    
    #print 'tx\tabundance\n',gene_tx_dict[gene_id]
    
    gene_datalist = all_gene_dict[gene_id]
    
    tx_node_dict = {} # map between tx and its subexon ids
    for row in gene_datalist:#'ENST00000371588'
        tx_id = row[5]
        if tx_id not in tx_node_dict:
            tx_node_dict[tx_id] = []

        tx_node_dict[tx_id].append(row[0])
    
    strand_dir1,loc_index_dict1, start_sites_dict1, end_sites_dict1, subexon_ids_dict1 = get_gene_data(gene_id,gene_datalist)
    
    if strand_dir1 == '+' :
        end_sites_list1 = sorted(end_sites_dict1.keys())
        start_sites_list1 = sorted(start_sites_dict1.keys())
        loc_list1 = sorted(loc_index_dict1.keys())
    else:
        end_sites_list1 = sorted(end_sites_dict1.keys(), reverse=True)
        start_sites_list1 = sorted(start_sites_dict1.keys(), reverse=True)
        loc_list1 = sorted(loc_index_dict1.keys(), reverse=True)

    set_global_vars_in_subpath(strand_dir1,loc_index_dict1,start_sites_dict1,end_sites_dict1,subexon_ids_dict1,{},end_sites_list1,start_sites_list1,loc_list1)
    
    gene_tx_dict[gene_id].sort(key=lambda x: x[1],reverse=True)
    
    tx_dict = {}

    for row in gene_tx_dict[gene_id]:
        tx_id = row[0]
        tx_abund = row[1]
        subexon_list = tx_node_dict[tx_id]
        binary_path, start_end_path = get_subpath_info(subexon_list)

        #if tx_id not in tx_dict:
        tx_dict[tx_id] = [tx_abund,subexon_list,binary_path,start_end_path]
        
    binary_full_path, start_end_full_path = get_subpath_info(subexon_ids_dict1.keys())

    
    
    return tx_dict,start_end_full_path,loc_list1,'found'


# In[7]:


def write_ground_truth_prob_tofile(out_DIR, gene_id,loc_prob_dict):
    
    out_filename = out_DIR + gene_id+'.csv'

    f = open(out_filename, 'w')
    writer = csv.writer(f)

    writer.writerow(['Index',"Location", "sample_1","sample_2","sample_3","sample_4","sample_5",
                    "sample_6","sample_7","sample_8","sample_9","sample_10","sample_11","sample_12",
                    "sample_13","sample_14","sample_15","sample_16","sample_17","sample_18","sample_19","sample_20"])


    for key in loc_prob_dict.keys():
        row = []
        loc = int(key.split('_')[0])
        idx = key.split('_')[1]
        row.append(loc)
        row.append(idx)
        for val in loc_prob_dict[key]:
            #print val[2],
            row.append(val[2])
        writer.writerow(row)


    f.close()


# In[ ]:

'''
#if __name__ == "__main__":

refined_gtf = '/data/israa/rsem_single_end_reads_sim/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_refined.gtf'
out_DIR = '/home/israa/Server/israa/genome_simulation/ground_truth_probabilities/'
DE_genes_file = '/home/israa/Server/sophia/simulated_data/DTU_single_cutoff/DTU_simulation_validation/genes_significant_isoform_switching.txt'
simulated_reads_DIR = '/home/israa/Server/sophia/simulated_data/DTU_single_cutoff/fullmodel/simulated_reads/'


all_gene_dict = get_all_genes_dict(refined_gtf) # all gene data dictionary

DE_gene_list = gene_list(DE_genes_file)

#len(DE_gene_list)    


# In[15]:


for gene_id in DE_gene_list:
    
    gene_id = 'ENSG00000072694'
    
    #file_path = out_DIR + gene_id+'.csv'
    #if os.path.exists(file_path):
     #   continue
        
    
    #gene_id = DE_gene_list[11]#'ENSG00000242372'
    print 'Computing ground truth probability for gene',gene_id,'...'

    all_samples_probabilities = []

    loc_prob_dict = {}

    for idx in range(1,21):

        abundance_file1 = simulated_reads_DIR+'/SRR6987574_sample'+str(idx)+'/SRR6987574_sample'+str(idx)+'.sim.isoforms.results'


        #gene_tx_dict = get_gene_tx_dict(abundance_file1)
        #gene_tx_dict[gene_id].sort(key=lambda x: x[1],reverse=True)
        #gene_tx_dict[gene_id]

        #gene_tx_dict = get_gene_tx_dict(abundance_file2)
        #gene_tx_dict[gene_id].sort(key=lambda x: x[1],reverse=True)
        #gene_tx_dict[gene_id]


        try: 
            tx_dict1,start_end_full_path1,loc_list1,val = get_gene_transcripts_info(gene_id,abundance_file1)
        
            if val == 'not found':
                print '>>>>>>here'
                break # break this loop and move to next gene
        except:
            break
            
        loc_probability1 = get_ground_truth_probability_vector(tx_dict1,start_end_full_path1,loc_list1)

        all_samples_probabilities.append(loc_probability1)

        sample = 'sample_'+str(idx)

        for row in loc_probability1:
            #print row
            key = str(row[0])+'_'+row[1]
            #print key
            if key not in loc_prob_dict:
                loc_prob_dict[key]= []

            loc_prob_dict[key].append([sample,row[1],row[2]])

   # write_ground_truth_prob_tofile(out_DIR, gene_id,loc_prob_dict)

    break

    #        site_prob_dict1 = {}
    #        for row in loc_probability1:
    #            site_prob_dict1[row[1]]=[row[0],row[2]]

    #        site_prob_dict2 = {}
    #        for row in loc_probability2:
    #            site_prob_dict2[row[1]]=[row[0],row[2]]

        #out_filename = out_DIR + gene_id+'.csv'


    #        for key in sorted(site_prob_dict1.keys(), key = stringSplitByNumbers):
    #            if key in site_prob_dict2:
    #                print site_prob_dict1[key][0],key,site_prob_dict1[key][1],site_prob_dict2[key][1]


# In[17]:


tx_dict1.keys()

'''