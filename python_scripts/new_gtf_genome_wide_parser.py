
# coding: utf-8

# In[1]:

#!/usr/bin/env python

"""new_gtf_genome_wide_parser.py: New gtf parses much faster and more effecient."""

__authos__      = "Israa Alqassem"
__copyright__   = "Copyright 2017, McSplicer"


import csv
import numpy as np
import time



def get_all_genes_dict(gtf_file):

    gene_dict = {}

    with open(gtf_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for line in reader:
            if line[2] == 'subexon':
                start_site = int(line[3])
                end_site = int(line[4])
                strand_dir = line[6]

                gene_id = ''
                chr_id = line[0]

                feature_list = line[8].split(';')
                for feature in feature_list:
                    tag_val = feature.split()
                    if len(tag_val) == 2:

                        tag = tag_val[0]
                        value = tag_val[1]

                        if tag=='SpliceEnd':
                            splice_end = value[1:2] # R, L, or B
                        elif tag=='NodeId':
                            subexon_id = int(value)
                        elif tag=='transcript_id':
                            trans_id = value[1:-1]
                        elif tag == 'gene_id':
                            gene_id =  value[1:-1] #remove double qouts

                if gene_id != '':
                    if gene_id not in gene_dict:
                        gene_dict[gene_id] = []


                    gene_dict[gene_id].append([subexon_id,strand_dir,splice_end,start_site,end_site,trans_id,chr_id])
                    #print 'gene_id',gene_id,'subexon_id',subexon_id,'strand_dir',strand_dir,'ss',start_site,'es',end_site,'splice_end',splice_end,'trans_id',trans_id

    return gene_dict




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
    #print 'index','location'
    for location in location_list:
        #if location not in loc_index_dict:
            #print index, location
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


def get_gene_data(gene_id,gene_datalist):

    #for gene_id in gene_dict.keys():

    subexon_ids_dict = {}
    start_sites = []
    end_sites = []

    #print '>>>>>',gene_id
    for row in gene_datalist:
        subexon_id = row[0]
        strand_dir = row[1]
        splice_end = row[2]
        start_site = row[3]
        end_site = row[4]
        trans_id = row[5]


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

        if strand_dir == '+': # Forward strand
            subexon_ids_dict[subexon_id] = [start_site,end_site]

            if splice_end=='R':
                end_sites.append(end_site)
            elif splice_end=='L':
                start_sites.append(start_site)

            elif splice_end=='B':
                start_sites.append(start_site)
                end_sites.append(end_site)

            elif splice_end !='-': # dash means internal exon, just ignore it, otherwise show error
                print('Error: Splice end value must be L, R, or B. Undefined splice end -> ' + splice_end)

        elif strand_dir == '-': # Reverse strand
            subexon_ids_dict[subexon_id] = [end_site,start_site]

            if splice_end=='R':
                end_sites.append(start_site)
            elif splice_end=='L':
                start_sites.append(end_site)

            elif splice_end=='B':
                start_sites.append(end_site)
                end_sites.append(start_site)

            elif splice_end !='-': # dash means internal exon, just ignore it, otherwise show error
                print('Error: Splice end value must be L, R, or B. Undefined splice end -> ' + splice_end)


    start_sites = list(np.unique(start_sites))
    end_sites = list(np.unique(end_sites))


    if strand_dir == '-':
        end_sites.sort(reverse=True)
        start_sites.sort(reverse=True)
    elif strand_dir == '+':
        end_sites.sort()
        start_sites.sort()


    loc_index_dict, start_sites_dict, end_sites_dict = create_location_dicts(start_sites,end_sites, strand_dir)



    return strand_dir,loc_index_dict, start_sites_dict, end_sites_dict, subexon_ids_dict

#start_time = time.time()
#print "--- %s seconds ---" % str('{0:0.2f}'.format(time.time() - start_time))
