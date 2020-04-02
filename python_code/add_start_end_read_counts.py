
# coding: utf-8

def preprocess_read_counts(subexon_ids_dict,subexon_trans_ids_dict,read_length,subpath_list,subpath_freq_list):
    '''
      This function add fake start and end node to solve the  identifiability issue
      subexon_trans_ids_dict: key: transcript ID, value: a list of subexon_ids
      ASSUMPTION: I assume the subexons are ordered by their coordinates, e.g., subexon 1 appears before subexon 2,
                  it's better to add a few lines to order the subexons based on their coordinates not ID.

     Returns new subpath_list and subpath_freq_list, where new read counts are added from
     the fake start node to the starting subexons in each transcript, and
     from the ending subexons in each transcript to the fake end node
     
     Read counts are added based on the number of reads per base pair and based on read length, i.e.,
        
        (Reads_count/[subexon_length-read_length+1]) *( read_length-1)
    
    '''
    
    #subexon_ids_dict,subexon_trans_ids_dict,read_length,subpath_list,subpath_freq_list
    #print 'subexon_ids_dict=',subexon_ids_dict
    #print 'subexon_trans_ids_dict=',subexon_trans_ids_dict
    #print 'subpath_list=',subpath_list
    #print 'subpath_freq_list=',subpath_freq_list
    
    
    
    internal_nodes = [sublist[1:-1] for sublist in subexon_trans_ids_dict.values() if len(subexon_trans_ids_dict.values()) > 2]
    internal_nodes = [i for sublist in internal_nodes for i in sublist]

    starting_subexons = [min(val) for key,val in subexon_trans_ids_dict.iteritems()]
    fake_first_node = min(starting_subexons)-1

    #print 'starting_subexons>>>',starting_subexons
    #print 'fake_first_node>>>',fake_first_node

    ending_subexons = [max(val) for key,val in subexon_trans_ids_dict.iteritems()]
    fake_last_node = max(ending_subexons)+1

    #print 'ending_subexons>>>',ending_subexons
    #print 'fake_last_node>>>>',fake_last_node


    subpath_freq_dict = {}

    for i in range(0,len(subpath_list)):

        subpath,freq = subpath_list[i],subpath_freq_list[i]
        key = '-'.join(str(node) for node in subpath)
        subpath_freq_dict[key] = freq

        first_subexon = subpath[0]
        f_start,f_end = subexon_ids_dict[first_subexon]
        f_subexon_len = abs(f_start-f_end)


        last_subexon = subpath[-1]   #last subexon in a path
        l_start,l_end = subexon_ids_dict[last_subexon]
        l_subexon_len = abs(l_start-l_end)


        #################################################################################################
        ######################################### Adding start node #####################################

        if first_subexon not in internal_nodes and first_subexon in starting_subexons:
            # Case a.1: 
            ## If the first subexon's length > f_subexon_len-read_length+1
            ## Add new read counts from the fake start node to the first subexon
            if len(subpath) == 1 and (f_subexon_len-read_length+1) > 0:

                new_freq = round(freq/(f_subexon_len-(read_length-1))*(read_length-1))
                subpath.append(fake_first_node)
                new_key = '-'.join(str(node) for node in sorted(subpath))
                subpath_freq_dict[new_key] = new_freq


            # Case a. 2: 
            ## If the first subexon's length < f_subexon_len-read_length+1
            ## recursively consider the following subexon(s)
            elif len(subpath) > 1 and (f_subexon_len-read_length+1) <= 0:

                total_len = -read_length+1
                node_list = []

                for j in range(0,len(subpath)):
                    node = subpath[j]
                    b1,b2 = subexon_ids_dict[node]
                    total_len += abs(b1-b2)
                    node_list.append(node)

                    if total_len > 0:
                        break

                if  total_len > 0 and len(node_list) > 0:
                    new_freq = round((freq/total_len)*(read_length-1))
                    node_list.append(fake_first_node)
                    new_key = '-'.join(str(node) for node in sorted(node_list))
                    subpath_freq_dict[new_key] = new_freq



        #################################################################################################
        ######################################### Adding end node #####################################

        if last_subexon not in internal_nodes and last_subexon in ending_subexons:
            # Case b.1: 
            ## If the last subexon's length > l_subexon_len-read_length+1
            ## Add new read counts from the last subexon to the end fake node
            if  len(subpath) == 1 and (l_subexon_len-read_length+1) > 0:
                new_freq = round(freq/(l_subexon_len-(read_length-1))*(read_length-1))
                subpath.append(fake_last_node)
                new_key = '-'.join(str(node) for node in subpath)
                subpath_freq_dict[new_key] = new_freq

            # Case b.2: 
            ## If the last subexon's length < l_subexon_len-read_length+1
            ## recursively consider the preceding subexon(s)
            elif len(subpath) > 1 and (l_subexon_len-read_length+1) <= 0:

                total_len = -read_length+1
                node_list = []
                for j in reversed(range(len(subpath))):
                    node = subpath[j]
                    node_list.append(node)
                    b1,b2 = subexon_ids_dict[node]
                    total_len += abs(b1-b2)

                    if total_len > 0:
                        break

                if  total_len > 0 and len(node_list) > 0:
                    new_freq = round((freq/total_len)*(read_length-1))
                    node_list.append(fake_last_node)
                    new_key = '-'.join(str(node) for node in sorted(node_list))
                    subpath_freq_dict[new_key] = new_freq


    new_subpath_list = []
    new_subpath_freq_list = []
    for key,val in subpath_freq_dict.iteritems():
        subpath = [int(node) for node in key.split('-')]
        new_subpath_list.append(subpath)
        new_subpath_freq_list.append(val)
    
    
    return new_subpath_list,new_subpath_freq_list





def add_starting_and_stopping_nodes(strand_dir,subexon_ids_dict,subexon_trans_ids_dict,start_sites_list,end_sites_list,loc_list,start_sites_dict,end_sites_dict,loc_index_dict):
    
    
    starting_flag = True
    ending_flag = True
    
    reverse_strand = False
    if strand_dir == '-' :
        reverse_strand = True
    
    ## Appending last node
    fake_last_node = ''
    if ending_flag:
        ending_subexons = [max(val) for key,val in subexon_trans_ids_dict.iteritems()]
        fake_last_node = max(ending_subexons)+1
        if strand_dir == '+' :
            start = end_sites_list[-1] + 200
            end = start + 5 #(read_length-1)
        else:
            start = end_sites_list[-1] - 200
            end = start - 5 #(read_length-1)
            
        start_sites_dict[start] = max(start_sites_dict.values())+1
        end_sites_dict[end] = max(end_sites_dict.values())+1
        loc_index_dict[start] = max(loc_index_dict.values())+1
        loc_index_dict[end] = max(loc_index_dict.values())+1
        subexon_ids_dict[fake_last_node]=[start,end]
        start_sites_list.append(start)
        end_sites_list.append(end)
        loc_list.append(start)
        loc_list.append(end)
        
        
    ## Appending first node
    fake_first_node = ''
    if starting_flag:
        starting_subexons = [min(val) for key,val in subexon_trans_ids_dict.iteritems()]
        fake_first_node = min(starting_subexons)-1
        if strand_dir == '+' :
            start = start_sites_list[0] - 200
            end = start + 5
        else:
            start = start_sites_list[0] + 200
            end = start - 5
        
        # Shift the indices by 1 to add new indices for the fake node
        for key,val in start_sites_dict.iteritems():
            start_sites_dict[key] = val + 1
        for key,val in end_sites_dict.iteritems():
            end_sites_dict[key] = val + 1
        for key,val in loc_index_dict.iteritems():
            loc_index_dict[key] = val + 2
        
        
        start_sites_dict[start] = min(start_sites_dict.values())-1
        end_sites_dict[end] = min(end_sites_dict.values())-1
        loc_index_dict[end] = min(loc_index_dict.values())-1
        loc_index_dict[start] = min(loc_index_dict.values())-1
        
        subexon_ids_dict[fake_first_node]=[start,end]
        start_sites_list.append(start)
        end_sites_list.append(end)
        loc_list.append(start)
        loc_list.append(end)
        
        # Sort again
        end_sites_list = sorted(end_sites_dict.keys(), reverse=reverse_strand)
        start_sites_list = sorted(start_sites_dict.keys(), reverse=reverse_strand)
        loc_list = sorted(loc_index_dict.keys(), reverse=reverse_strand)

    return subexon_ids_dict,start_sites_dict,end_sites_dict,loc_index_dict,start_sites_list,end_sites_list,loc_list

