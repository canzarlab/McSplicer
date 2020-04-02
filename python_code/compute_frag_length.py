def get_min_max_frag_length(subpath_l,subpath_r,left_frag_len_list,right_frag_len_list,read_length=100):
    
    '''
        To find the min fragment length for the left read we try to push it the right most position
        while the read still spans over the same subexons. The opposite is correct for the right read.
    
    
        P.S. This code contains lots of redunduncy and need to be cleaned
    
    '''
    
    
    # If left and right reads don't overlap
    if len(list(set(subpath_l).intersection(subpath_r))) == 0:
        #min_frag_length = 1+sum(left_frag_len_list[1:])+sum(right_frag_len_list[0:-1])+1
        #max_frag_length = (read_length - 1)+sum(left_frag_len_list[1:])+sum(right_frag_len_list[0:-1])+(read_length - 1)
        
        ## I should add a variable tmp_total same way we do the other part
        
        
        min_frag_length  = 0

        if len(subpath_l) > 1:
            if sum(left_frag_len_list[1:]) < read_length-1:
                min_frag_length += read_length
            else:
                min_frag_length += (1+sum(left_frag_len_list[1:]))
        else:
            min_frag_length += read_length
            
            
        if len(subpath_r) > 1:
            if sum(right_frag_len_list[0:-1]) < read_length-1: # sum(list[0:-1]) will exclude the last element
                min_frag_length += read_length
            else:
                min_frag_length += (sum(right_frag_len_list[0:-1])+1)
        else:
            min_frag_length += read_length
            
           
        max_frag_length = 0
        
        in_between_l = sum(left_frag_len_list[1:-1])
        
        if len(subpath_l) > 1:
            if left_frag_len_list[0] < read_length - in_between_l - 1:
                max_frag_length += left_frag_len_list[0]
            else: 
                max_frag_length += (read_length - in_between_l - 1)
        else:
            max_frag_length += left_frag_len_list[0]

            
        in_between_r = sum(right_frag_len_list[1:-1])
        
        if len(subpath_r) > 1:
            if right_frag_len_list[-1] < read_length - in_between_r - 1:
                max_frag_length += right_frag_len_list[-1]
            else: 
                max_frag_length += (read_length - in_between_r - 1 )
        else:
            max_frag_length += right_frag_len_list[0]

        max_frag_length += (sum(left_frag_len_list[1:])+sum(right_frag_len_list[0:-1]))
            
            
            
        
        return min_frag_length, max_frag_length
    
    
    min_dict_left = {}
    max_dict_left = {} 
    
    if len(subpath_l) == 1: 
        min_dict_left[subpath_l[0]] = read_length
        max_dict_left[subpath_l[0]] = min(read_length,left_frag_len_list[0])### TO DO; double check this
        
    else:
        
        tmp_total = sum(left_frag_len_list[1:-1])
        for _idx in range(len(subpath_l)):
            subexon = subpath_l[_idx]
            if _idx == 0: #first element in a subpath
                min_dict_left[subexon] = 1
                tmp_total += min_dict_left[subexon]
            elif _idx == len(subpath_l)-1: #last element
                min_dict_left[subexon] = read_length - tmp_total
                tmp_total += min_dict_left[subexon]
            else: # middle exon
                min_dict_left[subexon] = left_frag_len_list[_idx]

        
        tmp_total = sum(left_frag_len_list[1:-1])
        for _idx in range(len(subpath_l)):
            subexon = subpath_l[_idx]
            if _idx == 0: # first subexon
                max_dict_left[subexon] = min(left_frag_len_list[_idx],read_length - tmp_total - 1)
                tmp_total += max_dict_left[subexon]
            elif  _idx == len(subpath_l)-1: # last subexon
                max_dict_left[subexon] = max(1,read_length - tmp_total)
                tmp_total += max_dict_left[subexon]
            else:
                max_dict_left[subexon] = left_frag_len_list[_idx]
        
    
    max_dict_right = {}
    min_dict_right = {}
    
    if len(subpath_r) == 1:
        min_dict_right[subpath_r[0]] = read_length
        max_dict_right[subpath_r[0]] = max(read_length,right_frag_len_list[0])
    else:        
        
        tmp_total = sum(right_frag_len_list[1:-1])
        for _idx in range(len(subpath_r)):
            subexon = subpath_r[_idx]
            if _idx == 0: #first element in a subpath
                max_dict_right[subexon] = 1
                tmp_total += max_dict_right[subexon]
            elif _idx == len(subpath_r)-1: #last element
                max_dict_right[subexon] = read_length - tmp_total
                tmp_total += max_dict_right[subexon]
            else: # middle exon
                max_dict_right[subexon] = right_frag_len_list[_idx]

                
        tmp_total = sum(right_frag_len_list[1:-1])
        for _idx in range(len(subpath_r)):
            subexon = subpath_r[_idx]
            if _idx == 0: #first element in a subpath
                if right_frag_len_list[_idx] <= read_length-1:
                    min_dict_right[subexon] = right_frag_len_list[_idx]
                else:
                    min_dict_right[subexon] = (read_length-1) - tmp_total

                tmp_total += min_dict_right[subexon]

            elif _idx == len(subpath_r)-1: #last element
                min_dict_right[subexon] = read_length - tmp_total
                tmp_total += min_dict_right[subexon]
            else: # middle exon
                min_dict_right[subexon] = right_frag_len_list[_idx]



            
    print '>>left min, max',min_dict_left,max_dict_left
    print '>>right min,max',min_dict_right,max_dict_right
        
        
    keys = min_dict_right.keys()+ min_dict_left.keys()
    keys = unique(keys)

    min_frag_length = 0
    max_frag_length = 0
    for key in keys:
        if key in min_dict_right and key in min_dict_left:
            
            if key == subpath_l[-1] and key == subpath_r[0]:
                
                #if (len(subpath_r) != 1:
                 #   min_frag_length += max(min_dict_right[key],min_dict_left[key],left_frag_len_list[-1])    
                #else:
                   # min_frag_length += max(min_dict_right[key],min_dict_left[key])    
                    
                if len(subpath_r) != 1:
                    min_frag_length += max(min_dict_right[key],min_dict_left[key],left_frag_len_list[-1]) 
                else:
                    min_frag_length += max(min_dict_right[key],min_dict_left[key])
                    
                max_frag_length += max(max_dict_right[key],max_dict_left[key],left_frag_len_list[-1])
            else:
                min_frag_length += max(min_dict_right[key],min_dict_left[key])    
                max_frag_length += max(max_dict_right[key],max_dict_left[key])
                
        elif key in min_dict_right:
            min_frag_length += min_dict_right[key]
            max_frag_length += max_dict_right[key]
        else:
            min_frag_length += min_dict_left[key]
            max_frag_length += max_dict_left[key]


    if subpath_l == subpath_r:
        min_frag_length = read_length # Left and right read overlap
        
    print min_frag_length, max_frag_length

    return min_frag_length, max_frag_length