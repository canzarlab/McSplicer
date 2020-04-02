
# coding: utf-8


import csv
import numpy as np



def read_bootstrap_output(filename,location):
    
    first_line = True
    loc_prob_dict = {}
    loc_idx_dict = {}
    p_hat = -1

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

            if loc == location and step == 0:
                p_hat = prob

            if loc not in loc_prob_dict:
                loc_prob_dict[loc] = []

            loc_prob_dict[loc].append(prob)
            
    return loc_prob_dict,p_hat
    


# In[10]:

gtf_file_list = ['exskip.refined.gtf','alttss.refined.gtf','altacc.refined.gtf','altdon.refined.gtf','alttes.refined.gtf']


#GTF File:  alttss.refined.gtf
ground_truth_prob_list = [0.6016260162601627,0.526032315978456, 0.473967684021544,0.5348101265822784,0.4692028985507246]

loc_list = [1500,1000,2300,1100,2535]


#GTF File:  alttss.refined.gtf
#Ground truth prob = [0.526032315978456, 0.473967684021544]
#GTF File:  altacc.refined.gtf
#Ground truth prob = [0.526032315978456, 0.473967684021544]
#GTF File:  altdon.refined.gtf
#Ground truth prob = [0.5348101265822784, 0.4651898734177215]
#GTF File:  alttes.refined.gtf
#Ground truth prob = [0.4692028985507246, 0.5307971014492754]
#GTF File:  exskip.refined.gtf
#Ground truth prob = [0.6016260162601627, 0.3983739837398374]


total_num_files = 100
#bootstrap_steps = 1000
gene_id = 'A'
c = .95 # 95% confidence interval
Z = 1.96



for idx in range(len(gtf_file_list)):
    
    gtf_file = gtf_file_list[idx]
    print '******',gtf_file.split('.')[0].upper(),'******'
    
    #gtf_file = gtf_file_list[0]
    loc = loc_list[idx]
    p_true = ground_truth_prob_list[idx]
    
    count = 0.

    #print '[i]\tp_true\t\t p_hat\t p_hat-(2.*std)\t p_hat+(2.*std)\t  flag'
    print '[i]\tp_true\t\t p_hat\t p_hat-d2\t p_hat-d1\t  flag'
    
    p_hat_list = []

    for i in range(total_num_files):
        
        print '[%d]'%(i+1),

        filename = '../data/bootstrap_out/' + gtf_file.split('.')[0].upper()+'/'+gene_id+'_'+str(i)+'.csv'
        loc_prob_dict, p_hat = read_bootstrap_output(filename,loc)
        flag = 'False'
        n = len(loc_prob_dict[loc]) 
        p_hat_list.append(p_hat)

        ##### Method 1 #####
        #std = np.std(loc_prob_dict[loc])
        #if p_true >= (p_hat - Z*std) and p_true <= (p_hat + Z*std):
            #flag = 'True'
            #count += 1
        #print p_true,p_hat,p_hat-(Z*std),p_hat+(Z*std), flag
        
        
        ##### Method 2 #####
        ## distribution of P_hat - Ptrue can be estimated by Pb_i - Phat
        ## Method 3: change pbi-p_hat to pbi-p_true
        delta_list = [pbi-p_true for pbi in sorted(loc_prob_dict[loc])] #
        delta_list = sorted(delta_list)
        d1 = delta_list[int((1.-c)/2.*n)]
        d2 = delta_list[int((1.+c)/2.*n)]
        if p_true >= p_hat-d2 and p_true <= p_hat-d1:
            flag = 'True'
            count += 1   
        print p_true,p_hat,p_hat-d2,p_hat-d1, flag
    
    print '>>>>>>%.3f%%' % float(count/total_num_files*100.)
    print 'p_hat_avg=%.3f' % np.mean(p_hat_list)
    

 

