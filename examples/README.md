McSplicer Example
=================


Here we run McSplicer on an example of ENOPH1 gene obtained from two Autism patients:

1. A patien with mutation (sample ID: ERR2902115_S36, patient ID: 14668)
2. A patien without the mutation (Sample ID: ERR2902117_S5, patient ID: 11414)


All autism patient samples are available [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7351/samples/?s_page=1&s_pagesize=25)

### Output: ###

The output csv file contains the following columns:

 * Bootstrap step
 * Splice site index
 * Gene strand	
 * Chromosome 	
 * Splice site genome position	
 * McSplicer splice site usage estimate
 
 If you choose to run McSplicer with ```--bootstraps n```, step 0 in the output file corresponds to the estimates based on input data, and the following n steps correspond to bootstrapping results.
 
 Splice site index column represents the index of 3' start or 5' end splice sites as they appear in a gene according to their chronological order, e.g., s<sub>0</sub>, s<sub>1</sub>,..., e<sub>0</sub>, e<sub>1</sub>,..


### References: ###

Jaganathan, Kishore, et al. "Predicting splicing from primary sequence with deep learning." Cell 176.3 (2019): 535-548.


&copy; 2020 McSplicer





