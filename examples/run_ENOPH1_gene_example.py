
import sys, os


 
## =======================
##  McSplicer example
## =======================



print 'Running McSplicer on ENOPH1 gene obtained from two Autism patients:\n\t1. A patiend with mutation (sample ID: ERR2902115_S36, patient ID: 14668)\
\n\t2. A patiend without the mutation (Sample ID: ERR2902117_S5, patient ID: 11414)'

## Example #1. ENOPH1 gene, exon skipping in the mutated sample, see the estimate of splice site at chromosome 4, splice site locus: 83378068
os.system("python2 ../python_code/McSplicer.py \
		--gtf ./refined_gtf/ERR2902115_S36_ENOPH1_refined.gtf \
		--count_file ./signature_counts/ERR2902115_S36.cnt \
		--gene_id ENSG00000145293\
		--out_dir ./out\
		--bootstraps 10 --read_len 151\
		--prefix mutated")

## Example #2. ENOPH1 gene, exon not skipped in the control sample, see the estimate of splice site at chromosome 4, splice site locus: 83378068
os.system("python2 ../python_code/McSplicer.py \
		--gtf ./refined_gtf/ERR2902117_S5_ENOPH1_refined.gtf \
		--count_file ./signature_counts/ERR2902117_S5.cnt \
		--gene_id ENSG00000145293 \
		--out_dir ./out \
		--bootstraps 10 \
		--read_len 151 \
		--prefix control")
