
import sys, os


 
## =======================
##  McSplicer example
## =======================

print('============================================================================================')
print('McSplicer example on ENOPH1 gene obtained from two Autism patients:\n\t1. A patiend with mutation (sample ID: ERR2902115_S36, patient ID: 14668)\
\n\t2. A patiend without the mutation (Sample ID: ERR2902117_S5, patient ID: 11414)')
print('============================================================================================')

print('\n*************First, run exonRefine to generate refined gtf*************\n')
os.system("../bin/exonRefine  ./gtf/ENOPH1_S5.gtf  -p ./gtf/ENOPH1_S5_refined")
os.system("../bin/exonRefine  ./gtf/ENOPH1_S36.gtf  -p ./gtf/ENOPH1_S36_refined")

print('\n*************Second, run sigCount to generate signature counts from read bam files*************\n')
os.system("../bin/sigcount ./bam/ENOPH1_S36.bam ./gtf/ENOPH1_S36_refined.gtf ./signature_counts/ENOPH1_S36")
os.system("../bin/sigcount ./bam/ENOPH1_S5.bam ./gtf/ENOPH1_S5_refined.gtf ./signature_counts/ENOPH1_S5")



print('\n*************Third, run McSplicer*************\n')

## Example #1. ENOPH1 gene, exon skipping in the mutated sample, see the estimate of splice site at chromosome 4, splice site locus: 83378068
os.system("python3 ../python_scripts/McSplicer.py \
		--gtf ./gtf/ENOPH1_S36_refined.gtf \
		--count_file ./signature_counts/ENOPH1_S36.cnt \
		--gene_id ENSG00000145293\
		--out_dir ./out\
		--read_len 151\
		--prefix mutated")

## Example #2. ENOPH1 gene, exon not skipped in the control sample, see the estimate of splice site at chromosome 4, splice site locus: 83378068
os.system("python3 ../python_scripts/McSplicer.py \
		--gtf ./gtf/ENOPH1_S5_refined.gtf \
		--count_file ./signature_counts/ENOPH1_S5.cnt \
		--gene_id ENSG00000145293 \
		--out_dir ./out \
		--read_len 151 \
		--prefix control")
