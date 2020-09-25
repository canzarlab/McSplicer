

## generate_genomewide_reads_polyester.R (Derek, modified by Israa)

Simple script to simulates single- or paired-end reads (genomewide) using polyester.

**Input**: <br/> 1. Ground truth file: tab seperated file contains the trascript ids (column header should be transcript_id) and the corresponding abundance estimates including read counts. This scripts uses the read count (column index 5) and multiply it by scale parameter, however it can be modified to use TPM values instead. make sure to modify the code to match your file format. <br/>
2. cDNA fasta sequence file, e.g., Homo_sapiens.GRCh38.91.cdna.fa.<br/>
3. Scale to multiply the read count per transcript. <br/>
4. Output directory, make sure it has enough space. <br/>
5. TRUE or FALSE (capitalized) whether to simulate paired-end reads or single-end reads. <br/>
6. Read length. <br/>
7. Mean fragment length. <br/>
8. The standard deviation of the fragment length. <br/>
9. Seed to make the results reproducable. <br/>

**Output**: <br/>
Check <a target="_blank" href="https://github.com/alyssafrazee/polyester#output">Polyester Output</a>.<br/>

**Addtional notes**: <br/>
You have to install the following R packages: polyester, Biostrings.<br/>
Please check the documentation of simulate_experiment function from polyester.<br/> 
In this script no fold change is specified. <br/>
Tested on R version 3.5.1.
```
Rscript sim_genomewide_reads.R ./data/tpm.csv ./reference/Homo_sapiens.GRCh38.91/Homo_sapiens.GRCh38.91.cdna.fa 10 /out_dir/polyester_sim_reads/ FALSE 100 250 25 43
```





Authors:
----------------------------
* Derek Reiman (dreima2@uic.edu)
* Israa Alqassem (alqassem.isra@gmail.com)

