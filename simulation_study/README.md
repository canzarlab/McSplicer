

## Simulation study



We used <a target="_blank" href="https://github.com/alyssafrazee/polyester">Polyester</a> package to simulate reads from a human transcriptome (Homo\_sapiens.GRCh38.91) with abundances obtained from the RNA-experiment GEO accession GSM3094221 and computed by RSEM. Ground truth abundances are available in SRR6987574\_groundtruth.txt. The randomly selected list of 1000 genes is available in gene\_list.txt.


We ran generate_RNAseq_reads_polyester.R script to simulates single-end RNA-seq reads for our benchmark analyses.

**Input**: <br/> 1. Ground truth file: tab seperated file contains the trascript ids (column header should be transcript_id) and the corresponding abundance estimates including read counts. This scripts uses the read count (column index 5) and multiply it by scale parameter, however it can be modified to use TPM values instead. make sure to modify the code to match your file format. <br/>
2. cDNA fasta sequence file, you can download Homo_sapiens.GRCh38.91.cdna.fa from:
<a target="_blank" href="http://www.ensembl.org/Homo_sapiens">Ensembl (Homo sapiens)</a>.<br/>
3. Total number of reads. We generated 3 datasets with 20 million, 50 million and 75 million reads. <br/>
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
Tested on R version 3.5

To simulate our 50M short single-end RNA-seq reads dataset:
```
Rscript generate_RNAseq_reads_polyester.R ./SRR6987574_groundtruth.txt ./reference/Homo_sapiens.GRCh38.91.cdna.fa 50000000 ./out_dir/ FALSE 100 250 25 43
```
To generate the 20M and 75M RNA-seq reads datasets, you need to change the total number of reads while keeping all remaining args as they are.


Authors:
----------------------------
* Israa Alqassem (alqassem.isra@gmail.com)
* Derek Reiman (dreima2@uic.edu)


