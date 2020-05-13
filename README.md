McSplicer
=========
**McSplicer** is a probabilistic model for quantifying splicing processes, rather than modeling an individual outcome of a process such as exon skipping. McSplicer is based on gene-wide usage of splice sites. We assume that potential 5' and 3' splice sites are given. This information can be obtained from annotation databases or estimated from RNA-seq data by using existing assemblers.  McSplicer estimates splice site usages using an inhomogeneous Markov chain, hence the name **M**arkov **c**hain **Splicer**. Using these parameters, we can describe the splicing process, and estimate the probabilities of various local splicing events.


   ![McSplicer](https://github.com/canzarlab/McSplicer/blob/master/Figures/McSplicer_summary.png) 



Using this software
-------------------

## Installation<a name="installation"></a>

```shell
git clone https://github.com/shimlab/McSplicer.git
```

You can execute the script:

```shell
python2 ./python_code/McSplicer.py --help
```


### Dependencies<a name="dependencies"></a>

McSplicer was implemnted and tested on python 2.7, and requires only few standard packages:
- numpy>=1.13.1
- pandas>=0.20.3

## Usage <a name="usage"></a>


Execute McSplicer script with `--help` option for a complete list of options.  

Sample data and usage examples can be found at `examples` subfolder.

Initially, you need as inputs the annotation in GTF format and aligned and indexed RNA-seq reads in bam file, then you can run McSplicer easily in 3 steps:

1. Run exonRefine to refine the set of exons into non-overlapping segments.

```shell
./bin/exonRefine  <annotation.gtf> --prefix OUTPUT_PREFIX
```

2. Run SigCount to parse short RNA-seq read alignments and generate signature counts, McSplicer works only with single-end reads at the moment, so you need to convert paired-end signuture counts to single-end counts.

    * For single-end reads:

               
               ./bin/sigcount_se <alignments.bam> <annotation_refined.gtf> <outfile-prefix>
               
  
    * For paired-end reads:

               
               ./bin/sigcount_pe <alignments.bam> <annotation_refined.gtf> <outfile-prefix>
	            python2 ./python_code/pe2seCnt.py <pe_signature_count.cnt> > <se_signature_count.cnt>
               
	      
	       
		
3. run McSplucer to get splice site usage estimates.

```shell
python2 ./python_code/McSplicer.py \
		--gtf REFINED_GTF \
		--count_file SE_SIGNATURE_COUNT \
		--gene_id GENE_ID\
		--out_dir OUTPUT_DIRECTORY\
		--bootstraps NUM_OF_BOOSTRAPS\
		--read_len READ_LENGTH\
		--prefix OUT_FILE_PREFIX

  ```




McSplicer developers:
----------------------------
* Israa Alqassem (alqassem.isra@gmail.com)
* Yash Kumar Sonthalia (yks01247@gmail.com)


&copy; 2017 McSplicer





