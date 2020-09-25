## -------- simulate RNA-seq reads using polyester

library(polyester, quietly = T, warn.conflicts = F)
library(Biostrings, quietly = T, warn.conflicts = F)

#Rscript generate_genomewide_reads_polyester.R /data/israa/mcsplicer_paper_genome_simulation/SRR6987574_fullmodel_groundtruth.txt /data/sophia/reference/Homo_sapiens.GRCh38.91/Homo_sapiens.GRCh38.91.cdna.fa 10000000 /data/israa/mcsplicer_paper_genome_simulation/10M_reads/ FALSE 100 250 25 43
  

## Args -----------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 9 ) {
  
  stop(paste('You have to provide the following arguments:',
             '1. A file containing a transcript abundances',
             '2. cDNA fasta sequence file',
             '3. Total number of reads',
             '4. Output directory',
             '5. TRUE or FALSE whether to simulate paired-end reads or not.',
             '6. Read length',
             '7. Mean fragment length',
             '8. The standard deviation of the fragment length',
             '9. Seed value to make the results reproducable', sep='\n'))
} 

abundance_file <- args[1]  # file containing tpm values
fasta_file <- args[2] 
no_reads <-  as.numeric(args[3])
out_dir <- args[4]
paired_end_reads <-  args[5]  # paired or single end reads
read_length <- strtoi(args[6])
mean_frag_len <- strtoi(args[7])
sd_frag_len <- strtoi(args[8])
seed_val <- strtoi(args[9])
## ---------------------------------------------------------------------------------------------

#abundance_file <- '/home/isra/algbio/israa/mcsplicer_paper_genome_simulation/SRR6987574_fullmodel_groundtruth.txt'
#fasta_file = '/home/isra/algbio/sophia/reference/Homo_sapiens.GRCh38.91/Homo_sapiens.GRCh38.91.cdna.fa'
#no_reads <- 10000000
#out_dir <- '/home/isra/algbio/israa/genome_simulation_polyester/Polyester_new_simulated_data/'
#paired_or_single <- FALSE
#read_length <- 100 
#mean_frag_len <- 250 
#sd_frag_len <- 25
#seed_val <- 43




## 1) Read and parse transcript abundance file; adjust according to your abundance file column names
abundances <-  read.csv(abundance_file, sep=" ")

# 1.a) get transcript ID and count columns from abundance file
transcript_counts <- abundances[, c("transcript_id", "isoform_count_gr1")]
colnames(transcript_counts)[colnames(transcript_counts)=="isoform_count_gr1"] <- "expected_count"
# 1.b) calculate scale factor
scale_factor <- no_reads/sum(transcript_counts$expected_count)
# 1.c) adjust expected counts (rounded to 2 decimal places) to give desired number of total reads
transcript_counts$expected_count <- round(transcript_counts$expected_count * scale_factor, 2)
# 1.d) generate a named vector containing the counts of transcripts not 0
transcript_counts_vector <- transcript_counts[!(transcript_counts$expected_count == 0), 2]
names(transcript_counts_vector) <- transcript_counts[!(transcript_counts$expected_count == 0), 1]
# 1.e) extract ids of transcripts that are to be simulated 
transcript_ids <- names(transcript_counts_vector)



## 2) Read fasta file
fasta <-  readDNAStringSet(fasta_file)
fasta_data <- strsplit(names(fasta), " ")




## 3) Calculate reads per transcript
message("Extracting counts of transcripts to simulate. This may take awhile.")
message("progress:")

# 3.a) initiate empty lists to add data to during for loop
total_reads <- 0
id_list <- vector()
counts_list <- vector()

# 3.b) generate a progress bar for this step since it takes a sig. amount of time
pb_bar <- txtProgressBar(min = 0, max = length(fasta_data), style = 1)
pb_perc <- txtProgressBar(min = 0, max = length(fasta_data), style = 3)

for (i in 1:length(fasta_data)){
  # get name of transcript to look for
  target_id <- fasta_data[[i]][1]
  
  # skip the transcript if it is not found in the abundance file or of its count equals to 0
  if(!(target_id %in% transcript_ids)){
    next 
  }
  
  count <- transcript_counts_vector[[target_id]]
  
  id_list <- append(id_list, i)
  counts_list <- append(counts_list, count)
  total_reads <- total_reads + count
  
  # update progress bar
  setTxtProgressBar(pb_bar, i)
  setTxtProgressBar(pb_perc, i)
}

# 3.c) close the progress bar
close(pb_perc)
close(pb_bar)

# 3.d) subset fasta file to only get the entries we need 
fasta_subset <- fasta[id_list]
writeXStringSet(fasta_subset, filepath= paste0(out_dir, "/fasta_simulation.fa"))

#collect garbage to free up memory after long loop
gc()

## 4) Finally, simulate reads
message("Starting simulation of sample")
message('Total number of reads to generate: ',total_reads)

# 4.a) start simulation time
old <- Sys.time()

simulate_experiment_countmat(fasta = paste0(out_dir, "/fasta_simulation.fa"),
                               readmat = as.matrix(counts_list),
                               outdir = out_dir,
                               paired = paired_end_reads,
                               seed = seed_val,
                               readlen = read_length,
                               distr='normal',
                               fraglen = mean_frag_len,
                               fragsd = sd_frag_len,
                               error_model = "uniform")


# 4.b) end simulation time
new <- Sys.time() - old 

message("Simulation was done in ", new)
print(paste("Output is written to",out_dir,sep=" "))
