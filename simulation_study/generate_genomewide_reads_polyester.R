library(polyester)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 9 ) {
  
  stop(paste('You have to provide the following arguments:',
             '1. A file containing a transcript abundances',
             '2. cDNA fasta sequence file',
             '3. Scale for the initial read counts',
             '4. Output directory',
             '5.TRUE or FALSE whether to simulate paired-end reads or single-end reads.',
             '6. Read length',
             '7. Mean fragment length',
             '8. The standard deviation of the fragment length',
             '9. seed value to make the results reproducable', sep='\n'))
} 


#/home/isra/algbio/sophia/simulated_data/DTU_single_cutoff/fullmodel/simulated_reads/SRR6987574_sample1
abundance_file <- args[1]  # file containing tpm values
fasta_file <- args[2] 
scale <-  as.numeric(args[3])
out_dir <- args[4]
paired_or_single <-  args[5]  # paired or single end reads
read_length <- strtoi(args[6])
mean_frag_len <- strtoi(args[7])
sd_frag_len <- strtoi(args[8])
seed_val <- strtoi(args[9])

#abundance_file <- '/home/isra/algbio/sophia/simulated_data/DTU_single_cutoff_pablo/fullmodel/SRR6987574_fullmodel_groundtruth.txt'
#fasta_file = '/home/isra/algbio/sophia/reference/Homo_sapiens.GRCh38.91/Homo_sapiens.GRCh38.91.cdna.fa'
#scale <- 2
#out_dir <- '/home/isra/algbio/israa/genome_simulation_polyester/rsim_reads_test2/'
#paired_or_single <- FALSE
#read_length <- 100 
#mean_frag_len <- 250 
#sd_frag_len <- 25
#seed_val <- 40


fasta = readDNAStringSet(fasta_file)

# You need to change the following block depends on the format of your abundance file
##########################################################################################################################
## Read transcript id and its corresponding abundance as a key value list to reduce the search time

tx_abundance_list = read.csv(abundance_file, sep=" ")
trans_count_list <- vector(length=nrow(tx_abundance_list))
names(trans_count_list) <- tx_abundance_list$transcript_id     ## Define the transcript ids, i.e., target names as a list names

print("Reading abundances.....")

pb_perc <- txtProgressBar(min = 0, max = nrow(tx_abundance_list), style = 3) # perecntage
pb_bar <- txtProgressBar(min = 0, max = nrow(tx_abundance_list), style = 1)  # visual

for (i in  1:nrow(tx_abundance_list)){
  target_id <- toString(tx_abundance_list[i,2])     ## in the abundance file, the 2nd colum has the transcript ids
  est_count <- as.numeric(tx_abundance_list[i,13])  ##
  trans_count_list[[target_id]] <- est_count
  setTxtProgressBar(pb_perc, i)
  setTxtProgressBar(pb_bar, i)
}

close(pb_perc)
close(pb_bar)
##########################################################################################################################

print("Filtering out transcripts not found in the abundance file.....This may take a while!")

fasta_data <- strsplit(names(fasta), " ") #split fasta data by transcript_id gene_name
readspertx_list <- vector(length=0)
id_list <- vector(length=0)   # to keep track of the indices of transcript in fasta file ( not all tx in fasta file exist in the abundance file)
total_reads <- 0
#count <- 0

for (i in 1:length(fasta_data)){
  
  target_id <- fasta_data[[i]][1]
  
  if (!(target_id  %in% names(trans_count_list))){  # skip the transcript which was not found in the abundance file
    #print(paste (target_id,'not found in the abundance file', sep = " "))
    next 
  }
  
  est_count <- trans_count_list[[target_id]]
  X <- round(est_count * as.numeric(scale)) 
  if (X == 0) {
    next  # Discard transcripts where estimated counts is 0, no. reads for such transcript is 0.
  }
  
  id_list<-append(id_list, i)
  total_reads <- total_reads + X
  readspertx_list <- append(readspertx_list, X)
  
  #print(paste(target_id,X,sep=':'))
  #count <- count + 1
  #if (count >= 100){
  # break
  #}
}

fasta_new <- fasta[id_list] 
writeXStringSet(fasta_new, 'fasta_new.fa')
fold_changes = matrix(rep(1), nrow=length(fasta_new), ncol=2) 
dir.create(out_dir, showWarnings=FALSE)
gc()  # this is a garbage collector, saw this in Derek's code and thought it's cool to have it

print("Started simulation.....")
print(paste('Total number of reads to generate = ',total_reads,sep=''))
old <- Sys.time()
simulate_experiment('fasta_new.fa',reads_per_transcript = readspertx_list,num_reps=c(1,1),
                    readlen = read_length,distr='normal',fraglen = mean_frag_len,fragsd = sd_frag_len,
                    fold_changes=fold_changes, outdir=out_dir, paired=paired_or_single, seed = seed_val) 
new <- Sys.time() - old 
print(new)
print(paste("Done! Output is written to",out_dir,sep=" "))

