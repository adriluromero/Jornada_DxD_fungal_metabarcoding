#Jornada_DxD_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#Section 4
#File name: 4_Taxonomic_tables
#Step: Use RStudio to process paired-end sequences and generate a taxonomic table.

#In this script, you will: 
# 1. Read-in output ITS2 sequences and perform quality filtering with dada2.
# 2. Assign taxonomy to the sequences using the most recent UNITE release.
# 3. Export the OTU and taxonomic tables. 


#Clear the environment
rm(list = ls(all.names = TRUE))

#Load the packages
require(dada2); packageVersion("dada2")
require(phyloseq)
require(ggplot2)
require(Biostrings)
require(vegan)
require(dplyr)
require(ape)
require(textshape)

#1. Read-in output ITS2 sequences and perform quality filtering with dada2.

#Set the path to where the 'fastq.gz' files are located.
pathF<- "/project/egcc/dada2/Read_1_ITSx"
pathR<- "/project/egcc/dada2/Read_2_ITSx"

#Read in the names of the fastq files.
pathF.names<- sort(list.files(pathF, pattern= ".fastq.gz", full.names = TRUE))
pathR.names<- sort(list.files(pathR, pattern= ".fastq.gz", full.names = TRUE))

#Extract sample names where everything before ".fastq.gz" is the sequence name.
pathF.sample.names <- sapply(strsplit(basename(pathF.names), ".fastq.gz"), `[`, 1)
pathR.sample.names <- sapply(strsplit(basename(pathR.names), ".fastq.gz"), `[`, 1)

#Make folders for filtered files.
filtFs1 <- file.path(pathF, "filtered", paste0(pathF.sample.names, "_filt.fastq.gz"))
filtRs1 <- file.path(pathR, "filtered", paste0(pathR.sample.names, "_filt.fastq.gz"))

#Fix names.
names(filtFs1) <- pathF.sample.names
names(filtRs1) <- pathR.sample.names

#Filter sequences.
out1<- filterAndTrim(pathF.names, filtFs1, pathR.names, filtRs1,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, minLen=100, verbose = TRUE, multithread = FALSE)
saveRDS(out1, "/project/egcc/dada2/Filtered/Filtered.seqs.rds")

#Sample '96_S96_131_L001_R1_001_filt.fastq.gz' & '96_S96_267_L001_R2_001_filt.fastq.gz' did not pass the filter. 

#Remove them from the read in file to get error rates. 
table(file.exists(filtFs1))
#FALSE  TRUE 
#    1   135

table(file.exists(filtRs1))
#FALSE  TRUE 
#    1   135

exists<- file.exists(filtFs1) & file.exists(filtRs1)

#Only keep the files that do exist. In this case, 135 sequences.
filtFs1<-filtFs1[exists]
filtRs1<-filtRs1[exists]

#Learn the error rates for both Forward (F) and Reverse (R) reads.
errF.1 <- learnErrors(filtFs1, multithread=FALSE)
#56657736 total bases in 241988 reads from 135 samples will be used for learning the error rates.

errR.1 <- learnErrors(filtRs1, multithread=FALSE)
#62864548 total bases in 241988 reads from 135 samples will be used for learning the error rates.

saveRDS(errF.1, "/project/egcc/dada2/Filtered/Forward.errors.rds")
saveRDS(errR.1, "/project/egcc/dada2/Filtered/Reverse.errors.rds")

#Examine sequence inference.
dada16S_F <- dada(filtFs1, err = errF.1, multithread=FALSE)
dada16S_R <- dada(filtRs1, err = errR.1, multithread=FALSE)
saveRDS(dada16S_F, "/project/egcc/dada2/Filtered/Forward.denoise.rds")
saveRDS(dada16S_R, "/project/egcc/dada2/Filtered/Reverse.denoise.rds")

#Merge forward and reverse reads.
mergers <- mergePairs(dada16S_F, filtFs1, dada16S_R, filtRs1, verbose=TRUE)
saveRDS(mergers, "/project/egcc/dada2/Filtered/Merged.Seqs.rds")

#Make a sequence table.
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/project/egcc/dada2/Filtered/Seq.tab.rds")

#Remove chimeras.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE) 

#Identified 83 bimeras out of 788 input sequences.

saveRDS(seqtab.nochim,"/project/egcc/dada2/Filtered/No.chimera.seq.tab.rds")

seqtab.nochim.length.filter <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 164:600]
saveRDS(seqtab.nochim.length.filter,"/project/egcc/dada2/Filtered/No.chimera.short.removed.seq.tab.rds")

#For this next step, remove row #133 which is sample 96_S96_131_L001_R1_001.fastq.gz from 'out1'
#this is because it has 0 reads out and cbind won't work with it. 

out2<- out1[-133,]

#How many sequences were retained at each QC step? 
getN <- function(x) sum(getUniques(x))
track <- cbind(out1, sapply(dada16S_F, getN), 
               sapply(dada16S_R, getN), 
               sapply(mergers,getN), 
               rowSums(seqtab.nochim),
               rowSums(seqtab.nochim.length.filter))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim", "length.filter")

View(track)

#The 'track' table shows that these samples had no sequences left after the merging step:
#104_S104_4_L001_R1_001.fastq.gz 
#126_S126_28_L001_R1_001.fastq.gz
#36_S36_50_L001_R1_001.fastq.gz

#2. Assign taxonomy to the sequences using the most recent UNITE release.

#Download the UNITE general FASTA release for Fungi from https://unite.ut.ee/repository.php 
#At the time of running this pipeline, the most recent UNITE release was created in 2021-05-10.
#This particular release was downloaded from https://doi.org/10.15156/BIO/1280049

#Upload and unpack the UNITE file (*.tgz) into your directory. 
#There will be a 'sh_general_release_dynamic_10.05.2021.fasta' and 'sh_general_release_dynamic_10.05.2021_dev.fasta' file. 
#Use the 'sh_general_release_dynamic_10.05.2021.fasta' file. 

unite.ref <- "/project/egcc/UNITE/sh_general_release_10.05.2021/sh_general_release_dynamic_10.05.2021.fasta"
taxa <- assignTaxonomy(seqtab.nochim.length.filter, unite.ref, multithread = FALSE, tryRC = TRUE)
#Running this step of assigning taxonomy can take from minutes to hours.

#Check the results. 
taxa.print <- taxa 
rownames(taxa.print) <- NULL 
head(taxa.print)

#which phyla are present?
unique(taxa.print[,2])

# [1] "p__Basidiomycota"          "p__Chytridiomycota"       
# [3]  NA                         "p__Mortierellomycota"     
# [5] "p__Ascomycota"             "p__Mucoromycota"          
# [7] "p__Glomeromycota"          "p__Olpidiomycota"         
# [9] "p__Calcarisporiellomycota"     

#3. Export the OTU and taxonomic tables. 
ps <- phyloseq(otu_table(seqtab.nochim.length.filter, taxa_are_rows=FALSE), 
               tax_table(taxa))

#Rename ASV sequences as ASV# 
#save as a string in case needed later-on as DNA
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 705 taxa and 135 samples ]
#tax_table()   Taxonomy Table:    [ 705 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 705 reference sequences ]

#Taxonomic filtering
rank_names(ps)
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

#Create a table with the number of features for each phyla

table(tax_table(ps)[, "Phylum"], exclude = NULL)

#            p__Ascomycota          p__Basidiomycota 
#                      343                       176 
#p__Calcarisporiellomycota        p__Chytridiomycota 
#                        1                        95 
#         p__Glomeromycota      p__Mortierellomycota 
#                        1                        21 
#          p__Mucoromycota          p__Olpidiomycota 
#                        1                         1 
#                     <NA> 
#                       65 

#Filter out potential artifacts
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "p__Calcarisporiellomycota", "p__Glomeromycota", "p__Mucoromycota", "p__Olpidiomycota"))

#Now, make these into data frames and save the files.
otu.table.not.rare<- data.frame(t(data.frame(otu_table(ps))))
#Look at the sequencing depth of each sample and sort column sums. 
OTU_Col_Sums <- sort(colSums(otu.table.not.rare)) 
tax.table<- data.frame(tax_table(ps))
identical(row.names(otu.table.not.rare), row.names(tax.table)) #row names match up perfect if TRUE
#[1] TRUE
otu.table.tax.table<- cbind(otu.table.not.rare, tax.table)
write.csv(otu.table.tax.table, "/project/egcc/dada2/otu.tax.table.csv")
write.csv(OTU_Col_Sums, "/project/egcc/dada2/otu.col.sums.csv")