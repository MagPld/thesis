if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

setwd("~/examensarbete/test_scripting")

library(dada2)
list.files()

samples <- scan("samples", what="character")

forward_reads <- paste0(samples,".R1.fastq")
reverse_reads <- paste0(samples,".R2.fastq")

filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")

plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)


start_time <- Sys.time()
plotQualityProfile(reverse_reads)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              truncQ=2, rm.phix=FALSE, truncLen = c(250,200))
end_time <- Sys.time()
end_time - start_time


plotQualityProfile(filtered_forward_reads)
plotQualityProfile(filtered_reverse_reads)

##Generate error model
start_time <- Sys.time()
err_forward_reads <- learnErrors(filtered_forward_reads)
end_time <- Sys.time()
end_time - start_time

plotErrors(err_forward_reads, nominalQ=TRUE)

err_reverse_reads <- learnErrors(filtered_reverse_reads)
plotErrors(err_reverse_reads, nominalQ=TRUE)

##Dereplication
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

##Infer ASVs
start_time <- Sys.time()
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo")
end_time <- Sys.time()
end_time - start_time

dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo")


merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE)


seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)


seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
sum(seqtab.nochim)/sum(seqtab)  #0.98


# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab


write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)



## downloading DECIPHER-formatted SILVA v138 reference
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

## loading DECIPHER
library(DECIPHER)
packageVersion("DECIPHER") # v2.6.0 when this was initially put together, though might be different in the binder or conda installation, that's ok!

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)


# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)


BiocManager::install("decontam")
library(decontam)
packageVersion("decontam") # 1.1.2 when this was put together


BiocManager::install("DESeq2")
library(DESeq2) ; packageVersion("DESeq2")




write(asv_fasta, "ASVs.fa")
write.table(asv_tab, "ASVs_counts.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax, "ASVs_taxonomy.tsv",
            sep="\t", quote=F, col.names=NA)


count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, check.names = F, sep="\t")

tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")

# and setting the color column to be of type "character", which helps later
sample_info_tab$color <- as.character(sample_info_tab$color)

sample_info_tab # to take a peek

#sample_info_tab_fix <- sample_info_tab[nrow(sample_info_tab):1, ]
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~type) 

deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
