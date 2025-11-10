library(dada2)
library(ggplot2)
library(phyloseq)
library(Biostrings)
library(decontam)

## read files
path<-"data/raw"
#list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs<-sort(list.files(path, pattern="_1.fastq", full.name=TRUE))
fnRs<-sort(list.files(path, pattern="_2.fastq", full.name=TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[',1)

# read quality
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

## looks good. let's filter
#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)
#
# learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plotErrors(errF, nominalQ=TRUE)
#
# dada sample inference and merge reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

head(mergers[[1]])
#make sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

## remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# sanity check
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## assigning taxonomy to silva
taxa <- assignTaxonomy(seqtab.nochim, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# load metadata
metadata <- read.csv("data/dataframe_2023_clean.csv")

metadata$Sample_ID

metadata$Lab_ID

rownames(metadata) <- metadata$Lab_ID

rownames(seqtab.nochim)[!metadata$Lab_ID %in% rownames(seqtab.nochim)]

rownames(seqtab.nochim)[!rownames(seqtab.nochim)%in% metadata$Lab_ID]

nrow(seqtab.nochim)

## handover to phyloseq
PS <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(metadata),
               tax_table(taxa))

PS

all(sample_names(PS)==rownames(PS@sam_data))

summary(as.factor(PS@sam_data$State))

######### remove contaminants
### Removing contaminants
library("decontam")


PS_neg <- subset_samples(PS, grepl("NTC",rownames(PS@otu_table)))

PS@sam_data$Control <- FALSE


PS@sam_data$Control[which(sample_names(PS)%in%sample_names(PS_neg))] <- TRUE

# sanity check
PS@sam_data$Lab_ID[PS@sam_data$Control==FALSE]
rownames(PS@sam_data)[PS@sam_data$Control==TRUE]

## ----see-depths---------------------------------------------------------------
df <- as.data.frame(sample_data(PS)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(PS)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

ggplot(data=df, aes(x=Index, y=LibrarySize, color=Control)) + geom_point()

contamdf.freq <- isContaminant(PS, method="prevalence", neg="Control", threshold=c(0.1), normalize=TRUE)

table(contamdf.freq$contaminant)

### taxa to remove
PS@tax_table[rownames(contamdf.freq[contamdf.freq$contaminant==TRUE,]),5]

## let's remove them now and negative controls
Keep <- rownames(contamdf.freq[contamdf.freq$contaminant==FALSE,])
PS <- prune_samples(sample_data(PS)$Control == FALSE, PS)
PS <- prune_taxa(Keep, PS)


# sanity check
PS@sam_data$Lab_ID[PS@sam_data$Control==FALSE]
rownames(PS@sam_data)[PS@sam_data$Control==TRUE]

# removing ugly handlers
tax_table(PS)[, colnames(tax_table(PS))] <- gsub(tax_table(PS)[, colnames(tax_table(PS))], pattern="[a-z]__", replacement="")

PS

### quality filtering: removing low prevalent taxa (less than 5%)
x = phyloseq::taxa_sums(PS)

# remove singletons (prevalnce filter bellow 3%)
KeepTaxap <- microbiome::prevalence(PS)>0.03

PS <- phyloseq::prune_taxa(KeepTaxap, PS)

PS

summary(sample_sums(PS)) # library size

# adjust taxonomy
PS@tax_table[which(is.na(PS@tax_table[,2])),2] <- "Unknown_phylum"
tax <- as.data.frame(tax_table(PS))
tax[is.na(tax$Genus),]$Genus <- paste0("Unknown_genus_in_",tax[is.na(tax$Genus),]$Family)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Order)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Class)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Phylum)
PS@tax_table <-tax_table(as.matrix(tax))


saveRDS(PS, file="tmp/Phyloseq_f.RDS")
