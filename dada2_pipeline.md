library(dada2)
library(ggplot2)
library(phyloseq)
library(scales)

setwd("C:/Users/lexim/Desktop/Binf_F22_R/Assignment3")

#### working directory is set as object called "path", which is used as a variable in the following code
path <- "C:/Users/lexim/Desktop/Binf_F22_R/Assignment3"


#### list all sequencing files
list.files(path)


#### Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#make an object with forward reads (fnFs) and reverse reads (fnRs)
#this can be changed if reads come in a different format
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


#### Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### Now we're going to inspect read quality profiles
#### Start by visualizing quality profiles of the forward & reverse reads
plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

#### visually look at the quality of the reads
#### 20% is about the expected cutoff point for quality

#### WHAT TO TRUNCATE DEPENDS ON THE SAMPLE
#### If you are using a less-overlapping primer set, like V1-V2 or V3-V4, your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them

#### Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#### filter and trim reads
#### pay attention to truncation length (truncLen)
#### THIS STEP IS ALWAYS DIFFERENT DEPENDING ON READS
#### if things aren't working later on, try changing the following parameters

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, trimLeft = c(18,20)) # On Windows set multithread=FALSE
head(out)
out

#### this step takes a long time
errF <- learnErrors(filtFs, multithread=FALSE, MAX_CONSIST=20) #for my own dataset, I'd like to add MAX_CONSIST=20 if it doesn't find convergence at 10
errR <- learnErrors(filtFs, multithread=FALSE, MAX_CONSIST=20)

#### you can do this instead of graphing
dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)

#### Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#### Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### Sample interference
#### this is where it picks apart unique sequences
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

#### you can switch to a server around this point if you need to !!

#### inspect the data to see what's up
dadaFs[[1]]

#### merge paired end reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#### Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### construct sequence table with ASVs
#### probably save this table as a csv or something for my own data
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#### Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#### remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#### The seqtab.nochim is a table with ASV seq and the count of how many times it occured. We would like to transpose it to be handy for later ( so that ASVs are in rows, and samples are in headers)
#### This transposes the seqtab.nochim data if you want to look at it as a column
flipped_seqtab.nochim<- as.data.frame(t(seqtab.nochim))

#### Can now just look at the reads that made it through each step in our pipeline.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#### If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#### Assign taxonomy finally, woo! You have to download your reference database! Always use the most recent. I like the Silva. 
#### can be downloaded here: https://benjjneb.github.io/dada2/training.html

#### TROUBLESHOOTING - make sure you have placed the file in the path, sometimes R studio will read this from where it is installed

#### keep this database file here
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#### SAVING AS YOU GO 
#### we should save our taxa file ( has our ASVs and their ID), and the seqtab.nochim file which lists our seqs without chimeras. Also we should save the transposed seqtab.nochim.

write.csv(taxa, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/L520_taxa.csv")
write.csv(seqtab.nochim, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/L520_seqtab_nochim.csv")
write.csv(flipped_seqtab.nochim, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/L520_flipped_seqtab_nochim.csv")

#### more convenient to have the ASV seq and the count in one sheet
#### this saves your flipped seqtab no chim file with your taxa data as one data sheet

OTUabund<-cbind(flipped_seqtab.nochim, taxa)
write.csv(OTUabund, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/L520_OTUabund2.csv")

#### note to self: we can recall already processed files and continue on from here:)


#### Time to PLOT !!

#### load these packages for plotting

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")


#### construct a dataframe from our file names
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(samples.out)
rownames(samdf) <- samples.out


#### here we're making a phyloseq object using dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

#### Make a table that includes ASVs too, because that's what we care about
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#### the most basic plot that will be our base
p=plot_bar(ps)
p

#### now we're colouring by phylum, which is better
p + geom_bar(aes(fill=Phylum), stat="identity")

#### the previous plot displays absolute abundance, but we want relative abundance 
relative<- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
relative

#### graph relative abundance by phylum
Phylum_graph <- plot_bar(relative, fill="Phylum") +
  ylab ("Relative Abundance")

Phylum_graph

#### this is prettier
Phylum_graph + 
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  labs(title="Microbial community composition of \n Lake 520 permafrost (16S rRNA gene)")

#### save environment!
saveRDS(ps, "L520_amplicon_ps.rds")

#### Tax_glom will conglomerate all the reads of an identical taxa together ( here, we want to count everything that is the same phylum)
ps_phylum <- tax_glom(ps, "Phylum")

#### Now we are just taking this, and turning it into relative abundance as we did above
ps1_phylum_relabun <- transform_sample_counts(ps_phylum, function(ASV) ASV/sum(ASV))

#### The psmelt function of phyloseq will make a dataframe of our phyloseq data for us, we also need it to be a factor
taxa_abundance_table_phylum <- psmelt(ps1_phylum_relabun)

### save this file
write.csv(taxa_abundance_table_phylum, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/taxa_abundance_table_phylum.csv")


#### you can open this file and go from here if needed!
#### taxa_abundance_table_phylum <- read.csv("taxa_abundance_table_phylum.csv")

taxa_abundance_table_phylum$Phylum<-factor(taxa_abundance_table_phylum$Phylum)

#### we want phyla with no sequences to not have a dot
taxa_abundance_table_phylum[taxa_abundance_table_phylum == 0] <- NA


#### plot relative abundance of phyla as columns
ggplot(data=taxa_abundance_table_phylum,mapping=aes(x=Sample,y=Abundance*100,)) +
  geom_col() +
  aes(fill=Phylum) + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="Microbial community composition of Lake 520 permafrost, by phyla (16S rRNA gene)", 
       x="Samples", 
       y="Phylum Relative Abundance (%)")  + 
  theme(axis.text.x = element_text(face="bold", size=10, angle=90),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.title.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold", size=10, angle=0),
        axis.title.y = element_text(face="bold")) 

#### Plot by order - we'll do the same thing as abpve but want to count everything that is the same order
ps_order <- tax_glom(ps, "Order")
ps1_order_relabun <- transform_sample_counts(ps_order, function(ASV) ASV/sum(ASV))
taxa_abundance_table_order <- psmelt(ps1_order_relabun)
write.csv(taxa_abundance_table_order, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/taxa_abundance_table_order.csv")
taxa_abundance_table_order$Order<-factor(taxa_abundance_table_order$Order)
taxa_abundance_table_order[taxa_abundance_table_order == 0] <- NA

#### order bar plot
ggplot(data=taxa_abundance_table_order,mapping=aes(x=Sample,y=Abundance*100,)) +
  geom_col() +
  aes(fill=Order) + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="Microbial community composition of Lake 520 permafrost, by order (16S rRNA gene)", 
       x="Samples", 
       y="Order Relative Abundance (%)")  + 
  theme(axis.text.x = element_text(face="bold", size=10, angle=90),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.title.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold", size=10, angle=0),
        axis.title.y = element_text(face="bold"))

#### I want to see if there's any overlap in ASVs between my sample and control and I think I want to look at this at the genus level. Produce a table identified to the genus level, and then sort by abundance to see if there are any ovelapping ASVs

ps_genus <- tax_glom(ps, "Genus")
ps1_genus_relabun <- transform_sample_counts(ps_genus, function(ASV) ASV/sum(ASV))
taxa_abundance_table_genus <- psmelt(ps1_genus_relabun)
write.csv(taxa_abundance_table_genus, file="C:/Users/lexim/Desktop/Binf_F22_R/Assignment3/taxa_abundance_table_genus.csv")

#### The plots before were super basic, but this is the actual plot we want to use in publications and at conferences. Plot rel abundance of phyla as abundant-dependent points
ggplot(data=taxa_abundance_table_phylum,mapping=aes(x=Sample, y=Phylum)) +
  geom_point() +
  (aes(size = Abundance, colour=Phylum)) +
  labs(title="Microbial community composition of Lake 520 permafrost by phyla", x="Samples", y="Phylum Relative Abundance (%)", size = "Relative Abundance(%)")  +
  scale_size_area() +
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(face="plain", color="Black", size=12, angle=90),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.text.y = element_text(face="bold", color="Black", size=10, angle=0)) 
