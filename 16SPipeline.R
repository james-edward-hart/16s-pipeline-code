#!/usr/local/bin/Rscript

#16S Rotation Project Code


#Package Loading (Excluding Phangorn)
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(100)



#File Loading (“raw_data” is file name and it is in current wd))
miseq_path <- file.path("UMN_filtered")
filt_path <- file.path("UMN_filtered", "filtered")



#Seperate into forward and reverse reads
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]



#filter blank samples from data (ignore this chunk)
#fnFs <- fnFs[1:16]
#fnRs <- fnRs[1:16]



#print sample of 3 qaulity profiles (images save to current wd)
#ii <- sample(length(fnFs), 3)
#for(i in ii) {
#	jpeg(paste0('ForwardQualityProfile',i,'.png'))
#	print(plotQualityProfile(fnFs[i]) + ggtitle(“Fwd”))
#	dev.off()
#}
#for(i in ii) {
#	filename <- paste0('ReverseQualityProfile',i, '.png')
#	png(filename)
#	print(plotQualityProfile(fnRs[i]) + ggtitle(“Rev”))
#	dev.off()
#}



#Trim Sequences (in this case at 250 to 200, and first 10 nucleotides, with max 2 expected errors per read)
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
		      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft=10, truncLen=c(250, 200),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE)
}



#dereplicate fastq files and give derep objects names
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- gsub("_S[0-9]+_R1_001.fastq.gz", "", basename(filtFs))
names(derepFs) <- sam.names
names(derepRs) <- sam.names
save.image("16SAnalysis.rdata")


#first dada2 run to discover errors (should run as a nohup process; used full data set as it is small)
ddF <- dada(derepFs, err=NULL, selfConsist=TRUE)
ddR <- dada(derepRs, err=NULL, selfConsist=TRUE)



#plot errors for inspection (files to wd)
#png('ddFerrors.png')
#	plotErrors(ddF)
#	dev.off()
#png('ddRerrors.png')
#	plotErrors(ddR)
#	dev.off()



#second dada2 run to discover errors but given initial error matrix
dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE)
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE)
save.image("16SAnalysis.rdata")



#merge forward and reverse reads, excluding imperfect matches
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)



#make sequence table, remove chimeras
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
seqtab <- removeBimeraDenovo(seqtab.all)



#Assign Taxonomy to sequences using bayesian method classified training sets can be obtained here
#https://www.dropbox.com/sh/mfcivbudmc21cqt/AAB1l-AUM5uKvjrR33ct-cTXa?dl=0
ref_fasta <- "silva_nr_v123_train_set.fa"
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")



#Construct phylogenetic tree (*********update this once phangorn working********)
seqs <- getSequences(seqtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)



#properly format sample data
samdf <- data.frame(SampleID=rownames(seqtab))
rownames(samdf) <- samdf$SampleID
samdf$Group <- gsub("_.","",samdf$SampleID)
#This combines groups 1-4 (high BMI groups)
samdf[1:12,2]<-1




#Create Phyloseq object (without phylogenetic tree)
ps <- phyloseq(tax_table(taxtab), sample_data(samdf),otu_table(seqtab, taxa_are_rows = FALSE))
save.image("16SAnalysis.rdata")
ps



#create Phylum table
table(tax_table(ps)[, "Phylum"], exclude = NULL)



#create new phyloseq obj without unknown phylum features
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))



#Create a table with phylum taxonomies and feature prevalence
prevdf = apply(X = otu_table(ps0),
                 MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})



#Define low Prevalence features and filter them from the data (Specific to your data)
filterPhyla = c("Acidobacteria", "Candidate_division_OP3","Chloroflexi", "Cyanobacteria", "Gemmatimonadetes", "Gracilibacteria")
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1



#Plot Prevalence vs abundance
#png("prevelance_AbundanceFig.png")
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
#dev.off()



#Define a prevalence threshold and filter (in this case 0.0625 was used, for 16 samples each feature must appear once)
prevalenceThreshold = 0.062 * nsamples(ps0)
prevalenceThreshold
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)



#merge all taxa of the same genus
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)



#Generate a violin plot of abundance before and after transformation to relative abundance
#Here you can adjust the taxonomic level of the plots by changing Facet (ex. Facet="Class")
#png("TransformedAbundancePlots-CombinedGroups.png")
plot_abundance = function(physeq,title = "",
			     Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Group",y = "Abundance",
                                 color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
                position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
grid.arrange(nrow = 2, plotBefore, plotAfter)
#dev.off()



#subset and plot Clostridiales abundance plot by family to look for explanation
psOrd = subset_taxa(ps3ra, Order == "Clostridiales")
#png("ClostridialesAbundancePlot.png")
plot_abundance(psOrd, Facet = "Family", Color = NULL)
#dev.off()
save.image("16SAnalysis.rdata")
