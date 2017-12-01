                                        # Example analysis pipeline to obtain test statistics
                                        # and p.values for shape analysis
                                        # of read density patterns.

                                        # Cleaning the environment and setting up working directory:

rm(list=ls(all=TRUE)); current_wd=getwd(); local_wd=current_wd

                                        # Loading necessary packages:
library(data.table) #Fast datasets loading
library(fdasrvf) #The dynamic time warping functions
library(parallel) #enables parallelization for linux based machines
                                        #sourcing essential functions
source("aux_functions.R")
                                        #loading the GTF file containing exon/region information
exons<-fread("exons.gtf") 
in.bam.headers=which(exons$chrom %in% paste("chr", c(1:22, "X", "Y", "M"), sep=""))
exons=exons[in.bam.headers];  genelist=exons[[1]]
                                        #Setting minimal number of reads per exon/region
read.cutoff=20 

#######################################################################
                                        # Script for got Obtaining the reads per exon
                                        # The script focuses on obtaining R objects "all.reads" and "all.counts", 
                                        # using the efficient samtools functions incorporated in R via package Rsamtools.
                                        # "all.reads" is a list of genes. For each list element corresponding to gene 
                                        # there is another list - list of exons specific for the gene. 
                                        # Each element of this sublist contains the last subsublist correspnding to samples. 
                                        # Each element of subsublist stores read positions mapped to the exon 
                                        # in the selected sample.
                                        #In short:
                                        #list(all.genes)
                                        #-> sublist(exons)
                                        #->->subsublist(samples)
                                        #->->->vector (read positions, e.g. first/last/mid bp)

                                        # samotools accessible from R:
library(Rsamtools)
                                        # Path to main folder with datasets (change to your own):
bam_wd="~/works/Current/data_Nature_HOAXko/" ;
                                        # Folder names are corresponding to conditions
                                        # control and HOXA1 KO (change to your own):
conditions=c("control_A","control_B","control_C", "HOXA1_A", "HOXA1_B", "HOXA1_C" )
condition.labels<-c(1,1,1,2,2,2) 
                                        # name of the sorted alignet bam file, outputed by e.g. bowtie or tophat. 
                                        # Note that the scripts assumes that all files in both conditions and 
                                        # all samples are named exactly the same. The files and conditions are distinguished 
                                        # only by the folder name/path:
genericBAMfile="reads.sorted.bam" 
                                        # Filtering reads that map only to known chromosomes
                                        # using Rsamtool
what <- c("pos","mapq")
all.counts=list(); all.reads=list() 
gene.range=c(1:length(genelist)) 
                                        #getting reads mapped to genes. see get.read() function
  all.reads <- mclapply(gene.range, function(gene) get.reads(gene, output = "reads"), 
                        mc.cores = detectCores()-1, mc.preschedule=FALSE)
names(all.reads)=genelist[gene.range]

                                        #getting counts of reads mapped to gene using mclapply.
                                        #See the get.reads() function
all.counts <- mclapply(gene.range, function(gene) get.reads(gene, output = "counts"), 
                       mc.cores = detectCores()-1, mc.preschedule=FALSE)
names(all.counts)=genelist[gene.range]


                                        #As the above steps can be quite lengthy,
                                        #we suggest saving the results for future use.
save(all.reads, file="all.reads.Rdata")
save(all.counts, file="all.counts.Rdata")
                                        #For the rest of the code to work,
                                        #the results can be accessed by:
load("all.reads.Rdata")
load("all.counts.Rdata")



#############################################################################
                                        #main analysis for exon and gene level normalization (choose one)
                                        #The following script assumes that the reads positions per exon for each gene are known
                                        # and stored in all.reads and all.counts R objects in the following format
                                        #"all.reads" is a list of genes. 
                                        # For each list element corresponding to gene there is another list 
                                        # - list of exons specific for the gene. 
                                        # Each element of this sublist contains the last subsublist correspnding to samples. 
                                        # Each element of subsublist stores read positions mapped to the exon in the selected sample.

                                        #In short:
                                        #list(all.genes)
                                        #-> sublist(exons)
                                        #->->subsublist(samples)
                                        #->->->vector (read positions, e.g. first/last/mid bp)

                                        #I you don't have the follwoing format
                                        #the previous script above, gives an example
                                        #on how to obtain those from a sorted bam files.

                                        #loading the all.reads and all.counts R objects.
load("all.reads.Rdata")
load("all.counts.Rdata")

gtf.rows=c(1:length(all.reads))
                                        #1. Exon level analysis
#----------------------------------------------------------------------------
print(paste("CALCULATING STATISTICS", Sys.time(), sep=" "))
stats<-mclapply(gtf.rows, function(i) main.analysis(i, condition.labels) , 
                mc.preschedule = F, mc.cores = detectCores()-1)
names(stats) = genelist[gtf.rows]

print(paste("SAVING RESULTS in ",getwd(),"stats.Rdata | ", Sys.time(), sep=""))
save(stats, file="stats.Rdata")

#----------------------------------------------------------------------------
                                        #2. Gene level analysis (whole spliced genes compared)
print(paste("CALCULATING STATISTICS", Sys.time(), sep=" "))
stats.genes<-mclapply(gtf.rows, function(i)  main.analysis.gene.level(i, condition.labels),
                      mc.preschedule = F, mc.cores = detectCores()-1)
names(stats.genes) = genelist[gtf.rows]

print(paste("SAVING RESULTS in ",getwd(),"stats.genes.Rdata | ", Sys.time(), sep=""))
save(stats.genes, file="stats.genes.Rdata")

#############################################################################
                                        #Viewing the results of the analysis.
load("stats.Rdata") 
load("stats.genes.Rdata")
head(stats)
head(stats.genes)
                                        #Obtaining the list of DE genes on given significance level:
threshold=0.05

stat.pv=unlist(lapply(stats.genes, function(x) x[2])); 
idx = which(stat.pv >0 & stat.pv<threshold);
DE.genes = unlist(lapply(  names(idx), function(nam) substr(nam,1,10)))

statistics.pv=colnames(stats[[1]])[seq(2,6,2)]
DE.list<-c(lapply(statistics.pv, function(statname) get.DE(statname)),list(DE.genes))
names(DE.list)<-c("L0", "L0 Aligned", "L2 Aligned" ,"L0 gene level" )

DE.list

#############################################################################
                                        #Plotting particular gene (exon level analysis only)

check.now=DE.list[[3]];
which.model=3
                                        #Loop to plot all genes in "check.now" vector
for (gene.name in check.now)
{ 
  gene=which(genelist==gene.name)
  print(gene)
  current.gene=all.reads[[gene]]
  check=na.pass(stats[[gene]][,2*(c(1,2,3)[which.model[1]])]);
  exon.idx=which(check>0 & check<threshold & apply(all.counts[[gene]],1,min)>read.cutoff)
  
  for (id.ex in exon.idx)
  {
    plot.exons(current.gene[[id.ex]], gene.name = gene.name, id.ex = id.ex, prefix="test", condition.labels=condition.labels)
  }
}  

##################################################################################



