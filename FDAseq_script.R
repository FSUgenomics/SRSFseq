                                        # Example analysis pipeline to obtain test statistics and p.values
                                        #for shape analysis of read density patterns.

                                        #Cleaning the environment and setting up working directory:

rm(list=ls(all=TRUE)); current_wd=getwd(); local_wd=current_wd

                                        #Packages:
library(data.table) #Fast datasets loading
library(fdasrvf) #The dynamic time warping functions
library(Rsamtools) #samotools accessible from R
library(parallel) #enables parallelization for linux based machines


                                        #loading essential files and functions
source("aux_functions.R")
                                        #The GTF file containing exon information
exons<-fread("exons.gtf") 
                                        #Setting minimal numbe of reads per exon
read.cutoff=20 
                                        #Path main folder with datasets:
bam_wd="~/works/Current/data_Nature_HOAXko/" ;
                                        #Folder names corresponding to conditions:
conditions=c("control_A","control_B","control_C", "HOXA1_A", "HOXA1_B", "HOXA1_C" )
condition.labels<-c(1,1,1,2,2,2) #This needs to be fixed:D

                                        #name of the sorted alignet bam file, outputed by e.g. bowtie or tophat. Note that the scripts assumes that all files in both conditions and all samples are named exactly the same. The files and conditions are distinguished only by the folder name:
genericBAMfile="reads.sorted.bam" 

                                        # Filtering reads that map only to known chromosomes
in.bam.headers=which(exons$chrom %in% paste("chr", c(1:22, "X", "Y", "M"), sep=""))
exons=exons[in.bam.headers];  genelist=exons[[1]]
what <- c("pos","mapq")
all.counts=list(); all.reads=list() #Initializing output variables
gene.range=c(1:length(genelist))[1:100] #Range of genes we want to scan - subest of all 1:length(genelist). Default is all, that is 1:length(genelist)

                                        #getting reads mapped to genes
  all.reads <- mclapply(gene.range, function(gene) get.reads(gene, output = "reads"), mc.cores = detectCores()-1, mc.preschedule=FALSE)
  names(all.reads)=genelist[gene.range]


#getting counts of reads mapped to gene using mclapply. See the get.reads() function
  all.counts <- mclapply(gene.range, function(gene) get.reads(gene, output = "counts"), mc.cores = detectCores()-1, mc.preschedule=FALSE)
  names(all.counts)=genelist[gene.range]

#As the above steps can be quite lengthy, we suggest saving the results for future use. For the rest of the code to work, the results can be accessed by:

  save(all.reads, file="all.reads.Rdata")
  save(all.counts, file="all.counts.Rdata")

load("all.reads.Rdata")
load("all.counts.Rdata")



#############################################################################
                                        #main analysis for exon and gene level normalization (choose one)


                                        #1. Exon level normalization analysis
gtf.rows=c(1:length(all.reads))

print(paste("CALCULATING STATISTICS", Sys.time(), sep=" "))
stats<-mclapply(gtf.rows, function(i) main.analysis(i, condition.labels) , mc.preschedule = F, mc.cores = detectCores()-1)
names(stats) = genelist[gtf.rows]

print(paste("SAVING RESULTS in ",getwd(),"stats.Rdata | ", Sys.time(), sep=""))
save(stats, file="stats.Rdata")

#----------------------------------------------------------------------------
                                        #2. Gene level analysis (whole spliced genes compared)
print(paste("CALCULATING STATISTICS", Sys.time(), sep=" "))
stats.genes<-mclapply(gtf.rows, function(i)  main.analysis.gene.level(i, condition.labels), mc.preschedule = F, mc.cores = detectCores()-1)
names(stats.genes) = genelist[gtf.rows]


  print(paste("SAVING RESULTS in ",getwd(),"stats.genes.Rdata | ", Sys.time(), sep=""))

  save(stats.genes, file="stats.genes.Rdata")

#############################################################################                                        #Viewing the results of the analysis.
load("stats.Rdata") 
load("stats.genes.Rdata")
head(stats)
head(stats.genes)

                                        #Obtaining the list of DE genes on given significance level:
threshold=1

stat.pv=unlist(lapply(stats.genes, function(x) x[2])); 
idx = which(stat.pv >0 & stat.pv<threshold);
DE.genes = unlist(lapply(  names(idx), function(nam) substr(nam,1,10)))

statistics.pv=colnames(stats[[1]])[seq(2,6,2)]
DE.list<-c(lapply(statistics.pv, function(statname) get.DE(statname)),list(DE.genes))
names(DE.list)<-c("L0", "L0 Aligned", "L2 Aligned" ,"L0 gene level" )

#############################################################################                                        #Plotting particular gene (exon level analysis only)

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



