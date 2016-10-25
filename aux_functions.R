
get.reads<-function(gene, output, qualityFilter=20)
{
  cat(paste(gene/length(genelist)*100 , "%", genelist[gene]),"\n")
  chr=exons[gene]$chrom
  Starts=as.numeric(unlist(strsplit(exons[gene]$exonStarts, ","))); 
  Ends=as.numeric(unlist(strsplit(exons[gene]$exonEnds, ",")));  
  exon.count=length(Starts)
  if (output=="reads")
  {
    gene.reads=list()
    for (exon in 1:exon.count)
    {
      which=RangesList(IRanges(Starts[exon], Ends[exon]))
      names(which)<-paste(chr)
      param <- ScanBamParam(which=which, what=what, mapqFilter = qualityFilter)
      exon.reads=list()
      for (c in 1:length(conditions))
      {
        bamfile=paste(bam_wd,conditions[c], "/",genericBAMfile, sep="")
        exon_bam=scanBam(bamfile, param = param)
        exon.reads[[c]]<-0;
        temp=exon_bam[[1]]$pos
        if (length(temp)>0) exon.reads[[c]]=na.omit(temp);
      }
      names(exon.reads)=conditions
      gene.reads[[exon]]=exon.reads 
    }
    names(gene.reads)=paste("exon",c(1:exon.count))
    return(gene.reads)
  }
  
  if (output=="counts")
  {
    gene.counts=array(0,c(exon.count,length(conditions)));
    rownames(gene.counts)<-paste("exon",c(1:exon.count));
    colnames(gene.counts)<-conditions
    for (exon in 1:exon.count)
    {
      which=RangesList(IRanges(Starts[exon], Ends[exon]))
      names(which)<-paste(chr)
      param <- ScanBamParam(which=which, what=what, mapqFilter = qualityFilter)
      gene.counts[exon,]=sapply(1:length(conditions), function(c) countBam(paste(bam_wd,conditions[c], "/",genericBAMfile, sep=""), param = param)$records)
    }
    return(gene.counts)
  }
}


main.analysis.gene.level<-function(gene, condition.labels)
{
  current.gene=all.reads[[gene]]
  
  temp<-gene.normalization(current.gene, all.counts[[gene]], genelist[gene])
  gene.gnormalized<-temp[[2]] 
  
  gene.stats=array(0,c(1,2))
  colnames(gene.stats)=c("L0.ANOVA","L0.pv")

if (length(gene.gnormalized)>0)
{ 
  gene.stats=array(0,c(1,2))
  colnames(gene.stats)=c("L0.ANOVA","L0.pv")

      # Setting up for intensity analysis
      domain=c(min(unlist(gene.gnormalized)), max(unlist(gene.gnormalized))) 
      samples=get.intensity(gene.gnormalized, domain); 
      samples.xy<-add.xy(samples, domain)
      #Calculating the statsitics and p.values
      gene.stats[1,]=get.statistics(samples.xy,domain, condition.labels)
}
  return(gene.stats)
}




get.DE<-function(statistics.pv)
{
  stat.names=names(stats)
  stat.pv=unlist(lapply(stats, function(x) x[,statistics.pv, drop=T])); 
  idx = which(stat.pv >0 & stat.pv<threshold); DE = unlist(lapply(  names(idx), function(nam) substr(nam,1,10)))
  return(DE)  
}


main.analysis<-function(gene, condition.labels)
{

  if (gene%%100==0) print(paste("GENE", gene, sep=" "))
  current.gene=all.reads[[gene]]
  exon.count=length(current.gene) # was guarding before against 0 length entries
  
  gene.stats=array(0,c(exon.count,6))
  colnames(gene.stats)=c("L0.ANOVA","L0.pv",
                         "L0aligned.ANOVA", "L0aligned.pv",
                         "L2aligned.ANOVA", "L2aligned.pv")
  rownames(gene.stats)=names(all.reads[[gene]])
  
  for (exon in 1:exon.count)
  {
    current.exon=current.gene[[exon]]  
    
    if (min(all.counts[[gene]][exon,]) >read.cutoff) 
    {
      # Setting up for intensity analysis
      domain=c(min(unlist(current.exon)), max(unlist(current.exon))) 
      samples=get.intensity(current.exon, domain); sample1<-samples[,which(condition.labels==1)];   sample2<-samples[,which(condition.labels==2)];
      # Adding SRVF curves and SRVF alignment versions of the curves
      sink("nowhere")
      km=get.means(samples,amp="L0");
      kmL2=get.means(samples,amp="L2");
      sink(NULL)
      fn=km$fn; fnL2=kmL2$qn
      
      curves.list=lapply(list(samples,fn,fnL2), function(curves) add.xy(curves, domain))
      
      #Calculating the statsitics and p.values
      gene.stats[exon,]=unlist(lapply(curves.list, function(curves.xy) get.statistics(curves.xy,domain, condition.labels)))
    }
    
  }
  return(gene.stats)
}

get.statistics<-function(curves.xy, domain, condition.labels)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------


  {
    curves.xy.1=curves.xy[which(condition.labels==1)];
    curves.xy.2=curves.xy[which(condition.labels==2)];
  
 
  
    mean1=get.means.L2(curves.xy.1,domain); mean2=get.means.L2(curves.xy.2,domain); 
    mean.global=get.means.L2(curves.xy,domain); 
    
    ss.total=get.var(mean.global, curves.xy, dom=domain)
    d=length(mean.global)/abs(diff(domain))
    ss.group= 3*sum((mean.global-mean1)^2/d) + 3*sum((mean.global-mean2)^2/d)
    ss.error=ss.total-ss.group
    return(c(ss.group*4/(ss.error*1), 1-pf(ss.group*4/(ss.error*1),1,length(condition.labels)-2)))
  
  
}

add.xy<-function(curves, domain)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
#               Simple poltting wrap-up. Uses mostly global variables 

{
  curves.xy=list()
  for (j in 1:dim(curves)[2])
  { temp.int=data.frame(   cbind(seq(domain[1],domain[2],length.out = length(curves[,j])),curves[,j])     )
  colnames(temp.int)=c('x','y')
  curves.xy[[j]]=temp.int}
  return(curves.xy)
}





gene.normalization<-function(current.gene, lengths, name)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
#               Simple poltting wrap-up. Uses mostly global variables 
{
  #lengths = lens[[gene]]; name=substr(filelist[gene],1,10)
  splice.idx=sort(unique(which(apply(lengths, 1, max)  > read.cutoff, arr.ind = T))) #originally was min instead of max
  if (length(splice.idx>0))
  {
    bed.row=exons[which(exons$alignID==name)]
    Starts=as.numeric(unlist(strsplit(bed.row$exonStarts, ","))); Ends=as.numeric(unlist(strsplit(bed.row$exonEnds, ",")));
    Starts=Starts[splice.idx]; Ends=Ends[splice.idx];
    shift=0;
    if (length(splice.idx) >1 ) shift=cumsum(c(0,Starts[2:length(Starts)] -Ends[1:(length(Ends)-1)]));
    
    spliced=apply(
              mapply(c,
                lapply(1:length(splice.idx), function(idx) 
                  lapply(current.gene[[splice.idx[idx]]], function (l.exon) l.exon-shift[idx])
                      ) 
                    ),1,function(x) unlist(x)
                 )
  
    domain=c(min(unlist(spliced)), max(unlist(spliced))) #Those are shifted coordinates, which need to be recalculated later
    gene.samples=get.intensity(spliced, domain, multiplier = 10*length(splice.idx));
    gene.samples<-add.xy(gene.samples, domain = domain)
    new.coords=cbind(Starts-shift, Ends-shift)
    new.current.exons=list()
    new.current.gene=list()
   
    for (ex in 1: dim(new.coords)[1])
      {
        for (samp.idx in 1:length(gene.samples))
        {
          
          temp.gene<-gene.samples[[samp.idx]]
          temp.exon.idx=which(temp.gene$x >new.coords[ex,1] & temp.gene$x < new.coords[ex,2])
          to.exon.list=data.frame(cbind(temp.gene$x[temp.exon.idx],temp.gene$y[temp.exon.idx]))
          colnames(to.exon.list)<-c('x','y')
          new.current.exons[[samp.idx]]=to.exon.list
        }
        new.current.gene[[ex]]=new.current.exons
      }
  #new.current.gene is already in xy coordinates for intensities
  #spliced  are just merged reads into spliced exons
  # gene.samples are merged intensities with xy coordinates
  return(list(new.current.gene, spliced, gene.samples))
  }
  return(NULL)
}


plot.exons<-function(current.exon, models.to.plot=c(1,2,3), gene.name,id.ex, prefix, condition.labels)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
#               Simple poltting wrap-up. Uses mostly global variables
{
    domain=c(min(unlist(current.exon)), max(unlist(current.exon)))
    samples=get.intensity(current.exon, domain);
    # Adding SRVF curves and SRVF alignment versions of the curves
    sink("nowhere")
      km=get.means(samples,amp="L0");
      kmL2=get.means(samples,amp="L2");
    sink(NULL)
    fn=km$fn; fnL2=kmL2$qn
    curves.list=lapply(list(samples,fn,fnL2), function(curves) add.xy(curves, domain))
    names(curves.list)<-c("ANOVA ", "ANOVA for aligned functions", "ANOVA for L2 aligned functions")
    
    #Plotting paramteres:      
    png(paste(local_wd,'/',prefix,gene.name,'_exons_',paste(id.ex, collapse=","),'.png',sep=""), width=2*1024, height=length(models.to.plot)*1024)
    layout(c(1,sort(rep(c(1:length(models.to.plot)+1),2))))
    fsize=10; PCH=1; CEX=5;  LWD=10;
    par(mar=c(10,15,10,12), mgp = c(15, 7, 0), bg = "white");
    
# The point patterns (NGS Data)  
    plot( current.exon[[1]], array(0, length(current.exon[[1]])) + rnorm(length(current.exon[[1]]),0,0.01),col=rgb(0,0,0,0,0),
         type="p", pch=PCH, cex=CEX,cex.lab=fsize-2, cex.axis=fsize-2, cex.main=fsize, cex.sub=fsize-2,
         xlab="Reference genome", ylab="", ylim=c(-0.1,.6), xlim=domain,main=paste("Expression pattern for: " ,gene.name ,id.ex))
    
    for (k in 1:length(condition.labels))
    {
      lines(current.exon[[k]], array(0, length(current.exon[[k]])) + rnorm(length(current.exon[[k]]),0,0.01)+(k-1)/10,type="p", pch=PCH, cex=CEX, col=rgb(condition.labels[k]-1,0,0,0.05), lwd=LWD)
    }
    legend("topright", legend=rev(unlist(lapply(current.exon, function(l) length(l))))[1:length(condition.labels)], cex=CEX-1, col=rev(condition.labels), pch=c(rep(PCH,length(condition.labels))) )

# Transformed intensities for analysis  
      lapply(models.to.plot, function(m.id) plot.curves(m.id, domain, curves.list[m.id], id.ex = id.ex, condition.labels=condition.labels))
            
dev.off();
}

plot.curves<-function(model.idx, domain, curves, id.ex, condition.labels)
{
  fsize=10; PCH=1; CEX=6;  LWD=10;
  par(mar=c(10,15,10,12), mgp = c(15, 7, 0), bg = "white");
  title<-names(curves)
  curves<-curves[[1]];
  curves1<-curves[which(condition.labels==1)];
  curves2<-curves[which(condition.labels==2)];
  #Cosmetic adjustments of ylim for plots
  mean1.L2=get.means.L2(curves1,domain); mean2.L2=get.means.L2(curves2,domain); 
  yrange=c(1.2*(-abs(min(mean1.L2, mean2.L2))),1.2*max(mean1.L2, mean2.L2 ) )
  plot(curves1[[1]], col=1, type="l",lty=1, lwd=LWD , main=title, pch=PCH, cex=CEX,cex.lab=fsize-2, cex.axis=fsize-2, cex.main=fsize, cex.sub=fsize-2 ,xlab="",ylim=yrange,ylab="" ) 
  lapply(curves1 ,function(dat) lines(dat, col=1 ,lty=1, lwd=LWD));
  lapply(curves2 ,function(dat) lines(dat, col=2 ,lty=1, lwd=LWD)); 
  legend("top", legend=paste("P.value = ", round(stats[[gene]][id.ex,(c(1,2,3)[model.idx])*2],4),sep=""), cex=CEX)
}


plot.gene<-function(gene)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
#               Simple poltting wrap-up. Uses mostly global variables
{
  current.gene=all.reads[[gene]]
  
  temp<-gene.normalization(current.gene, all.counts[[gene]], genelist[gene])
  gene.gnormalized<-temp[[2]] 
  
  if (length(gene.gnormalized)>0)
  { 
      # Setting up for intensity analysis
      domain=c(min(unlist(gene.gnormalized)), max(unlist(gene.gnormalized))) 
      samples=get.intensity(gene.gnormalized, domain); 
      # Adding SRVF curves and SRVF alignment versions of the curves
      samples.xy<-add.xy(samples, domain); sample1<-samples.xy[1:3];   sample2<-samples.xy[4:6];
      
  fsize=3; PCH=1; CEX=2;  LWD=3;
  png(paste(local_wd,'/Figures/gene',genelist[gene],'.png',sep=""), width=1920, height=1080)
  par(mar=c(6,7,3,5), mgp = c(5, 2, 0));
  
  yrange=c(min(samples), max(samples))
  # Intensity functions
  plot(samples.xy[[1]], col=1, type="l",lty=1, lwd=LWD , main="", pch=PCH, cex=CEX,cex.lab=fsize, cex.axis=fsize, cex.main=fsize, cex.sub=fsize, ylab="density" ,xlab="" , ylim=yrange) 
  lapply(sample1 ,function(dat) lines(dat, col=1 ,lty=1, lwd=LWD));
  lapply(sample2 ,function(dat) lines(dat, col=2 ,lty=1, lwd=LWD)); 
  
  dev.off()
  }
  return(NULL)
}


get.exon.list<-function(Starts, Ends, midpoints)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
{
  exon.list=list()
  n.ex=length(Starts)
  for (n in 1:n.ex)
  {
    exon.list[[n]]=midpoints[which(midpoints<Ends[n] & midpoints >Starts[n])]
  }
  return(exon.list)
}

  normalize01<-function(single.sample, domain)
    #  ARGUMENTS: ---------------------------------------------------------------
  #  VALUES: ------------------------------------------------------------------
  #  DESCRIPTION: --------------------------------------------------------------
  {single.sample<-(single.sample-domain[1])/(domain[2]-domain[1])}


get.intensity<-function(sample, support, multiplier=1) 
  #  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------

{
  n=length(sample);
  #N=max(min(8*median(unlist(sapply(sample, function(s) length(s)))), 256)*multiplier, 20*multiplier)
  #ORRRR:
  N=128*multiplier;
  intensity=array(1/diff(support), c(N,n));  
  bwi=array(0,c(1,n))
  for (i in 1:n)
  {   if (length(sample[[i]])>=read.cutoff)  bwi[i] = bw.nrd0(sample[[i]]);  }
  bwm=mean(bwi)
  
  for (i in 1:n)
  {
    
#Manual BW:
    #if (length(sample[[i]])>=read.cutoff) intensity[[i]]= density(sample[[i]],bw=sqrt(abs(diff(support)))/(log(log(length(sample[[i]])+2*exp(1))*multiplier+exp(1))), from=support[1], to=support[2],n=N);
#Automatic BW:    
    if (length(sample[[i]])>=read.cutoff) intensity[,i]= density(sample[[i]],bw=bwm/2, from=support[1], to=support[2],n=N)$y;
    #*length(sample[[i]])/width
  }
  return(intensity)
}


get.means<-function(f, amp="L0")
  
  #  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
{
  n=dim(f)[2];
  t=seq(0,1, length.out = dim(f)[1])
  
  if(amp=="L0") 
  {
  km<-time_warping(f, time=t, lambda = 0, method = "mean", showplot = F,
                     smooth_data = FALSE, sparam = 10, parallel = F, cores = 4)
    return(km)
  }
  if(amp=="L2")
  {
    f<-srsf_to_f(f,time=t)
    km<-time_warping(f, time=t, lambda = 0, method = "mean", showplot = F,
                     smooth_data = FALSE, sparam = 10, parallel = F, cores = 4)
    return(km)
  }
  
  
  
  
}
#   
get.var<-function(mean1, intensity1, dom)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
{
  d=length(mean1)/abs(diff(dom));
  
  interpolated=array(0,c(length(mean1), length(intensity1)))
  for (i in 1: length(intensity1))
  {
    if (length(intensity1[[i]]$x > read.cutoff)) #This should not be needed if I have forced constant intensoties in get.intensity()
      interpolated[,i]=approx(intensity1[[i]]$x, intensity1[[i]]$y + 0.0001, xout=seq(dom[1],dom[2],length.out = length(mean1)), method = "linear")$y 
  }
  return(sum(unlist(apply(interpolated, 2, function(intens) (sum(   (mean1-intens)^2/d) )  )) ))
}

get.means.L2<-function(intensity1, dom)
#  ARGUMENTS: ---------------------------------------------------------------
#  VALUES: ------------------------------------------------------------------
#  DESCRIPTION: --------------------------------------------------------------
{
  l.out=512;
  interpolated=array(0,c(l.out, length(intensity1)))
  for (i in 1: length(intensity1))
  {
    if (length(intensity1[[i]]$x > read.cutoff)) #This should not be needed if I have forced constant intensoties in get.intensity()
      interpolated[,i]=approx(intensity1[[i]]$x, intensity1[[i]]$y + 0.0001, xout=seq(dom[1],dom[2],length.out = l.out), method = "linear")$y 
  }
  means=apply( interpolated, 1, function(x) mean(x)) #We might need to add normalization here - just in case.
  return(means)}



