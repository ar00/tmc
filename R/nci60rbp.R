
# date: 15/05/2014
# aim:
# - plot n. targets per RBP and vice versa (not target and RBP are detected 
#  in NCI-60 protein profiles)
# - distribution of average RBP protein log(Intensity,10) on available lines 
# - distribution of average RBP mRNA log(Intensity,10) on available lines
# - distribution of SD RBP protein log(Intensity,10) on available lines
# - distribution of SD RBP mRNA log(Intensity,10) on available lines
# return: 
# - ./nci60rbp.pdf


library(GenomicRanges)


load("../../results/data.rda")
load("../../results/summarizedExperiment.rda")
load("../../annotation/CISBP-RNA.rda")
load("../../results/annotation.rda")
load("../../results/annotationHG19.rda")


# obtain RBP->target interactions (binary value 1/0)
rbp.db=read.delim("../../results/CISBP-RNA/rbp.target.binary.txt",
 header=TRUE,as.is=TRUE,sep="\t")


# obtain NCI-60 protein data
exprm = as.matrix(assays(sset)[[2]])
ind = apply(exprm, 1, function(x) all(is.na(x)))
# remove N/A rows inserted in the construction of summarizedExperiment object
exprm = exprm[!ind,]
# replace 0 values with N/A according to the specification by Prof. Kuster
for(i in 1 : dim(exprm)[1]){
 for(j in 1 : dim(exprm)[2]){
  if(exprm[i,j]==0){ exprm[i,j]=NA }
 }
}

# remove N/A rows from NCI-60 protein data
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]
# rename protein data
prot.exprm=exprm
# obtain RBP protein data
rbp.prot=exprm[is.na(match(rownames(exprm),
 annotation$ID[is.na(match(annotation$HGNC.symbol,rbp.pwm))==FALSE]))==FALSE,]


# obtain RNA data
exprm = as.matrix(assays(sset)[[1]])
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]
rna.exprm=exprm
# obtain RBP RNA data
rbp.rna=exprm[is.na(match(rownames(exprm),
 annotation$ID[is.na(match(annotation$HGNC.symbol,rbp.pwm))==FALSE]))==FALSE,]



# subset the RBP-target interaction matrix for RBPs expressed in NCI-60
rbp.ID=annotation$ID[is.na(match(annotation$HGNC.symbol,colnames(rbp.db)))==FALSE]
rbp.ID.prot = rownames(prot.exprm)[is.na(match(rownames(prot.exprm),
 rbp.ID))==FALSE]
rbp.HGNC.prot = annotation$HGNC.symbol[is.na(match(annotation$ID,rbp.ID.prot))==FALSE]
rbp.db = rbp.db[,match(rbp.HGNC.prot,colnames(rbp.db))]


pdf("nci60rbp.pdf",height=8)
par(mfrow=c(4,2))

# histogram of number of inferred targets per RBP
hist(as.numeric(colSums(rbp.db)),breaks=25,prob=T,
 main="Detected protein",cex.main=.75,xlab="N. targets per RBP")
# histogram of number of inferred RBPs per target
hist(as.numeric(rowSums(rbp.db)),breaks=25,prob=T,
 main="Detected protein",cex.main=.75,xlab="N. RBPs per target")

# number of N/A values per RBP in NCI-60 protein profiles
rbpNAprot = as.numeric(apply(rbp.prot,1,function(x)(length(x[is.na(x)==T]))))
# number of N/A values per RBP in NCI-60 mRNA profiles
rbpNArna = as.numeric(apply(rbp.rna,1,function(x)(length(x[is.na(x)==T]))))
hist(rbpNAprot,prob=T,main="",xlab="N/A frequency in RBP protein data")
barplot(rbpNArna,ylim=c(0,1),xlab="RBP",ylab="mRNA N/A frequency")




# mean and SD of RBP expression data in NCI-60 protein profiles 
rbpMeans=as.numeric(apply(rbp.prot,1,function(x)(mean(x[is.na(x)==FALSE]))))
rbpSDs=as.numeric(apply(rbp.prot,1,
 function(x)(sd(x[is.na(x)==FALSE]))))
# mean and SD of protein expression data in NCI-60 protein profiles
protMeans=as.numeric(apply(prot.exprm,1,function(x)(mean(x[is.na(x)==FALSE]))))
protSDs=as.numeric(apply(prot.exprm,1,
 function(x)(sd(x[is.na(x)==FALSE]))))

# KS-test between RBP vs non-RBP protein expression values
x.rbp=as.numeric(apply(rbp.prot,1,function(x)(mean(x[is.na(x)==FALSE]))))
x.non.rbp=as.numeric(apply(prot.exprm[is.na(match(rownames(prot.exprm),
 rownames(rbp.prot)))==TRUE,],1,function(x)(mean(x[is.na(x)==FALSE]))))
ks.test(x.rbp,x.non.rbp,alternative="less")

hist(protMeans,xlab="Mean protein (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,.5+max(max(hist(protMeans,plot=F)$density),max(hist(rbpMeans,plot=F)$density))))
hist(rbpMeans,xlab="Mean protein (log(Intensity,10))",
 ylim=c(0,.5+max(max(hist(protMeans,plot=F)$density),
 max(hist(rbpMeans,plot=F)$density))),prob=T,breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected protein","Detected RBP"))

hist(protSDs,xlab="SD protein (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,.5+max(max(hist(protSDs,plot=F)$density),max(hist(rbpSDs,plot=F)$density))))
hist(rbpSDs,xlab="SD protein (log(Intensity,10))",prob=T,
 ylim=c(0,.5+max(max(hist(protSDs,plot=F)$density),
 max(hist(rbpSDs,plot=F)$density))),breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected protein","Detected RBP"))



# mean and SD of RBP expression data in NCI-60 protein profiles
rbpMeans=as.numeric(apply(rbp.rna,1,function(x)(mean(x[is.na(x)==FALSE]))))
rbpSDs=as.numeric(apply(rbp.rna,1,
 function(x)(sd(x[is.na(x)==FALSE]))))
# mean and SD of protein expression data in NCI-60 protein profiles
rnaMeans=as.numeric(apply(rna.exprm,1,function(x)(mean(x[is.na(x)==FALSE]))))
rnaSDs=as.numeric(apply(rna.exprm,1,
 function(x)(sd(x[is.na(x)==FALSE]))))

# KS-test between RBP vs non-RBP mRNA expression values
x.rbp=as.numeric(apply(rbp.rna,1,function(x)(mean(x[is.na(x)==FALSE]))))
x.non.rbp=as.numeric(apply(rna.exprm[is.na(match(rownames(rna.exprm),
 rownames(rbp.rna)))==TRUE,],1,function(x)(mean(x[is.na(x)==FALSE]))))
ks.test(x.rbp,x.non.rbp,alternative="less")


hist(rnaMeans,xlab="Mean mRNA (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,.5+max(max(hist(rnaMeans,plot=F)$density),max(hist(rnaMeans,plot=F)$density))))
hist(rbpMeans,xlab="Mean mRNA (log(Intensity,10))",
 ylim=c(0,.5+max(max(hist(rbpMeans,plot=F)$density),
 max(hist(rbpMeans,plot=F)$density))),prob=T,breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected mRNA","Detected RBP"))

hist(rnaSDs,xlab="SD mRNA (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,4+max(max(hist(rnaSDs,plot=F)$density),max(hist(rnaSDs,plot=F)$density))))
hist(rbpSDs,xlab="SD mRNA (log(Intensity,10))",prob=T,
 ylim=c(0,4+max(max(hist(rbpSDs,plot=F)$density),
 max(hist(rbpSDs,plot=F)$density))),breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected mRNA","Detected RBP"))

dev.off()


















