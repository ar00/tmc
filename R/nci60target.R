
# date: 15/05/2014
# aim:
# - distribution of average target protein log(Intensity,10) on available lines
# - distribution of average target mRNA log(Intensity,10) on available lines
# - distribution of SD target protein log(Intensity,10) on available lines
# - distribution of SD target mRNA log(Intensity,10) on available lines
# return:
# - ./nci60target.pdf


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
# obtain RBP target protein data
target.ID=annotation$ID[is.na(match(annotation$HGNC.symbol,
 rownames(rbp.db)))==FALSE]
target.prot=exprm[is.na(match(rownames(exprm),target.ID))==FALSE,]

# obtain RNA data
exprm = as.matrix(assays(sset)[[1]])
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]
rna.exprm=exprm
# obtain RBP RNA data
rbp.rna=exprm[is.na(match(rownames(exprm),
 annotation$ID[is.na(match(annotation$HGNC.symbol,rbp.pwm))==FALSE]))==FALSE,]
# obtain RBP target protein data
target.ID=annotation$ID[is.na(match(annotation$HGNC.symbol,
 rownames(rbp.db)))==FALSE]
target.rna=exprm[is.na(match(rownames(exprm),target.ID))==FALSE,]



# subset the RBP-target interaction matrix for RBPs expressed in NCI-60
rbp.ID=annotation$ID[is.na(match(annotation$HGNC.symbol,colnames(rbp.db)))==FALSE]
rbp.ID.prot = rownames(prot.exprm)[is.na(match(rownames(prot.exprm),
 rbp.ID))==FALSE]
rbp.HGNC.prot = annotation$HGNC.symbol[is.na(match(annotation$ID,rbp.ID.prot))==FALSE]
rbp.db = rbp.db[,match(rbp.HGNC.prot,colnames(rbp.db))]


pdf("nci60target.pdf")
par(mfrow=c(2,2))

# mean and SD of target expression data in NCI-60 protein profiles 
targetMeans=as.numeric(apply(target.prot,1,function(x)(mean(x[is.na(x)==FALSE]))))
targetSDs=as.numeric(apply(target.prot,1,
 function(x)(sd(x[is.na(x)==FALSE]))))
# mean and SD of protein expression data in NCI-60 protein profiles
protMeans=as.numeric(apply(prot.exprm,1,function(x)(mean(x[is.na(x)==FALSE]))))
protSDs=as.numeric(apply(prot.exprm,1,
 function(x)(sd(x[is.na(x)==FALSE]))))

hist(protMeans,xlab="Mean protein (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,.5+max(max(hist(protMeans,plot=F)$density),max(hist(targetMeans,plot=F)$density))))
hist(targetMeans,xlab="Mean protein (log(Intensity,10))",
 ylim=c(0,.5+max(max(hist(protMeans,plot=F)$density),
 max(hist(targetMeans,plot=F)$density))),prob=T,breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected protein","Detected protein target"))

hist(protSDs,xlab="SD protein (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,2+max(max(hist(protSDs,plot=F)$density),max(hist(targetSDs,plot=F)$density))))
hist(targetSDs,xlab="SD protein (log(Intensity,10))",prob=T,
 ylim=c(0,2+max(max(hist(protSDs,plot=F)$density),
 max(hist(targetSDs,plot=F)$density))),breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected protein","Detected protein target"))



# mean and SD of target expression data in NCI-60 protein profiles
targetMeans=as.numeric(apply(target.rna,1,function(x)(mean(x[is.na(x)==FALSE]))))
targetSDs=as.numeric(apply(target.rna,1,
 function(x)(sd(x[is.na(x)==FALSE]))))
# mean and SD of protein expression data in NCI-60 protein profiles
rnaMeans=as.numeric(apply(rna.exprm,1,function(x)(mean(x[is.na(x)==FALSE]))))
rnaSDs=as.numeric(apply(rna.exprm,1,
 function(x)(sd(x[is.na(x)==FALSE]))))

hist(rnaMeans,xlab="Mean mRNA (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,.5+max(max(hist(rnaMeans,plot=F)$density),max(hist(targetMeans,plot=F)$density))))
hist(targetMeans,xlab="Mean mRNA (log(Intensity,10))",
 ylim=c(0,.5+max(max(hist(rnaMeans,plot=F)$density),
 max(hist(targetMeans,plot=F)$density))),prob=T,breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected mRNA","Detected mRNA target"))

hist(rnaSDs,xlab="SD mRNA (log(Intensity,10))",prob=T,breaks=25,main="",
 ylim=c(0,2+max(max(hist(rnaSDs,plot=F)$density),max(hist(targetSDs,plot=F)$density))))
hist(targetSDs,xlab="SD mRNA (log(Intensity,10))",prob=T,
 ylim=c(0,2+max(max(hist(rnaSDs,plot=F)$density),
 max(hist(targetSDs,plot=F)$density))),breaks=25,add=T,col=rgb(1,0,0,.25))
legend("topright",fill=c(rgb(1,1,1),rgb(1,0,0,.25)),bty="n",
 legend=c("Detected mRNA","Detected mRNA target"))

dev.off()


















