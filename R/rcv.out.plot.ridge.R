
# date: 14/05/2014
# aim: it is based on simple linear regression (no penalty parameter tuneling)
# - plot histogram of RMSPE ratios, M0 RMSPE, M1 RMSPE
# - plot histogram of correlation between preditions and 
#   observations achieved by M1 in testing subset
# - plot histogram of correlation between preditions and
#   observations achieved by M0 in testing subset
# - plot histogram of correlation ratios
#
# return: 
# - ../results/rcv.out.ridge.penalty.pdf


library("GenomicRanges")
library("gplots")


load("../../results/summarizedExperiment.rda")
load("../../results/data.rda")
load("../../results/annotation.rda")


# proteomics data
exprm = as.matrix(assays(sset)[[2]])
# remove matrix rows for IDs artificially added
# to costruct the summarizedExperiment object
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]

# replace 0 values with N/A according to the specification by Prof. Kuster
for(i in 1 : dim(exprm)[1]){
 for(j in 1 : dim(exprm)[2]){
  if(exprm[i,j]==0){ exprm[i,j]=NA }
 }
}

# remove IDs whose protein levels are N/A overall cell lines
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]

# count non-N/A values per gene overall cell lines
nonNA=apply(exprm,1,function(x) length(x[is.na(x)==FALSE]))


# summary of extra-samples prediction error performances 
cvf = read.table("../rcv.summary.ridge.penalty.txt",
 header=T,sep="\t",as.is=TRUE)


pdf("../rcv.out.ridge.penalty.pdf")
par(mfrow=c(2,3))
# plot histogram of average RMSPE ratios, M0 RMSPE, M1 RMSPE
tmp = cbind(cvf$HGNC.symbol,cvf$rmspe.ratio.avg)
tmp = unique(tmp)
rmspe.ratio.avg=as.numeric(tmp[,2])
hist(rmspe.ratio.avg,breaks=50,prob=T,main="",xlab="Average M0/M1 RMSPE ratio")

tmp = cbind(cvf$HGNC.symbol,cvf$m0.rmspe.avg)
tmp = unique(tmp)
m0.rmspe.avg=as.numeric(tmp[,2])
hist(m0.rmspe.avg,breaks=50,prob=T,main="",xlab="Average M0 RMSPE")

tmp = cbind(cvf$HGNC.symbol,cvf$m1.rmspe.avg)
tmp = unique(tmp)
m1.rmspe.avg=as.numeric(tmp[,2])
hist(m1.rmspe.avg,breaks=50,prob=T,main="",xlab="Average M1 RMSPE")

# histogram of CV for parameters tuneling
tmp = cbind(cvf$HGNC.symbol,cvf$cvm)
tmp = unique(tmp)
cvm=as.numeric(tmp[,2])
hist(cvm,breaks=50,prob=T,main="Parameters tuneling",
 xlab="Mean cross-validated error")

tmp = cbind(cvf$HGNC.symbol,cvf$cvsd)
tmp = unique(tmp)
cvsd=as.numeric(tmp[,2])
hist(cvsd,breaks=50,prob=T,main="Parameters tuneling",
 xlab="Standard error of cvm")

tmp = cbind(cvf$HGNC.symbol,cvf$m1.coef.n)
tmp = unique(tmp)
m1.coef.n=as.numeric(tmp[,2])
hist(m1.coef.n,breaks=50,prob=T,main="Parameters tuneling",
 xlab="N. of non-zero coefficients")

tmp = cbind(cvf$HGNC.symbol,cvf$lambda.min)
tmp = unique(tmp)
lambda.min=as.numeric(tmp[,2])
hist(lambda.min,breaks=50,prob=T,main="Parameters tuneling",
 xlab="Lambda giving min cvm")
dev.off()










