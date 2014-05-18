
# date: 16/04/2014
# aim: report RBP-target interactions as 1 (>=1 binding site in 
#  either 5' or 3'UTR of the gene) or 0 (otherwise). 
# return: 
# - ../data/rbp.target.binary.txt


# report the RBP-target interactions 
# as 1 (>=1 RBP binding site in either 5' or 3'UTR) or 0 
f=read.delim("../data/rbp.match.qv.20.txt",
 header=TRUE,as.is=TRUE,sep="\t")
zrs=rep(NA,dim(f)[1])
int=data.frame(rbp=zrs, target=zrs)
int[,1]=f$rbp
int[,2]=f$target.hgnc
int=unique(int)
rbp=unique(int$rbp)
target=unique(int$target)
m=matrix(0,nrow=length(target),ncol=length(rbp),
 dimnames=list(target,rbp))
for(i in 1 : dim(int)[1]){
 m[match(int$target[i],rownames(m)),
  match(int$rbp[i],colnames(m))] = 1
 show(i)
}

write.table(m,"../data/rbp.target.binary.txt",row.names=T,col.names=T,quote=F,sep="\t")











