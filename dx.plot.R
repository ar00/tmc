# This function makes several diagnostic plots for selectd models.
# Arguments are: 
# arg1: NCI60 model id (e.g. ID107)
# arg2: file of model goodness-of-fit measures (self-regualted models are excluded)  
dx.plot <- function(x,fn,fitf){

 f=read.table(fitf,header=T,as.is=T,sep="\t")
 # M0,M1 fitted values are necessary for diagnostic plots
 models = data.frame(ID=f$ID,cl=f$Cell.line, m0.fitted=f$m0.fitted,
  m1.fitted=f$m1.fitted,response=f$Y)


 targetID=x # ID of interest to dispoaly
 hgnc=annotation$HGNC.symbol[match(targetID,annotation$ID)]

 # RNA levels 
 X0=as.numeric(t(rna.exprm[match(targetID,rownames(rna.exprm)),]))

 # PROT levels
 Y=as.numeric(t(prot.exprm[match(targetID,rownames(prot.exprm)),]))

 # replace NA values with closest non-NA values 
 if(length(Y[is.na(Y)==T])>0){Y=na_impute_for_vector(Y)}
 if(length(X0[is.na(X0)==T])>0){X0=na_impute_for_vector(X0)}

 # standardize variables
 X0 = (X0-mean(X0[is.na(X0)==F]))/sd(X0[is.na(X0)==F])
 Y = (Y-mean(Y))/sd(Y)
 
 # calculate RNA vs PROT correlation
 prot.rna.cor=signif(cor.test(Y,X0,method="spearman")$estimate,digits=2)

 # subset model performances for the gene of interest
 # in order to acquire the fitted values for each cell line 
 m=models[is.na(match(models$ID,targetID))==F,]
 lines = unique(m$cl)

 # calculate the standardized M1 residuals
 fitted.m1.by.line = list(); response.by.line = c();
 residuals.m1.by.line = list(); std.residuals.m1.by.line = list();
 for(l in 1 : length(lines)){
  fitted.m1.by.line[[lines[l]]] = m$m1.fitted[is.na(match(m$cl,lines[l]))==F]
  response.by.line[lines[l]] = m$response[is.na(match(m$cl,lines[l]))==F][1]
  residuals.m1.by.line[[lines[l]]] = fitted.m1.by.line[[lines[l]]]-response.by.line[lines[l]] 
 }
 rmse.m1=sqrt(sum(unlist(residuals.m1.by.line)^2)/dim(m)[1])
 for(l in 1 : length(lines)){
  std.residuals.m1.by.line[[lines[l]]] = residuals.m1.by.line[[lines[l]]] / rmse.m1
 }

 # test normality of M1 residuals
 norm.std.residuals.m1.pv = signif(ad.test(unlist(std.residuals.m1.by.line))$p.value,digits=2)

 # calculate the standardized M0 residuals
 fitted.m0.by.line = list(); response.by.line = c();
 residuals.m0.by.line = list(); std.residuals.m0.by.line = list();
 for(l in 1 : length(lines)){
  fitted.m0.by.line[[lines[l]]] = m$m0.fitted[is.na(match(m$cl,lines[l]))==F]
  response.by.line[lines[l]] = m$response[is.na(match(m$cl,lines[l]))==F][1]
  residuals.m0.by.line[[lines[l]]] = fitted.m0.by.line[[lines[l]]]-response.by.line[lines[l]]
 }
 rmse.m0=sqrt(sum(unlist(residuals.m0.by.line)^2)/dim(m)[1])
 for(l in 1 : length(lines)){
  std.residuals.m0.by.line[[lines[l]]] = residuals.m0.by.line[[lines[l]]] / rmse.m0
 }

 # test normality of M0 residuals
 norm.std.residuals.m0.pv = signif(ad.test(unlist(std.residuals.m0.by.line))$p.value,digits=2)

 # plot M1 fitted values (X-axis) vs standardized residuals of M1 (Y-axis)  
 plot(unlist(fitted.m1.by.line),unlist(std.residuals.m1.by.line),ylim=c(-2.5,2.5),
  main=paste(targetID,hgnc,sep="|"),cex=.5,pch=19,cex.main=.85,
  xlab="fitted.m1",ylab="std.residuals.m1",las=2,col=rgb(1,0,0,.25))
 abline(h=0,lty=2,lwd=1.5)
 mtext(side=3, line=.5, text=paste("RMSPE ratio=",
  fn$rmspe.ratio.avg[match(hgnc,fn$HGNC.symbol)],sep=""),
  cex=.65,adj=.25)
 
 # plot response (X-axis) vs M1 fitted (Y-axis)
 plot(rep(response.by.line[1], length(fitted.m1.by.line[[1]])), fitted.m1.by.line[[1]], 
  cex=.5,pch=19,col=rgb(1,0,0,.25),axes=T, ylab="fitted.m1", xlab="response",
  main=paste(targetID,hgnc,sep="|"),cex.main=.85,las=2,
  xlim=c(min(response.by.line)-.5,max(response.by.line)+.5),
  ylim=c(min(unlist(fitted.m1.by.line)),max(unlist(fitted.m1.by.line))))
 for(l in 2 : length(lines)){
  points(rep(response.by.line[l],length(fitted.m1.by.line[[l]])), fitted.m1.by.line[[l]],
   cex=.5, pch=19,col=rgb(1,0,0,.25))
 }

 # M1 QQ-plot
 qqnorm(unlist(std.residuals.m1.by.line),xlab="normal.scores",ylab="std.residuals.m1",
  cex.main=.75,las=2,col=rgb(1,0,0,.25),main=paste(targetID,hgnc,sep="|"),
  cex=.5,pch=19)
 qqline(unlist(std.residuals.m1.by.line), col = rgb(1,0,0,.25))
 mtext(side=3, line=.65, text=paste("norm.test.pv=",norm.std.residuals.m1.pv,
  sep=""),cex=.5,adj=.25)

 # plot M0 fitted values (X-axis) vs residuals of M0 (Y-axis)
 plot(unlist(fitted.m0.by.line),unlist(std.residuals.m0.by.line),ylim=c(-2.5,2.5),
  main=paste(targetID,hgnc,sep="|"),cex=.5,pch=19,cex.main=.75,
  xlab="fitted.m0",ylab="std.residuals.m0",las=2,col=rgb(1,0,0,.25))
 abline(h=0,lty=2,lwd=1.5)
 mtext(side=3,line=.65,text=paste("RMSPE=",
  fn$m0.rmspe.avg[match(hgnc,fn$HGNC.symbol)],sep=""),cex=.5,adj=.25)

 # plot response (X-axis) vs M0 fitted (Y-axis)
 plot(rep(response.by.line[1], length(fitted.m0.by.line[[1]])), fitted.m0.by.line[[1]],
  cex=.5,pch=19,col=rgb(1,0,0,.25),axes=T, ylab="fitted.m0", xlab="response",
  main=paste(targetID,hgnc,sep="|"),cex.main=.85,las=2,
  xlim=c(min(response.by.line)-.5,max(response.by.line)+.5),
  ylim=c(min(unlist(fitted.m0.by.line)),max(unlist(fitted.m0.by.line))))
 for(l in 2 : length(lines)){
  points(rep(response.by.line[l],length(fitted.m0.by.line[[l]])), fitted.m0.by.line[[l]],
   cex=.5, pch=19,col=rgb(1,0,0,.25))
 }

 # plot RNA vs PROT
 plot(X0[order(X0)],Y[order(X0)],col=rgb(1,0,0,1),cex.axis=.75,cex=.5,cex.main=.75,
  las=2,main=paste(targetID,hgnc,sep="|"),xlab="RNA",ylab="PROT",pch=19)
 mtext(paste("cor(RNA,PROT)",prot.rna.cor,sep="="),
  side=3,line=.5,cex=.65,adj=.25)

}
