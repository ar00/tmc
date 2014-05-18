
# date: 29/03/2014
# aim: R objects containing RBPs per motifs (PWMs) and vice versa 
#  from the CISBP-RNA database (version 0.6)
# output: 
# - ../data/CISBP-RNA.rda


library(hash)

# obtain the mapping of identifiers between Ensembl gene ids and HGNC symbols
ids = read.delim("../data/Ensembl.BioMart.v75.human.txt",
 header=T,as.is=T,sep="\t")

# obtain human cisBP-RNA data
cisbpRNAf=read.table("../data/RBP_Information_all_motifs.txt",
 header=T,sep="\t",as.is=T)

## extract only CISBP-RNA data supported by SELEX and RNAcompete
lgc=array(FALSE, dim(cisbpRNAf)[1])
lgc[is.na(match(cisbpRNAf$Motif_Type, "SELEX"))==FALSE] = TRUE
lgc[is.na(match(cisbpRNAf$Motif_Type, "RNAcompete"))==FALSE] = TRUE
cisbpRNAf=cisbpRNAf[lgc==TRUE,]

## obtain the RBPs provided with a motif in cisBP-RNA
rbp.cisbprna = unique(ids[is.na(match(ids[,1], cisbpRNAf$DBID))==FALSE,3])
rbp.cisbprna = rbp.cisbprna[nchar(rbp.cisbprna)>0]
rbp.cisbprna.ensembl = unique(ids[is.na(match(ids[,3], rbp.cisbprna))==FALSE,1])

## obtain the list of motifs in cisBP-RNA
motif.ids = cisbpRNAf$Motif_ID[is.na(match(cisbpRNAf$DBID, ids[,1]))==FALSE]

## obtain the RBPs corresponding to each motif in cisBP-RNA
rbp.x.pwm = hash()
for(i in 1 : length(motif.ids)){
 tmp.hgnc = unique(ids[match(cisbpRNAf$DBID[is.na(match(cisbpRNAf$Motif_ID, 
  motif.ids[i]))==FALSE], ids[,1]),3])
 tmp.hgnc = tmp.hgnc[is.na(tmp.hgnc)==FALSE]
 rbp.x.pwm[motif.ids[i]] = tmp.hgnc 
}

## obtain the PWMs corresponding to each RBP in cisBP-RNA
pwm.x.rbp = hash()
for(i in 1 : length(rbp.cisbprna.ensembl)){
 tmp = unique(cisbpRNAf$Motif_ID[is.na(match(cisbpRNAf$DBID, 
  rbp.cisbprna.ensembl[i]))==FALSE])
 pwm.x.rbp[ids[match(rbp.cisbprna.ensembl[i], ids[,1]),3]] = tmp
}

rbp.pwm = keys(pwm.x.rbp) 
pwm = keys(rbp.x.pwm)

save(rbp.pwm,pwm,rbp.x.pwm,pwm.x.rbp,file="../data/CISBP-RNA.rda")




















