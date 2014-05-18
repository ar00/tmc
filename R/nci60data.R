
# date: 26/03/2014
# aim: save the initial data from NCI60 and 
# basic Ensembl annotation in to an R session
# output: ../data/data.rda


# genomic annotations from Ensembl Biomart
biomart = read.table("../data/Ensembl.BioMart.v75.human.txt",header=T,as.is=T,sep="\t")
# NCI60 rna processed data
rna = read.delim("../data/Table_S8_NCI60_transcriptome.txt",
 header=T,comment.char="#",as.is=T,sep="\t")
rna[,3] = gsub(" \\(includes.*\\)","",rna[,3])
# NCI60 protein processed data
prot = read.delim("../data/Table_S3_NCI60_proteome.txt",header=T,as.is=T,sep="\t")
save(biomart,rna,prot,file="../data/data.rda")

