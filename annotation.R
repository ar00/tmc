
### date: 19/03/2014
### aim: 
### - retrieve IPI-HGNC relations directly from NCI60 RNA processed data
### - retrieve genomic ccordinates
### output: 
### - ./annotation/ipi.hgnc.nci60.txt
### - annotation.rda


# load Ensembl BioMart genomic annotations and NCI60 processed data
load("data.rda")

# retrieve genomic annotation for NCI60 RNA profiled genes
ipiHgnc = unique(rna[,2:3]) 
l = dim(ipiHgnc)[1]
annotation = data.frame(ID = character(1), IPI = character(l), HGNC.symbol = character(l), Ensembl.Gene.ID = character(l), 
 Chromosome.Name =  character(l), Gene.Start..bp. = integer(l), Gene.End..bp. = integer(l), Strand = character(l))
annotation$ID = paste("ID",1:l,sep="")
annotation$IPI = ipiHgnc[,1] 
annotation$HGNC.symbol = ipiHgnc[,2]
annotation$Ensembl.Gene.ID = biomart$Ensembl.Gene.ID[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Chromosome.Name = biomart$Chromosome.Name[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Gene.Start..bp. = biomart$Gene.Start..bp.[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Gene.End..bp. = biomart$Gene.End..bp.[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Strand = biomart$Strand[match(ipiHgnc[,2],biomart$HGNC.symbol)]
# remove genes with incomplete genomic annotation
ind = apply(annotation[,4:8], 1, function(x) all(is.na(x)))
annotation = annotation[ !ind, ]

write.table(annotation,file="./annotation/nci60annotation.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
save(annotation,file="annotation.rda")









