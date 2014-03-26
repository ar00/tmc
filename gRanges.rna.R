
# date: 26/03/2014
# aim: define genomic ranges for NCI60 processed RNA data. 
#  Metadata are summarized gene expression values in log10 scale. 
# output: gRanges.rna.rda



library(GenomicRanges)

load("data.rda")
load("annotation.rda")

# obtain common genes to annotation file and RNA processed data file
cg=intersect(rna[,3],annotation[,3])
# filter out uncommon genes  
rna = rna[is.na(match(rna[,3],cg))==FALSE,]
annotation = annotation[is.na(match(annotation$HGNC.symbol,cg))==FALSE,]
# filter out rows corresponding to unmapped probes
rna = rna[is.na(rna[,3])==FALSE,] 
# filter out duplicated annotations of gene with multiple IPIs
annotation = annotation[duplicated(annotation$HGNC.symbol)==FALSE,]


# summarize probe-level intensity values (log10) to gene-level intensity values (average)
g = unique(rna[,3])
expr = data.frame()
for(i in 1 : length(g)){
 tmp = rna[is.na(match(rna[,3],g[i]))==FALSE,9:67] # columns 9-67 correspond to NCI60 cell lines
 if(is.null(dim(tmp)[1])){ expr = rbind(expr, tmp) }
 if(!is.null(dim(tmp)[1])){ expr = rbind(expr,colMeans(tmp)) }
 show(i)
}
colnames(expr) = colnames(rna)[9:67]
rownames(expr) = g


# sort objects by HGNC symbols alphabetic order
annotation = annotation[order(annotation$HGNC.symbol),]
expr = expr[order(rownames(expr)),]


# define genomic ranges
gr.rna = GRanges(seqnames = annotation$Chromosome.Name,ranges = IRanges(start=annotation$Gene.Start..bp.,
 end = annotation$Gene.End..bp., names = annotation$ID),strand = annotation$Strand,
  HGNC.symbol = annotation$HGNC.symbol, expr)

save(gr.rna,file="gRanges.rna.rda")




