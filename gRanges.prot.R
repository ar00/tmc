
# date: 26/03/2014
# aim: define genomic ranges for NCI60 processed protein data.
#  Metadata are protein expression levels according to lbel-free quantification (LFQ) in NCI60 cell lines.
# output: gRanges.protein.rda

library("GenomicRanges")

load("data.rda")
load("annotation.rda")


# obtain common genes to annotation file and protein processed data file
cg = intersect(prot$Accession,annotation$IPI)
# filter out uncommon genes
prot = prot[is.na(match(prot$Accession,cg))==FALSE,]
annotation = annotation[is.na(match(annotation$IPI,cg))==FALSE,]
# filter out annotations of duplicated IPIs
annotation = annotation[duplicated(annotation$IPI)==FALSE,]

# sort objects by IPI alphabetic order
annotation = annotation[order(annotation$IPI),]
expr = prot[order(prot$Accession),]

expr = expr[,grepl("LFQ",colnames(expr))]
rownames(expr) = prot$Accession
colnames(expr) = gsub("LFQ_","",colnames(expr))

# define genomic ranges
gr.prot = GRanges(seqnames = annotation$Chromosome.Name,ranges = IRanges(start=annotation$Gene.Start..bp.,
 end = annotation$Gene.End..bp., names = annotation$ID),strand = annotation$Strand,
  IPI = annotation$IPI, expr)

save(gr.prot,file="gRanges.prot.rda")


