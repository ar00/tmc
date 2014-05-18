
# date: 19/03/2014
# aim: 
# - get IPI-HGNC relations directly from NCI60 files 
#   with RNA and protein processed data
# - get genomic ccordinates for genes profiled at 
#   either level in the NCI60 panel
# output: 
# - ../data/annotation.rda


# load Ensembl BioMart genomic annotations and 
# NCI60 processed data at either RNA or protein level
load("../data/data.rda")

# get mapping between gene symbols and UCSC transcripts,
# along with UCSC transcripts UTR coordinates (chrom,start,end)
kgXgeneSymbol = read.table("../data/kgXref.geneSymbol.max.53utr.txt",
 header=TRUE,as.is=TRUE,sep="\t")

# get IPI-HGNC relations for genes profiled at the protein level 
symbol.start = as.numeric(sapply(prot$Fasta.headers, function(x)
 (regexpr("Gene_Symbol",x)[1])))
ss = substr(prot$Fasta.headers[symbol.start != -1], 
 start=symbol.start[symbol.start != -1], stop=1000000L)
tmp = gsub("Gene_Symbol=","",as.vector(sapply(ss, function(x)
 (strsplit(x," ")[[1]][1]))))
symbol = as.vector(sapply(tmp, function(x)(strsplit(x,";")[[1]][1])))
zrs = rep(0,length(symbol))
ipiHgncProt = data.frame(IPI=zrs,HGNC.symbol=zrs)
ipiHgncProt[,1] = prot$Accession[symbol.start != -1]
ipiHgncProt[,2] = symbol
ipiHgncProt = ipiHgncProt[is.na(ipiHgncProt[,1])==FALSE & 
 is.na(ipiHgncProt[,2])==FALSE,]

# get IPI-HGNC relations for genes profiled at the mRNA level
ipiHgncRna = unique(rna[,2:3])
ipiHgncRna = ipiHgncRna[is.na(ipiHgncRna[,1])==FALSE & 
 is.na(ipiHgncRna[,2])==FALSE,]
colnames(ipiHgncRna) = c("IPI","HGNC.symbol")
# get IPI-HGNC relations for genes profiled at either level
ipiHgnc = unique(rbind(ipiHgncProt,ipiHgncRna))


# get genomic coordinates for NCI60 RNA profiled genes
l = dim(ipiHgnc)[1]
annotation = data.frame(ID = character(1),IPI = character(l),
 HGNC.symbol = character(l),Ensembl.Gene.ID = character(l),
 Chromosome.Name =  character(l),Gene.Start..bp. = integer(l), 
 Gene.End..bp. = integer(l), Strand = character(l),
 UCSC.max5utr = character(l),UCSC.max5utr.start = integer(l),
 UCSC.max5utr.end = integer(l),UCSC.max3utr = character(l),
 UCSC.max3utr.start = integer(l),UCSC.max3utr.end = integer(l))
annotation$ID = paste("ID",1:l,sep="")
annotation$IPI = ipiHgnc[,1] 
annotation$HGNC.symbol = ipiHgnc[,2]
annotation$Ensembl.Gene.ID = biomart$Ensembl.Gene.ID[match(ipiHgnc[,2],
 biomart$HGNC.symbol)]
annotation$Chromosome.Name = biomart$Chromosome.Name[match(ipiHgnc[,2],
 biomart$HGNC.symbol)]
annotation$Gene.Start..bp. = biomart$Gene.Start..bp.[match(ipiHgnc[,2],
 biomart$HGNC.symbol)]
annotation$Gene.End..bp. = biomart$Gene.End..bp.[match(ipiHgnc[,2],
 biomart$HGNC.symbol)]
annotation$Strand = biomart$Strand[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$UCSC.max5utr = kgXgeneSymbol$ucsc.5utr[match(ipiHgnc[,2],
 rownames(kgXgeneSymbol))]
annotation$UCSC.max5utr.start = kgXgeneSymbol$ucsc.5utr.start[match(ipiHgnc[,2], rownames(kgXgeneSymbol))]
annotation$UCSC.max5utr.end = kgXgeneSymbol$ucsc.5utr.end[match(ipiHgnc[,2],
 rownames(kgXgeneSymbol))]
annotation$UCSC.max3utr = kgXgeneSymbol$ucsc.3utr[match(ipiHgnc[,2],
 rownames(kgXgeneSymbol))]
annotation$UCSC.max3utr.start = kgXgeneSymbol$ucsc.3utr.start[match(ipiHgnc[,2],rownames(kgXgeneSymbol))]
annotation$UCSC.max3utr.end = kgXgeneSymbol$ucsc.3utr.end[match(ipiHgnc[,2],
 rownames(kgXgeneSymbol))]
# remove genes with incomplete genomic annotation
ind = apply(annotation[,4:8], 1, function(x) all(is.na(x)))
annotation = annotation[ !ind, ]


save(annotation,file="../data/annotation.rda")









