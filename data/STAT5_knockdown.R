library(GEOquery)
library(data.table)
library(biomaRt)
library(GenomicRanges)
library(Homo.sapiens)
library(magrittr)

# Get data
gds <- getGEO("GSE68993")
show(gds)
data.expr<-data.table(featureNames(gds), exprs(gds[[1]]))
probeset<-getGEO("GPL16686", )
allIDs<-unique(Table(probeset)$GB_ACC)

# Get the probeset data into right format
probeData<-data.frame(matrix(unlist(strsplit(Table(probeset)$SPOT_ID, ":|-")), nrow=53981, byrow=T), Table(probeset)$RANGE_STRAND )
names(probeData)<-c("chr","start","end","strand")
probeData.clean<-probeData[probeData$end != "unknown",]

# Generate the granges object
probeData.clean$start<-as.numeric(as.character(probeData.clean$start))
probeData.clean$end<-as.numeric(as.character(probeData.clean$end))
probeData.clean$chr<-droplevels(probeData.clean$chr)
probeDataGranges<-makeGRangesFromDataFrame(probeData.clean)



hg19<-genes(Homo.sapiens, columns="SYMBOL")
intersect(subjectHits, hg19)
a<-findOverlaps(probeDataGranges, hg19)
probeIndex<-data.frame(a)$queryHits
geneSymbolIndex<-data.frame(a)$subjectHits
finalAnnoTable<-data.frame(probeDataGranges[probeIndex,], as.character(hg19$SYMBOL[geneSymbolIndex]))

spotIDNew<-paste0(finalAnnoTable$seqnames, ":", finalAnnoTable$start, "-", finalAnnoTable$end, sep="")
spotIDNew.df<-data.frame(spotIDNew, as.character(finalAnnoTable[,6]))
names(spotIDNew.df)<-c("SPOT_ID", "GENE_SYMBOL")

spotIDNew.df<-merge(data.frame(Table(probeset)), spotIDNew.df)
spotIDNew.df$ID<-as.character(spotIDNew.df$ID)

# Add Gene symbol to expression
names(data.expr)[1]<-"ID"
data.expr.anno<-merge(data.expr, spotIDNew.df[,c(2,9)])

data.expr.anno.ctl<-apply(data.expr.anno[,2:4], 1, mean)
data.expr.anno.kdSTAT5<-apply(data.expr.anno[,8:10], 1, mean)
data.expr.anno.final<-data.frame(data.expr.anno.ctl, data.expr.anno.kdSTAT5, data.expr.anno$GENE_SYMBOL)

write.table(data.expr.anno.final, "data_stat5_knockdown.txt", sep="\t", row.names=F)

