library(GEOquery)
library(data.table)

# Script to get data for JAK2 knockdown experiment


# Get data
gds <- getGEO("GSE54645")
show(gds)
data.expr<-data.table(featureNames(gds), exprs(gds[[1]]))
probeset<-getGEO("GPL4685")

# Extract required information and merge
ProbeAndGeneSym<-data.table(Table(probeset)[,c("ID", "Gene Symbol")])
setkey(ProbeAndGeneSym, ID)
setkey(data.expr, V1)
data.expr.anno<-data.expr[ProbeAndGeneSym,nomatch=NA]

# Get names of samples
pData(gds[[1]])[1]


# Average the values corresponding to replicates
data.expr.anno.ave<-cbind(data.expr.anno[,c(1,14)],
  apply(data.expr.anno[,2:4], 1, mean),
  apply(data.expr.anno[,5:7], 1, mean),
  apply(data.expr.anno[,8:10], 1, mean),
  apply(data.expr.anno[,11:13], 1, mean)
)

# Change the name
names(data.expr.anno.ave)<-c("Probename", "Gene Symbol", "sh-GFP","sh-JAK2-shRNA1", "sh-JAK2-shRNA2", "shRNA-control")

write.table(data.expr.anno.ave, "data_jak2_knockdown.txt", sep="\t", row.names=F)
