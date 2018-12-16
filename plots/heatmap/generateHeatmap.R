setwd("~/Documents/AM205/project/AM205-Project/")

library("RColorBrewer")
predictions<-read.csv("data/PD-L1_predicted_direction.csv", sep="\t")
predictions<-predictions[,-16]

pdf("heatmap.pdf")
heatmap(t(as.matrix(predictions[,2:40])), Rowv = F, Colv=NA , col=brewer.pal(3,"Spectral"), 
        scale="none", labCol = predictions[,1])

legend("topleft", c("Unchanged","Increase", "Decrease"), fill = brewer.pal(3,"Spectral"))
dev.off()
s