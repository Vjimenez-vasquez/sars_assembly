a <- read.csv("lineage_report.csv", header=TRUE)
b <- gsub("Consensus_","",a$taxon)
a$strain <- gsub("_threshold_0.6_quality_25","",b)
d <- read.csv("profundidad_ns.txt", header=FALSE, sep=" ")
names(d) <- c("strain","Mean_Coverage","%Non N")
e <- merge(a,d,by="strain", all.x=TRUE)
f <- e[,c(1,3,6,18,19)]
write.csv(f,"nexstrain2gisaid.csv",row.names=FALSE)
q("no") ;
