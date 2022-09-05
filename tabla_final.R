run <- 1055
a <- read.csv("lineage_report.csv", header=TRUE)
b <- gsub("Consensus_","",a$taxon)
a$strain <- gsub("_threshold_0.6_quality_25","",b)
d <- read.csv("profundidad_ns.txt", header=FALSE, sep=" ")
names(d) <- c("strain","Mean_Coverage","ch")
d$ch <- d$ch*100
e <- merge(a,d,by="strain", all.x=TRUE)
f <- e[,c(1,3,6,18,19)]
names(f) <- c("taxon","lineage","scorpio_call","MeanCoverage","NonN")
g <- data.frame(f$taxon,corrida=rep(run,length(f$taxon)),f$lineage,f$scorpio_call,f$MeanCoverage,f$NonN)
names(g) <- c("taxon","corrida","lineage","scorpio_call","MeanCoverage","NonN")
write.csv(g,"nexstrain2gisaid_vero.csv",row.names=FALSE)
g$taxon <- gsub("_.*","",g$taxon)
write.csv(g,"nexstrain2gisaid.csv", row.names=FALSE)
q("no") ;
