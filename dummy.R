circ <- as.data.frame(read.table("C:/Users/ainho/eccDNA/BED/SRR6315423_unknown_circle.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))

dummy_circ <- circ[,1:3]

n<-nrow(dummy_circ)

dummy_circ[4:5]<-rep(0,n)
dummy_circ[6]<-rep(100,n)
dummy_circ[7:11] <-rep(0,n)

write.table(dummy_circ, file = "dummy_circ.bed", sep="\t",row.names=FALSE,col.names = FALSE)

