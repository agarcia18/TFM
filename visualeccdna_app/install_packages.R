using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  n<-length(need)
  if(n>0){
    libsmsg<-if(n>2) paste(paste(need[1:(n-1)],collapse=", "),",",sep="") else need[1]
    print(libsmsg)
    if(n>1){
      install.packages(need)
      lapply(need,require,character.only=TRUE)
    }
  }
}

bioc_using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  n<-length(need)
  if(n>0){
    libsmsg<-if(n>2) paste(paste(need[1:(n-1)],collapse=", "),",",sep="") else need[1]
    print(libsmsg)
    if(n>1){
      libsmsg<-paste(libsmsg," and ", need[n],sep="")
    }
      BiocManager::install(need)
      lapply(need,require,character.only=TRUE)
    }
  }


using("shinydashboard","dplyr","ggplot2","stringr","plotly","DT")
install.packages("BiocManager",version=3.13)
bioc_using("TxDb.Hsapiens.UCSC.hg38.knownGene","GenomicRanges","org.Hs.eg.db","RITANdata","RITAN","Biostrings","BSgenome.Hsapiens.UCSC.hg38")
