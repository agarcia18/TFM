local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.rstudio.com/"
  r["BioCsoft"] <- "https://bioconductor.org/packages/3.12/bioc"
  r["BioCann"] <- "https://bioconductor.org/packages/3.12/data/annotation"
  r["BioCexp"] <- "https://bioconductor.org/packages/3.12/data/experiment"
  options(repos = r)
})

