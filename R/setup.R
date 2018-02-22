# chiptune/R/setup.R: setup script for qugacomp R scripts
# Tested on bioconductor/release_base2:R3.4.3_Bioc3.6

#
# Get path to the directory of this script
#
initial.options <- commandArgs(trailingOnly=FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#
# Package install and load
#

# QuGAcomp
if (!require("QuGAcomp")) {
  update.packages(checkBuilt=TRUE, ask=FALSE, repos="https://cran.ism.ac.jp/")
  # GenomicRanges from Bioconductor
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  # Install curl package
  install.packages("curl", repos="https://cran.ism.ac.jp/")
  library("curl")
  # QuGAcomp from source file
  qugacomp.package.url <- "https://github.com/dritoshi/QuGAcomp/raw/master/QuGAcomp_0.99.2.tar.gz"
  qugacomp.package.loc <- file.path(".", basename(qugacomp.package.url))
  curl_download(qugacomp.package.url, qugacomp.package.loc)
  install.packages(qugacomp.package.loc)
  file.remove(qugacomp.package.loc)
}
library("QuGAcomp")

# pforeach
if (!require("pforeach")) {
  install.packages("devtools", repos="https://cran.ism.ac.jp/")
  devtools::install_github("hoxo-m/pforeach")
}
library("pforeach")

# Corrplot
if (!require("corrplot")) {
  install.packages("corrplot", repos="https://cran.ism.ac.jp/")
}
library("corrplot")

sessionInfo()
