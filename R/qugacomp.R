# A script to exec quantitative genome annotation comparison (QuGAcomp) https://github.com/dritoshi/QuGAcomp
# Tested on bioconductor/release_base2:R3.4.3_Bioc3.6

#
# Package update and install
#
update.packages()

# Corrplot
install.packages("corrplot")

# GenomicRanges from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

# QuGAcomp from source file
qugacomp.package.url <- "https://github.com/dritoshi/QuGAcomp/raw/master/QuGAcomp_0.99.2.tar.gz"
qugacomp.package.loc <- file.path(".", basename(qugacomp.package.url))
download.file(qugacomp.url, qugacomp.package.loc, method="curl")
install.packages(qugacomp.package.loc)

#
# Load libraries
#
packages.required <- c("corrplot", "QuGAcomp")
lapply(packages.required, require, character.only=TRUE)

#
# Download files
#

# TODO: download genome info file from github, and run script to get bed files

#
# Load files
#
genome.length.file <- "./hg19.info"
experiments <- c("DRX013180", "SRX003937", "SRX022562")

#
# QuGAcomp
#
for (exp in experiments) {
  bed.file <- file.path(".", paste(exp, ".bed", sep=""))
  gr <- loadBedFile(bed.file, genome.length.file)
  fat <- fat(gr, 200)
  unistd <- unifyStrand(fat)
  cov <- coverage(unistd)
  bin500 <- lapply(cov, function(x)rleBinning(x, 500))
  bin500 <- flatRleList(bin500)
  assign(paste(exp, ".bin500", sep=""), bin500)
}
