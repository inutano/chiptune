# A script to exec quantitative genome annotation comparison (QuGAcomp) https://github.com/dritoshi/QuGAcomp
# Tested on bioconductor/release_base2:R3.4.3_Bioc3.6

#
# Get path to the directory of this script
#
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#
# Package update, install, and load
#
update.packages(checkBuilt=TRUE, ask=FALSE)

# Corrplot
if (!require("corrplot")) {
  install.packages("corrplot")
}
library("corrplot")

# QuGAcomp
if (!require("QuGAcomp")) {
  # GenomicRanges from Bioconductor
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  # Install curl package
  install.packages("curl")
  library("curl")
  # QuGAcomp from source file
  qugacomp.package.url <- "https://github.com/dritoshi/QuGAcomp/raw/master/QuGAcomp_0.99.2.tar.gz"
  qugacomp.package.loc <- file.path(".", basename(qugacomp.package.url))
  curl_download(qugacomp.package.url, qugacomp.package.loc)
  install.packages(qugacomp.package.loc)
}
library("QuGAcomp")

#
# Download files
#

# TODO: download genome info file from github, and run script to get bed files


#
# List up IDs of target experiments
#
experiments <- c("DRX013180", "SRX003937", "SRX022562")

#
# Path to the data directory
#
bed.data.dir <- file.path(".", script.basename, "..", "data", system("ls -t data | head -1", intern=TRUE), "bed")
genome.length.file <- file.path(bed.data.dir, "hg19.info")

#
# Load files and binning peak data
#
for (exp in experiments) {
  bed.file <- file.path(bed.data.dir, paste(exp, "bed", sep="."))
  gr <- loadBedFile(bed.file, genome.length.file)
  fat <- fat(gr, 200)
  unistd <- unifyStrand(fat)
  cov <- coverage(unistd)
  bin500 <- lapply(cov, function(x)rleBinning(x, 500))
  bin500 <- flatRleList(bin500)
  assign(paste(exp, "bin500", sep="."), bin500)
}

#
# QuGAcomp for all pairs
#
for (i in 1:NROW(experiments)) {
  left <- experiments[i]
  for (j in i:NROW(experiments)) {
    if (j > i) {
      right <- experiments[j]
      assign(
        paste("quga", left, right, sep="."),
        qugacomp(
          get(paste(left, "bin500", sep=".")),
          get(paste(right, "bin500", sep="."))
        )
      )
    }
  }
}

#
# Plot
#

# Prepare matrix
num <- NROW(experiments)
mat.cor <- matrix(0, nrow=num, ncol=num)
rownames(mat.cor) <- experiments
colnames(mat.cor) <- experiments
diag(mat.cor) <- rep(1,num)

# Calculate Pearson's correlation coefficient
for (i in 1:NROW(experiments)) {
  left <- experiments[i]
  for (j in i:NROW(experiments)) {
    if (j > i) {
      right <- experiments[j]
      mat.cor[i, j] <- mat.cor[j, i] <- pearsonCoef(get(paste("quga", left, right, sep=".")))
    }
  }
}

# Draw plot and save to pdf
mat.cor.max <- max(mat.cor[upper.tri(mat.cor, diag=F)])
mat.cor.min <- min(mat.cor[upper.tri(mat.cor, diag=F)])

print(mat.cor.max)
print(mat.cor.min)

print(mat.cor)

pdf("corrplot.pdf", width=5, height=5)
corrplot(mat.cor,
  method="circle",
  diag=FALSE,
  cl.lim=c(mat.cor.min, 1),
  cl.ratio=0.5
)
dev.off()
