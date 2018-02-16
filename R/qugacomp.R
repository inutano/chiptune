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
update.packages(checkBuilt=TRUE, ask=FALSE, repos="https://cran.ism.ac.jp/")

# Corrplot
if (!require("corrplot")) {
  install.packages("corrplot", repos="https://cran.ism.ac.jp/")
}
library("corrplot")

# QuGAcomp
if (!require("QuGAcomp")) {
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
}
library("QuGAcomp")

# pforeach
if (!require("pforeach")) {
  install.packages("devtools", repos="https://cran.ism.ac.jp/")
  devtools::install_github("hoxo-m/pforeach")
}
library("pforeach")

#
# Download files
#
system(paste("bash ", file.path(".", script.basename, "..", "bin", "update_datalist.sh")))

#
# Path to the data directory
#
bed.data.dir <- file.path(".", script.basename, "..", "data", system("ls -t data | head -1", intern=TRUE), "bed")
genome.length.file <- file.path(bed.data.dir, "hg19.info")

#
# List up IDs of target experiments
#
#experiments <- c("DRX013180", "SRX003937", "SRX022562")
experiments <- system(paste("ls", bed.data.dir, "| grep '.bed$' | sed -e 's:.bed$::'"), intern=TRUE)
exps.num <- NROW(experiments)

#
# Load files and binning peak data
#
bin500.list <- pforeach (i = 1:exps.num) ({
  exp <- experiments[i]
  bed.file <- file.path(bed.data.dir, paste(exp, "bed", sep="."))
  list(
    flatRleList(
      lapply(
        coverage(
          unifyStrand(
            fat(
              loadBedFile(bed.file, genome.length.file), 200)
          )
        ), function(x)rleBinning(x, 500)
      )
    )
  )
})
print("All bed files loaded. Exec QuGAcomp and correlation calculation..")

#
# QuGAcomp for all pairs
#

# Comparison and calculation for each cell of the matrix
mat.cor <- pforeach(i = 1:exps.num, .combine=cbind) ({
  foreach (j = 1:exps.num) %do% {
    if (i > j) {
      pearsonCoef(qugacomp(bin500.list[[i]], bin500.list[[j]]))
    }
  }
})

# Configure matrix
mat.cor <- matrix(sapply(mat.cor, as.numeric), nrow=exps.num, ncol=exps.num)
storage.mode(mat.cor) <- "numeric"
mat.cor[lower.tri(mat.cor)] <- t(mat.cor)[lower.tri(mat.cor)]

rownames(mat.cor) <- experiments
colnames(mat.cor) <- experiments
diag(mat.cor) <- rep(1,exps.num)

print("Calculation done. Preparing plotting..")

#
# Plot
#

# Draw plot and save to pdf
mat.cor.max <- max(mat.cor[upper.tri(mat.cor, diag=F)])
mat.cor.min <- min(mat.cor[upper.tri(mat.cor, diag=F)])

analysis.dir.path <- file.path(".", script.basename, "..", "analysis", "corrplot")
today <- format(as.POSIXlt(Sys.time(), "GMT-9"), "%Y%m%d-%H%M")
output.pdf.name <- paste("corrplot", today, "pdf", sep=".")
output.pdf.path <- file.path(analysis.dir.path, output.pdf.name)
pdf(output.pdf.path)
corrplot(mat.cor,
  method="circle",
  diag=FALSE,
  cl.lim=c(mat.cor.min, 1),
  cl.ratio=0.2
)
dev.off()
