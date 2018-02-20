# chiptune/R/qugacomp.R: exec quantitative genome annotation comparison (QuGAcomp) https://github.com/dritoshi/QuGAcomp
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
source(file.path(".", script.basename, "setup.R"))

#
# Load bin data
#
source(file.path(".", script.basename, "bed2bin.R"))

#
# QuGAcomp for all pairs
#
createCorrMatrix <- function(exps.vec, bin.list){
  mat.length <- NROW(exps.vec)

  # Comparison and calculation for each cell of the matrix
  mat <- pforeach(i = 1:mat.length, .combine=cbind) ({
    foreach (j = 1:mat.length) %do% {
      if (i > j) {
        pearsonCoef(qugacomp(bin.list[[i]], bin.list[[j]]))
      }
    }
  })

  # Convert indices (list) to numeric
  mat <- matrix(sapply(mat, as.numeric), nrow=mat.length, ncol=mat.length)
  storage.mode(mat) <- "numeric"

  # Fill lower triangle
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  # Apply names and diagonal
  rownames(mat) <- exps.vec
  colnames(mat) <- exps.vec
  diag(mat) <- rep(1, mat.length)

  # Return matrix
  mat
}

corr.matrix.rds.file <- file.path(bed.data.dir, "corr.matrix.rds")
if (!file.exists(corr.matrix.rds.file)) {
  saveRDS(createCorrMatrix(experiments, bin.list), corr.matrix.rds.file)
}
mat.cor <- readRDS(corr.matrix.rds.file)

print("Correlation calculated. Preparing plotting..")

#
# Plot
#

# Draw plot and save to pdf
mat.cor.max <- max(mat.cor[upper.tri(mat.cor, diag=FALSE)])
mat.cor.min <- min(mat.cor[upper.tri(mat.cor, diag=FALSE)])

# Prepare directory to save corrplot
analysis.dir.path <- file.path(".", script.basename, "..", "analysis", "corrplot")
dir.create(analysis.dir.path, showWarnings=FALSE, recursive=TRUE)

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
invisible(dev.off())

print(paste("The plot saved at", output.pdf.path))
