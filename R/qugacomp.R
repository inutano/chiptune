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
# Timestamp
#
today <- format(as.POSIXlt(Sys.time(), "GMT-9"), "%Y%m%d-%H%M")

#
# Package install and load
#
source(file.path(".", script.basename, "setup.R"))

#
# Load bin data
#
source(file.path(".", script.basename, "bed2bin.R"))

#
# Create a directory to save results
#
analysis.dir <- file.path(data.dir, "analysis", today)
dir.create(analysis.dir, showWarnings=FALSE, recursive=TRUE)

corrplot.dir <- file.path(analysis.dir, "corrplot")
dir.create(corrplot.dir, showWarnings=FALSE, recursive=TRUE)

matrix.dir <- file.path(analysis.dir, "matrix")
dir.create(matrix.dir, showWarnings=FALSE, recursive=TRUE)


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

corr.matrix.rds.file <- file.path(rds.dir, "corr.matrix.rds")
if (!file.exists(corr.matrix.rds.file)) {
  bin.list <- readRDS(bin.rds.file)
  saveRDS(createCorrMatrix(experiments, bin.list), corr.matrix.rds.file)
}
mat.cor <- readRDS(corr.matrix.rds.file)

mat.cor.file <- file.path(matrix.dir, "all.matrix.tsv")
write.table(mat.cor, file=mat.cor.file, sep="\t", quote=FALSE)
print(paste("Correlation calculated. Matrix saved at", mat.cor.file))

#
# Plot all for all
#

# max/min
mat.cor.max <- max(mat.cor[upper.tri(mat.cor, diag=FALSE)])
mat.cor.min <- min(mat.cor[upper.tri(mat.cor, diag=FALSE)])

# Draw plot and save to pdf
output.pdf.name <- paste("all", "corrplot", "pdf", sep=".")
output.pdf.path <- file.path(corrplot.dir, output.pdf.name)
pdf(output.pdf.path)
corrplot(mat.cor,
  method="circle",
  diag=FALSE,
  cl.lim=c(mat.cor.min, 1),
  cl.ratio=0.2
)
invisible(dev.off())

print(paste("The plot saved at", output.pdf.path))

#
# Plot per TF
#

# Get the list of TFs
tfs.vec <- readLines(file.path(data.dir, "tfs_download.txt"))

# Get the metadata table
metadata <- read.delim(file.path(data.dir, "data.reduced.tsv"), header=FALSE)

# Create matrix for each TF
x <- pforeach(i = 1:NROW(tfs.vec)) ({
  # Create submatrix
  tf <- tfs.vec[i]
  tf.exps <- metadata[metadata$V2 == tf,]$V1
  tf.mat <- mat.cor[rownames(mat.cor) %in% tf.exps, colnames(mat.cor) %in% tf.exps]

  # Save in tsv
  tf.mat.file.path <- file.path(matrix.dir, paste(tf, "matrix.tsv", sep="."))
  write.table(tf.mat, file=tf.mat.file.path, sep="\t", quote=FALSE)

  # Plot
  tf.mat.max <- max(tf.mat[upper.tri(tf.mat, diag=FALSE)])
  tf.mat.min <- min(tf.mat[upper.tri(tf.mat, diag=FALSE)])

  # Output PDF file path
  tf.output.pdf.name <- paste(tf, "corrplot", "pdf", sep=".")
  tf.output.pdf.path <- file.path(corrplot.dir, tf.output.pdf.name)

  # Draw plot and save
  pdf(tf.output.pdf.path)
  corrplot(tf.mat,
    method="circle",
    diag=FALSE,
    cl.lim=c(tf.mat.min, 1),
    cl.ratio=0.2
  )
  invisible(dev.off())
})
