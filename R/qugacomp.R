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
# Create directories to save results
#
data.dir <- file.path(".", script.basename, "..", "data", system(paste("ls -t ", file.path(".", script.basename, "..", "data"), " | head -1", sep=""), intern=TRUE))
metadata.dir <- file.path(data.dir, "metadata")

rds.dir <- file.path(data.dir, "rds")
tfs.rds.dir <- file.path(rds.dir, "eachTF")
dir.create(tfs.rds.dir, showWarnings=FALSE, recursive=TRUE)

matrix.dir <- file.path(data.dir, "matrix")
tfs.matrix.dir <- file.path(matrix.dir, "eachTF")
dir.create(tfs.matrix.dir, showWarnings=FALSE, recursive=TRUE)


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

all.matrix.rds.file <- file.path(rds.dir, "all.matrix.rds")
if (!file.exists(all.matrix.rds.file)) {
  # Load bin file, run bed2bin if not rds file exists
  bin.rds.file <- file.path(rds.dir, "bin.rds")
  if (!file.exists(bin.rds.file)) {
    source(file.path(".", script.basename, "bed2bin.R"))
  }
  bin.list <- readRDS(bin.rds.file)
  all.mat <- createCorrMatrix(experiments, bin.list)
  saveRDS(all.mat, all.matrix.rds.file)

  # Save matrix
  all.mat.file <- file.path(matrix.dir, "all.matrix.tsv")
  write.table(all.mat, file=all.mat.file, sep="\t", quote=FALSE)
}

#
# Save submatrix for each tf
#

# Get the list of TFs
tfs.vec <- readLines(file.path(metadata.dir, "tfs_downloaded.txt"))

# Get the metadata table
metadata <- read.delim(file.path(metadata.dir, "exps_downloaded.reduced.tsv"), header=FALSE)

# Create matrix for each TF
x <- pforeach(i = 1:NROW(tfs.vec)) ({
  # Create submatrix
  tf <- tfs.vec[i]

  # Create or load matrix from RDS
  tf.rds.file <- file.path(tfs.rds.dir, paste(tf, "rds", sep="."))
  if (!file.exists(tf.rds.file)) {
    tf.exps <- metadata[metadata$V2 == tf,]$V1
    tf.mat <- all.mat[rownames(all.mat) %in% tf.exps, colnames(all.mat) %in% tf.exps]
    saveRDS(tf.mat, tf.rds.file)

    # Save in tsv
    tf.mat.file.path <- file.path(tfs.matrix.dir, paste(tf, "matrix.tsv", sep="."))
    write.table(tf.mat, file=tf.mat.file.path, sep="\t", quote=FALSE)
  }
})

print(paste("Correlation coefficient calculated. Matrix saved at", matrix.dir))
