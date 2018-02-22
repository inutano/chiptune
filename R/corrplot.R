# chiptune/R/corrplot.R: draw correlation plot from correlation coefficient matrix
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
# Create directories to save results
#
data.dir <- file.path(".", script.basename, "..", "data", system("ls -t data | head -1", intern=TRUE))

rds.dir <- file.path(data.dir, "rds")
tfs.rds.dir <- file.path(rds.dir, "eachTF")

analysis.dir <- file.path(data.dir, "analysis", today)
corrplot.dir <- file.path(analysis.dir, "corrplot")
tfs.corrplot.dir <- file.path(corrplot.dir, "eachTF")
dir.create(tfs.corrplot.dir, showWarnings=FALSE, recursive=TRUE)

#
# Load matrix
#
all.matrix.rds.file <- file.path(rds.dir, "all.matrix.rds")
if (!file.exist(all.matrix.rds.file)) {
  source(file.path(".", script.basename, "qugacomp.R"))
}

all.mat <- readRDS(all.matrix.rds.file)
all.mat.max <- max(all.mat[upper.tri(all.mat, diag=FALSE)])
all.mat.min <- min(all.mat[upper.tri(all.mat, diag=FALSE)])

#
# Plot all for all
#
output.pdf.name <- paste("all", "corrplot", "pdf", sep=".")
output.pdf.path <- file.path(corrplot.dir, output.pdf.name)
pdf(output.pdf.path)
corrplot(all.mat,
  method="circle",
  diag=FALSE,
  cl.lim=c(all.mat.min, 1),
  cl.ratio=0.2
)
invisible(dev.off())

#
# Plot per TF
#

# Get the list of TFs
tfs.vec <- readLines(file.path(data.dir, "tfs_download.txt"))

# Get the metadata table
metadata <- read.delim(file.path(data.dir, "data.reduced.tsv"), header=FALSE)

# Draw plot for each TF
x <- pforeach(i = 1:NROW(tfs.vec)) ({
  tf <- tfs.vec[i]

  # Create or load matrix from RDS
  tf.rds.file <- file.path(tfs.rds.dir, paste(tf, "rdfs", sep="."))
  if(!file.exist(tf.rds.file)) {
    print(paste("Data for", tf, "is missing in", tfs.rds.dir))
  } else {
    tf.mat <- readRDS(tf.rds.file)
    tf.mat.max <- max(tf.mat[upper.tri(tf.mat, diag=FALSE)])
    tf.mat.min <- min(tf.mat[upper.tri(tf.mat, diag=FALSE)])

    # Output PDF file path
    tf.output.pdf.name <- paste(tf, "corrplot", "pdf", sep=".")
    tf.output.pdf.path <- file.path(tfs.corrplot.dir, tf.output.pdf.name)

    # Draw plot and save
    pdf(tf.output.pdf.path)
    corrplot(tf.mat,
      method="circle",
      diag=FALSE,
      cl.lim=c(tf.mat.min, 1),
      cl.ratio=0.2
    )
    invisible(dev.off())
  }
})

print(paste("The correlation plots were saved at", corrplot.dir))
