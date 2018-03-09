# chiptune/R/boxplot.R: draw boxplot from correlation coefficient matrix
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

analysis.dir <- file.path(data.dir, "analysis", today)
boxplot.dir <- file.path(analysis.dir, "boxplot")
tfs.boxplot.dir <- file.path(boxplot.dir, "eachTF")
dir.create(tfs.boxplot.dir, showWarnings=FALSE, recursive=TRUE)

#
# Create data frame for each TF
#

# Get the list of TFs
tfs.vec <- readLines(file.path(metadata.dir, "tfs_downloaded.txt"))

tfs.mean.sd.mat <- pforeach(i = 1:NROW(tfs.vec), .combine=rbind) ({
  tf <- tfs.vec[i]

  # Create or load matrix from RDS
  tf.rds.file <- file.path(tfs.rds.dir, paste(tf, "rds", sep="."))
  if(!file.exists(tf.rds.file)) {
    c(tf, NA, NA)
  } else {
    tf.mat <- readRDS(tf.rds.file)
    tf.values.vec <- tf.mat[upper.tri(tf.mat, diag=FALSE)]
    c(tf, mean(tf.values.vec), sd(tf.values.vec))
  }
})
tfs.mean.sd.df <- as.data.frame(tfs.mean.sd.mat)
colnames(tfs.mean.sd.df) <- c("TF", "mean", "sd")

tfs.mean.sd.df <- tfs.mean.sd.df[!is.na(tfs.mean.sd.df$sd),]
tfs.mean.sd.df$mean <- as.numeric(as.character(tfs.mean.sd.df$mean))
tfs.mean.sd.df$sd <- as.numeric(as.character(tfs.mean.sd.df$sd))

p <- ggplot(tfs.mean.sd.df, aes(x = as.factor(TF)))
p <- p + geom_boxplot(aes(
  lower = mean - sd,
  upper = mean + sd,
  middle = mean,
  ymin = mean - 3*sd,
  ymax = mean + 3*sd,
  ),
  stat="identity"
)
ggsave(file="test.pdf", plot=p)
