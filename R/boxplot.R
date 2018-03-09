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
# Load matrix all x all
#
all.matrix.rds.file <- file.path(rds.dir, "all.matrix.rds")
if (!file.exists(all.matrix.rds.file)) {
  source(file.path(".", script.basename, "qugacomp.R"))
}

all.mat <- readRDS(all.matrix.rds.file)

#
# Create submatrix of in/out of the same TF for each experiment
#
tfs.vec <- readLines(file.path(metadata.dir, "tfs_downloaded.txt"))
metadata.df <- read.delim(file.path(metadata.dir, "exps_downloaded.reduced.tsv"), header=FALSE)
colnames(metadata.df) <- c("expid","tf","celltypeclass","celltype")

# for testing with small number of items
tfs.vec <- tfs.vec[1:5]
tfs.vec <- tfs.vec[!is.na(tfs.vec[summary(metadata.df$tf)[tfs.vec] > 2])]

all.mean.sd <- pforeach(i = 1:NROW(tfs.vec), .combine=rbind) ({
  tf = tfs.vec[i]
  exps.sameTF <- metadata.df[metadata.df$tf == tf,]$expid
  exps.diffTF <- metadata.df[metadata.df$tf != tf,]$expid

  foreach(exp_i = 1:NROW(exps.sameTF), .combine=rbind) %do% {
    expid <- as.character(exps.sameTF[exp_i])
    exps.sameButNotMe <- exps.sameTF[exps.sameTF != expid]

    vs.sameTF <- all.mat[rownames(all.mat) %in% exps.sameButNotMe, colnames(all.mat) == expid]
    vs.diffTF <- all.mat[rownames(all.mat) %in% exps.diffTF, colnames(all.mat) == expid]

    rbind(
      c(expid, tf, "same", mean(vs.sameTF), sd(vs.sameTF)),
      c(expid, tf, "diff", mean(vs.diffTF), sd(vs.diffTF))
    )
  }
})
all.mean.sd <- as.data.frame(all.mean.sd)
colnames(all.mean.sd) <- c("expid", "tf", "vs", "mean", "sd")

pforeach(i = 1:NROW(tfs.vec)) ({
  tf = tfs.vec[i]
  df <- all.mean.sd[all.mean.sd$tf == tf,]

  df$mean <- as.numeric(as.character(df$mean))
  df$sd <- as.numeric(as.character(df$sd))

  # Output PDF file path
  tf.output.pdf.name <- paste(tf, "boxplot", "pdf", sep=".")
  tf.output.pdf.path <- file.path(tfs.boxplot.dir, tf.output.pdf.name)

  p <- ggplot(df, aes(x = as.factor(expid)))
  p <- p + geom_boxplot(aes(
    colour = vs,
    lower = mean - sd,
    upper = mean + sd,
    middle = mean,
    ymin = mean - 3*sd,
    ymax = mean + 3*sd,
    ),
    stat="identity"
  )
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file=tf.output.pdf.path, plot=p, width=1080, unit="mm")
})

# Get the list of TFs
# tfs.vec <- readLines(file.path(metadata.dir, "tfs_downloaded.txt"))
#
# tfs.mean.sd.mat <- pforeach(i = 1:NROW(tfs.vec), .combine=rbind) ({
#   tf <- tfs.vec[i]
#
#   # Create or load matrix from RDS
#   tf.rds.file <- file.path(tfs.rds.dir, paste(tf, "rds", sep="."))
#   if(!file.exists(tf.rds.file)) {
#     c(tf, NA, NA)
#   } else {
#     tf.mat <- readRDS(tf.rds.file)
#     tf.values.vec <- tf.mat[upper.tri(tf.mat, diag=FALSE)]
#     c(tf, mean(tf.values.vec), sd(tf.values.vec))
#   }
# })
# tfs.mean.sd.df <- as.data.frame(tfs.mean.sd.mat)
# colnames(tfs.mean.sd.df) <- c("TF", "mean", "sd")
#
# tfs.mean.sd.df <- tfs.mean.sd.df[!is.na(tfs.mean.sd.df$sd),]
# tfs.mean.sd.df$mean <- as.numeric(as.character(tfs.mean.sd.df$mean))
# tfs.mean.sd.df$sd <- as.numeric(as.character(tfs.mean.sd.df$sd))
#
# p <- ggplot(tfs.mean.sd.df, aes(x = as.factor(TF)))
# p <- p + geom_boxplot(aes(
#   lower = mean - sd,
#   upper = mean + sd,
#   middle = mean,
#   ymin = mean - 3*sd,
#   ymax = mean + 3*sd,
#   ),
#   stat="identity"
# )
# ggsave(file="test.pdf", plot=p)
