# chiptune/R/search.R: find similar ChIP-seq experiments
# usage:
#   Rscript search.R <input.bed>
# Input bed file must have one-line header, and the strand column (6th) should not be "." (no-strand)
# Tested on bioconductor/release_base2:R3.4.3_Bioc3.6

#
# Parse argument
#
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Require input bed file to search similar experiments.")
} else {
  if (!file.exists(args[1])) {
    stop("Input bed file not found.")
  }
}
input.bed.path <- args[1]


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
# Create a directory to save results
#
data.dir <- file.path(".", script.basename, "..", "data", system(paste("ls -t ", file.path(".", script.basename, "..", "data"), " | head -1", sep=""), intern=TRUE))
metadata.dir <- file.path(data.dir, "metadata")
rds.dir <- file.path(data.dir, "rds")

search.dir <- file.path(data.dir, "search", today)
dir.create(search.dir, showWarnings=FALSE, recursive=TRUE)

#
# Package install and load
#
source(file.path(".", script.basename, "setup.R"))

#
# Source bed2bin and load rds file
#
source(file.path(".", script.basename, "bed2bin.R"))
bin.rds.file <- file.path(rds.dir, "bin.rds")
bin.list <- readRDS(bin.rds.file)

#
# Load given bed data and bin
#
input.bin <-loadAndBin(input.bed.path, genome.length.file, fat.num, bin.length)

#
# Calculate correlation and append metadata
#

# Get the metadata table
metadata <- read.delim(file.path(metadata.dir, "exps_downloaded.reduced.tsv"), header=FALSE)
metadata <- data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)

# Calculate/Append
result <- pforeach(i = 1:NROW(bin.list), .combine=rbind) ({
  list(
    pearsonCoef(qugacomp(input.bin, bin.list[[i]])),
    metadata[i,2],
    metadata[i,3],
    metadata[i,4]
  )
})

# Put row/col names and sort by coefficient
result <- as.data.frame(result)
rownames(result) <- metadata$V1
colnames(result) <- c("Correlation", "Antigen", "CellTypeClass", "CellType")

result <- transform(result, Correlation = as.numeric(Correlation), Antigen = as.character(Antigen), CellTypeClass = as.character(CellTypeClass), CellType = as.character(CellType))

result <- result[order(-result$Correlation),]
result <- cbind(ExpID = rownames(result), result)

# Save the result in a tsv file with rownames
result.file <- file.path(search.dir, paste(basename(input.bed.path), "tsv", sep="."))
write.table(result, file=result.file, sep="\t", quote=FALSE, row.names=FALSE)
print(result)
print(paste("Search result saved at", result.file))
