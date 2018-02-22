# chiptune/R/bed2bin.R: load bed files and binning
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
# Path to the data directory (use latest-downloaded data)
#
data.dir <- file.path(".", script.basename, "..", "data", system("ls -t data | head -1", intern=TRUE))
bed.data.dir <- file.path(data.dir, "bed")
genome.length.file <- file.path(bed.data.dir, "hg19.info")

rds.dir <- file.path(data.dir, "rds")
dir.create(rds.dir, showWarnings=FALSE, recursive=TRUE)

#
# List up IDs of target experiments (all the bed files in the directory)
#
#experiments <- c("DRX013180", "SRX003937", "SRX022562")
experiments <- system(paste("ls", bed.data.dir, "| grep '.bed$' | sed -e 's:.bed$::'"), intern=TRUE)
exps.num <- NROW(experiments)

#
# Configurations
#
fat.num <- 200
bin.length <- 500

#
# Load files and binning peak data
#

loadAndBin <- function(bed.file.path, glength.file, fat.num){
  flatRleList(
    lapply(
      coverage(
        unifyStrand(
          fat(
            loadBedFile(
              bed.file.path,
              glength.file
            ),
            fat.num
          )
        )
      ),
      function(x)rleBinning(x, bin.length)
    )
  )
}

pLoadAndBin <- function(data.dir, exps.vec, glength.file, fat.num){
  pforeach (i = 1:NROW(exps.vec)) ({
    list(
      loadAndBin(
        file.path(data.dir, paste(exps.vec[i], "bed", sep=".")),
        glength.file,
        fat.num
      )
    )
  })
}

bin.rds.file <- file.path(rds.dir, "bin.rds")
if (!file.exists(bin.rds.file)) {
  saveRDS(pLoadAndBin(bed.data.dir, experiments, genome.length.file, fat.num), bin.rds.file)
  print(paste(exps.num, "bed files loaded and saved at", bin.rds.file))
}
