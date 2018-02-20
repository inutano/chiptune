# chiptune/R/doanload.R: download bed files
# Tested on bioconductor/release_base2:R3.4.3_Bioc3.6

#
# Get path to the directory of this script
#
initial.options <- commandArgs(trailingOnly=FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#
# Download files
#
system(paste("bash ", file.path(".", script.basename, "..", "bin", "update_datalist.sh")))
