#!/bin/bash
# require: curl, lftp

#
# Variables
#

# qvalue of macs2 peak call, bed files are avaiable for 05 (1E-05), 10 (1E-10), 20 (1E-20)
QVAL=20

#
# Functions
#
date_cmd() {
  local arg="${1}"
  if [[ "$(uname)" == 'Darwin' ]]; then
    gdate ${arg}
  elif [[ "$(expr substr $(uname -s) 1 5)" == 'Linux' ]]; then
    date ${arg}
  else
    echo "Your platform ($(uname -a)) is not supported." 2>dev/null
    exit 1
  fi
}

#
# Global variables
#
BASEDIR="$(cd $(dirname ${0}) && pwd -P)"

data_dir="${BASEDIR}/../data/$(date_cmd +%Y%m%d-%H%M)"
mkdir -p "${data_dir}"

bed_dir="${data_dir}/bed"
mkdir -p "${bed_dir}"

#
# Get genome length file from github
#
hg19_chrominfo_url="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz"
curl "${hg19_chrominfo_url}" | gunzip | awk -F '\t' '$1 ~ /^chr(.|..)$/ { print $1 "\t" $2 }' > "${bed_dir}/hg19.info"

#
# Prepare the list of experiments
#

# The original list of processed experiments on the remote FTP server
experimentList_url="http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab"

# Select ChIP-seq experiments of human TFs
datalist_path="${data_dir}/data.tsv"
curl -s ${experimentList_url} | awk -F'\t' '$2 == "hg19" && $3 == "TFs and others"' > "${datalist_path}"

# The TF experiment ranking
ranking_path="${data_dir}/ranking.txt"
cat "${datalist_path}" | cut -f 4 | sort | uniq -c | sort -nr > "${ranking_path}"

# Top 10 TFs
top10tfs_path="${data_dir}/top10.txt"
cat "${ranking_path}" | awk '$1 > 4' | awk '$0=$2' | head -10 > "${top10tfs_path}"

# Get each 3 experiment IDs of top 10 TFs
top10_each3expids_path="${data_dir}/top10_each3expids.tsv"
cat "${top10tfs_path}" | while read tf; do
  cat "${datalist_path}" | awk -v tf="${tf}" -F '\t' '$4 == tf' | head -3
done > "${top10_each3expids_path}"

# Create lftp script and exec to download bed files
lftp_script_path="${data_dir}/download_bed.lftp"
FTP_base="ftp://ftp.biosciencedbc.jp/data/chip-atlas/data/hg19/eachData/bed${QVAL}/"
cat "${top10_each3expids_path}" |\
  awk -F '\t' -v qval="${QVAL}" -v ftp="${FTP_base}" -v outdir="${bed_dir}" 'BEGIN{ print "open " ftp } { print "pget -n 8 -O " outdir " " $1 "." qval ".bed" }' > "${lftp_script_path}"
lftp -f "${lftp_script_path}"

# edit bed files
ls "${bed_dir}" | grep 'bed$' | while read bed; do
  id=$(echo "${bed}" | sed -e 's:\..*$::g')
  cat "${bed_dir}/${bed}" | awk -F '\t' 'BEGIN{ OFS="\t"; print "header" } $6 = "+"' \
    > "${bed_dir}/${id}.bed"
  rm -f "${bed_dir}/${bed}"
done

# Remove bed files with only 0/1 peak
find "${bed_dir}" | awk '/\.bed$/' | xargs wc -l | awk '$1 < 3 { print $2 }' | xargs rm
