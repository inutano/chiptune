FROM bioconductor/release_base2:R3.6.2_Bioc3.10
RUN apt-get update -y && apt-get install -y curl lftp
ADD . /app
RUN cd /app && Rscript --vanilla R/setup.R
CMD ["bash"]
