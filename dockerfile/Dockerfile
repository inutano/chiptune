#
# Dockerfile for inutano/chiptune
#

FROM bioconductor/release_base2:R3.4.3_Bioc3.6

RUN apt-get update -y && apt-get install -y curl lftp

RUN git clone https://github.com/inutano/chiptune /chiptune && \
    cd /chiptune && \
    Rscript --vanilla R/setup.R

CMD ["bash"]
