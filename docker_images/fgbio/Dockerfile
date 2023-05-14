FROM debian:bullseye-20211220-slim AS build-base

ARG fgbio_version=2.1.0
ENV FGBIO /usr/local/bin/fgbio.jar

ARG samtools_version=1.16.1


# run update and install necessary tools
RUN apt-get update -y --fix-missing && apt-get install -y \
wget build-essential zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev python make zlib1g vim default-jre python3 python3-pip git openjdk-11-jdk curl


# SAMTOOLS

# Pulling SAMTools from its repository, unpacking the archive and installing
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
&& tar jxf samtools-${samtools_version}.tar.bz2 \
&& cd samtools-${samtools_version} \
&& ./configure --prefix=/samtools-${samtools_version}/ \
&& make \
&& make install \
&& cp bin/* /usr/local/bin/ \
&& cd .. \
&& rm -r samtools-${samtools_version} \
&& rm samtools-${samtools_version}.tar.bz2

# FGBIO

RUN wget https://github.com/fulcrumgenomics/fgbio/releases/download/${fgbio_version}/fgbio-${fgbio_version}.jar -O /usr/local/bin/fgbio.jar && \
 chmod 0644 /usr/local/bin/fgbio.jar

# BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && make && \
    mv bwa /usr/bin

