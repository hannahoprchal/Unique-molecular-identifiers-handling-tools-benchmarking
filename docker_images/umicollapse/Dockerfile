FROM debian:bullseye-20211220-slim AS build-base


ARG samtools_version=1.16.1

ARG fastp_version=0.23.1

# run update and install necessary tools
RUN apt-get update -y --fix-missing && apt-get install -y \
wget build-essential zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev python make zlib1g vim default-jre python3 python3-pip git openjdk-11-jdk curl


# UMICOLLAPSE
RUN git clone https://github.com/Daniel-Liu-c0deb0t/UMICollapse.git && \
    cd UMICollapse && \
    mkdir lib && cd lib && \
    wget https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar && \
    cd .. && \
    ./build.sh && \
    mv umicollapse /usr/bin

# FASTP
ADD http://opengene.org/fastp/fastp.${fastp_version} .
RUN mv fastp.${fastp_version} fastp && \
    chmod a+x ./fastp && \
    mv fastp /usr/bin

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

# BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && make && \
    mv bwa /usr/bin

