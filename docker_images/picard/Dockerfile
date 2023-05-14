#similar: https://hub.docker.com/r/broadinstitute/picard/dockerfile
FROM debian:bullseye-20211220-slim AS build-base

ARG picard_version=2.27.5
ENV PICARD /usr/local/bin/picard.jar

# run update and install necessary tools
RUN apt-get update -y --fix-missing && apt-get install -y \
wget build-essential zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev python make zlib1g vim default-jre python3 python3-pip git openjdk-11-jdk curl


# PICARD

RUN wget https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar -O /usr/local/bin/picard.jar && \
 chmod 0644 /usr/local/bin/picard.jar