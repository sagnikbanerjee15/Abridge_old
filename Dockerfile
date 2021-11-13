FROM ubuntu
MAINTAINER Sagnik Banerjee <sagnikbanerjee15@gmail.com>

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive 

ARG ABRIDGE_VERSION=1.0.0
ARG ZPAQ_VERSION=715

# Update base image and install software
RUN apt-get -y update
RUN apt-get install -y --no-install-recommends git python3 less vim wget time build-essential zlib1g-dev libncurses5-dev gcc g++ make unzip zip cmake ca-certificates
RUN apt-get clean

# Create directories for installation
RUN mkdir /software 

# Download and install dependencies
RUN cd /software &&\
	git clone https://github.com/google/brotli.git && \
	cd brotli && \
	mkdir out && \
	cd out && \
	../configure-cmake &&\
	make && \
	make test && \
	make install
	
RUN cd /software &&\
	mkdir zpaq && \
	cd zpaq &&\
	wget http://mattmahoney.net/dc/zpaq${ZPAQ_VERSION}.zip && \
	unzip zpaq715.zip && \
	make && \
	make install


# Downloading the current git repo - change this to a specific version later
RUN cd /software && \
	git clone https://github.com/sagnikbanerjee15/Abridge.git && \
	cd Abridge/src && \
	make && \
	make install
	
ENV PATH /softwares/Abridge:${PATH}