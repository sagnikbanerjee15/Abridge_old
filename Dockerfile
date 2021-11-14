FROM centos
MAINTAINER Sagnik Banerjee <sagnikbanerjee15@gmail.com>

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive 

ARG ABRIDGE_VERSION=1.0.0
ARG ZPAQ_VERSION=715
ARG SAMTOOLS_VERSION=1.14

# Update base image and install software
RUN yum -y update
RUN rpm -Uvh https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
RUN yum -y install git python3 less vim wget time zlib-devel ncurses-devel gcc gcc-c++ make unzip zip cmake ca-certificates bzip2-devel xz-devel libcurl-devel htop autoconf automake binutils bison flex gettext libtool make patch pkgconfig redhat-rpm-config rpm-build rpm-sign ctags elfutils patchutils p7zip p7zip-plugins 
RUN yum clean all

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
	unzip zpaq${ZPAQ_VERSION}.zip && \
	make && \
	make install

RUN cd /software && \
	mkdir samtools && \
	cd samtools && \
	wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	cd samtools-${SAMTOOLS_VERSION} && make && make install &&\
	cd .. && rm samtools-${SAMTOOLS_VERSION}.tar.bz2

# Downloading the current git repo - change this to a specific version later
RUN cd /software && \
	git clone https://github.com/sagnikbanerjee15/Abridge.git && \
	cd Abridge/src && \
	make && \
	make install
	
ENV PATH /software/Abridge:/software/Abridge/scripts:${PATH}