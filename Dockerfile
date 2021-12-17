FROM ubuntu
MAINTAINER Sagnik Banerjee <sagnikbanerjee15@gmail.com>

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

ARG ABRIDGE_VERSION=1.0.0
ARG ZPAQ_VERSION=715
ARG SAMTOOLS_VERSION=1.14

# Update base image and install software
RUN apt-get -y update
RUN apt-get -y install git python3 less vim wget time zlib1g zlib1g-dev lzma-dev libncurses5-dev libcurl4-nss-dev liblzma-dev libncursesw5-dev make unzip zip build-essential gcc g++ cmake ca-certificates libbz2-dev xz-utils htop autoconf automake binutils bison flex gettext libtool make patch pkg-config p7zip-full p7zip-rar 
RUN apt-get clean all

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

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
# Downloading the current git repo - change this to a specific version later
RUN cd /software && \
	git clone https://github.com/sagnikbanerjee15/Abridge.git &&\
	cd Abridge/src &&\
	make &&\
	chmod -R 777 /software/Abridge

ENV PATH /software/Abridge:/software/Abridge/src:/software/Abridge/scripts:${PATH}


#./scripts/run_abridge.py  --keep_intermediate_error_files  --compress  --genome example/arath_1.fa  --inputsamfilenames example/SRR13711353_SE_100000.sam  --output_directory example/SRR13711353_SE_100000_compress_level_1_ignore_scores_0_ignore_quality_scores_0_ignore_soft_clippings_0_ignore_mismatches_0_ignore_unmapped_reads_0_save_all_quality_scores_0_save_exact_quality_scores_0  --level 1 