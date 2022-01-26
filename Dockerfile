FROM rust
MAINTAINER Sagnik Banerjee <sagnikbanerjee15@gmail.com>

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

# Update base image and install software
RUN apt-get -y update
RUN apt-get -y install git python3 less vim wget time zlib1g zlib1g-dev lzma-dev \
	libncurses5-dev libcurl4-nss-dev liblzma-dev libncursesw5-dev make unzip zip build-essential \
	gcc g++ cmake ca-certificates libbz2-dev xz-utils htop autoconf automake binutils bison flex \
	gettext libtool make patch pkg-config p7zip-full p7zip python
RUN apt-get clean all

# Create directories for installation
RUN mkdir /software 

###################################################################################################################################################################################################
# STAR
###################################################################################################################################################################################################

ARG STAR_VERSION=2.7.9a
RUN mkdir -p /software/STAR
RUN cd /software/STAR && wget --no-check-certificate https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.zip && unzip ${STAR_VERSION}.zip 
RUN cd /software/STAR/STAR-${STAR_VERSION}/source && make STAR
ENV PATH="${PATH}:/software/STAR/STAR-${STAR_VERSION}/source"

###################################################################################################################################################################################################
# HiSAT2
###################################################################################################################################################################################################

ARG HISAT2_VERSION=2.2.1
RUN mkdir -p /software/HISAT2 && cd /software/HISAT2 && \
	wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
	unzip download
ENV PATH="${PATH}:/software/HISAT2/hisat2-${HISAT2_VERSION}"
	
###################################################################################################################################################################################################
# NCBI SRA toolkit
###################################################################################################################################################################################################

ARG SRATOOLKIT_VERSION=2.11.3
RUN mkdir -p /software/sratoolkit && cd /software/sratoolkit && \
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKIT_VERSION}/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz && \
	tar -xvzf sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz
ENV PATH="${PATH}:/software/sratoolkit/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64/bin"
	
###################################################################################################################################################################################################
# Bowtie1
###################################################################################################################################################################################################

ARG BOWTIE1_VERSION=1.3.1
RUN mkdir /software/bowtie1 && cd /software/bowtie1 && \
	wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-${BOWTIE1_VERSION}-linux-x86_64.zip && \
	unzip bowtie-${BOWTIE1_VERSION}-linux-x86_64.zip
ENV PATH="${PATH}:/software/bowtie1//bowtie-${BOWTIE1_VERSION}-linux-x86_64"

###################################################################################################################################################################################################
# Bowtie2
###################################################################################################################################################################################################

ARG BOWTIE2_VERSION=2.4.5
RUN mkdir /software/bowtie2 && cd /software/bowtie2 && \
	wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
	unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
ENV PATH="${PATH}:/software/bowtie2/bowtie2-${BOWTIE2_VERSION}-linux-x86_64"

###################################################################################################################################################################################################
# Brotli
###################################################################################################################################################################################################
RUN mkdir -p /software/brotli cd /software/brotli &&\
	git clone https://github.com/google/brotli.git && \
	cd brotli && \
	mkdir out && \
	cd out && \
	../configure-cmake &&\
	make && \
	make test && \
	make install

###################################################################################################################################################################################################
# ZPAQ
###################################################################################################################################################################################################
ARG ZPAQ_VERSION=715
RUN mkdir -p /software/zpaq && cd /software/zpaq &&\
	mkdir zpaq && \
	cd zpaq &&\
	wget http://mattmahoney.net/dc/zpaq${ZPAQ_VERSION}.zip && \
	unzip zpaq${ZPAQ_VERSION}.zip && \
	make && \
	make install

###################################################################################################################################################################################################
# Samtools
###################################################################################################################################################################################################
ARG SAMTOOLS_VERSION=1.14
RUN mkdir -p /software/samtools && cd /software/samtools && \
	mkdir samtools && \
	cd samtools && \
	wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	cd samtools-${SAMTOOLS_VERSION} && make && make install &&\
	cd .. && rm samtools-${SAMTOOLS_VERSION}.tar.bz2

###################################################################################################################################################################################################
# ABRIDGE
###################################################################################################################################################################################################

ARG ABRIDGE_VERSION=1.0.0
RUN mkdir -p /software/abridge && cd /software/abridge && \
	wget https://github.com/sagnikbanerjee15/Abridge/archive/refs/tags/ABRIDGE_v${ABRIDGE_VERSION}.tar.gz &&\
	tar -xvzf ABRIDGE_v${ABRIDGE_VERSION}.tar.gz && \
	cd Abridge-ABRIDGE_v${ABRIDGE_VERSION}/src &&\
	make &&\
	chmod -R 777 /software/abridge/Abridge-ABRIDGE_v${ABRIDGE_VERSION}
ENV PATH="${PATH}:/software/Abridge-ABRIDGE_v${ABRIDGE_VERSION}:/software/Abridge-ABRIDGE_v${ABRIDGE_VERSION}/src:/software/Abridge-ABRIDGE_v${ABRIDGE_VERSION}/scripts"

###################################################################################################################################################################################################
# LCQS	
###################################################################################################################################################################################################

RUN mkdir /software/LCQS && cd /software/LCQS && \
	git clone https://github.com/SCUT-CCNL/LCQS.git && \
	cd LCQS && \
	make
ENV PATH="${PATH}:/software/LCQS/LCQS"

###################################################################################################################################################################################################
# FCLQC	
###################################################################################################################################################################################################

RUN mkdir -p /software/FCLQC && cd /software/FCLQC
RUN cd /software/FCLQC && git clone https://github.com/Minhyeok01/FCLQC.git
RUN cd /software/FCLQC/FCLQC && cargo build --release




