FROM ubuntu:16.04

MAINTAINER GuangyuWang <guangywang@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# Configure environment
ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
ENV SHELL /bin/bash
ENV CONDA_VER "4.1.11"

# Create conda dir
RUN mkdir -p $CONDA_DIR

# System dependencies
RUN apt-get clean all && \
	apt-get update
RUN apt-get install -y build-essential curl wget git pkg-config imagemagick gcc libgtk2.0-dev && \
	apt-get clean && \
	apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Miniconda3
RUN cd /tmp && \
    mkdir -p $CONDA_DIR && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-$CONDA_VER-Linux-x86_64.sh && \
    # echo "9ea57c0fdf481acf89d816184f969b04bc44dea27b258c4e86b1e3a25ff26aa0 *Miniconda3-3.19.0-Linux-x86_64.sh" | sha256sum -c - && \
    /bin/bash Miniconda3-$CONDA_VER-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda3-$CONDA_VER-Linux-x86_64.sh 


RUN $CONDA_DIR/bin/conda update --yes conda && \
###config is new added!!!
$CONDA_DIR/bin/conda config --add channels conda-forge && \
$CONDA_DIR/bin/conda config --add channels bioconda && \
$CONDA_DIR/bin/conda config --add channels r

RUN $CONDA_DIR/bin/conda install --yes \
	nomkl \
	glib \
	openblas \
	numpy \
	scipy \
	pillow \
	matplotlib \
	cython \
	pandas \
	gensim \
	protobuf \
	libgcc \
	imagehash \
	r-base=3.4* r-bcp r-cluster r-mass r-essentials r-spatstat libiconv bioconductor-hicrep r-changepoint tzlocal imagehash && \
	$CONDA_DIR/bin/conda install -c conda-forge r-spatstat libiconv && \
	$CONDA_DIR/bin/conda install r-optparse && \
	$CONDA_DIR/bin/conda install gtk2 && \
	$CONDA_DIR/bin/conda install -y -c menpo opencv3 && \
	python -c "import cv2;print(cv2.__version__)" && \
	$CONDA_DIR/bin/conda  clean -tipsy

 
RUN wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu55_55.1-7_amd64.deb
RUN dpkg -i libicu55_55.1-7_amd64.deb
RUN apt-get -f install
RUN rm libicu55_55.1-7_amd64.deb && \
    apt-get clean

WORKDIR /bin/
COPY /src/ /bin/
