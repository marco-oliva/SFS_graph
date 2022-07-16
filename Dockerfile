# Get the base Ubuntu image from Docker Hub
FROM ubuntu:20.04

# Update apps on the base image
RUN apt-get -y update

# Install GCC and dependencies
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get install -y wget bzip2 autoconf automake make cmake gcc g++ perl zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev libncurses5-dev git


# Start building
COPY . /usr/src/SFS_graph

# Get Themisto
WORKDIR /usr/src/
RUN git clone https://github.com/algbio/themisto.git \
    && cd themisto \
    && git submodule init \
    && git submodule update \
    && cd build \
    && cmake .. -DMAX_KMER_LENGTH=31 \
    && make -j \

# Get MONI
WORKDIR /usr/src/
RUN git clone https://github.com/marco-oliva/moni.git \
    && cd moni \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_INSTALL_PERFIX=/sfs/bin .. \
    && make -j \
    && make install

# Get binaries
WORKDIR /sfs/bin
RUN cp /usr/src/themisto/build/bin/themisto .
RUN cp /usr/src/SFS_graph/SFS_graph.py .
RUN chmod +x SFS_graph.py

ENV PATH /sfs/bin:$PATH