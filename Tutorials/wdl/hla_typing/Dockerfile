# Install samtools and bedtools on auxiliary image

FROM ubuntu:18.04 as auxiliary

# set the environment variables
ENV samtools_version 1.9
ENV bedtools_version 2.29.2

# run update and install necessary tools from package manager
RUN apt-get update -y && apt-get install -y \
    build-essential \
    cmake \
    zlib1g-dev \
    libhdf5-dev \
    libnss-sss \
    curl \
    autoconf \
    bzip2 \
    python3-dev \
    python3-pip \
    python3-biopython \
    pigz \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    unzip


WORKDIR /usr/bin/
RUN curl -SL https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    > samtools-${samtools_version}.tar.bz2 && \
    tar -xjvf samtools-${samtools_version}.tar.bz2 && \
    cd /usr/bin/samtools-${samtools_version} && \
    ./configure && \
    make && \
    make install

# install bedtools
WORKDIR /usr/bin
RUN curl -SL https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools-${bedtools_version}.tar.gz > bedtools-${bedtools_version}.tar.gz && \
    tar -xzvf bedtools-${bedtools_version}.tar.gz && \
    cd /usr/bin/bedtools2 && \
    make && \
    ln -s /usr/bin/bedtools2/bin/bedtools /usr/bin/bedtools


# Making final image
FROM ubuntu:18.04

COPY --from=auxiliary /usr/bin/samtools-1.9/samtools /usr/bin/samtools-1.9/samtools
COPY --from=auxiliary /usr/bin/bedtools2/bin/bedtools /usr/bin/bedtools
ENV PATH="${PATH}:/usr/bin/samtools-1.9"

# set the environment variables
ENV kallisto_version 0.44.0
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# run update and install necessary tools from package manager
RUN apt-get update -y && apt-get install -y \
    build-essential \
    cmake \
    zlib1g-dev \
    libhdf5-dev \
    libnss-sss \
    curl \
    autoconf \
    bzip2 \
    python3-dev \
    python3-pip \
    python3-biopython \
    pigz \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    unzip

# install python libraries
RUN pip3 install numpy && \
    pip3 install scipy && \
    pip3 install pandas

# install kallisto
RUN mkdir -p /usr/bin/kallisto \
    && curl -SL https://github.com/pachterlab/kallisto/archive/v${kallisto_version}.tar.gz \
    | tar -zxvC /usr/bin/kallisto && \
    mkdir -p /usr/bin/kallisto/kallisto-${kallisto_version}/build && \
    cd /usr/bin/kallisto/kallisto-${kallisto_version}/build && cmake .. && \
    cd /usr/bin/kallisto/kallisto-${kallisto_version}/ext/htslib && autoreconf && \
    cd /usr/bin/kallisto/kallisto-${kallisto_version}/build && make && \
    cd /usr/bin/kallisto/kallisto-${kallisto_version}/build && make install


# git lfs
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    apt-get install -y git-lfs && \
    git lfs install --system --skip-repo

# install arcasHLA
WORKDIR /home/
RUN git clone https://github.com/RabadanLab/arcasHLA.git arcasHLA-master

ENV PATH="${PATH}:/home/arcasHLA-master/"

RUN arcasHLA reference --update
