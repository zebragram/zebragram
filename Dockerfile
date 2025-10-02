FROM ubuntu:22.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    ninja-build \
    python3 \
    git \
    libboost-stacktrace-dev \
    libssl-dev \
    libomp-dev \
    libgmp-dev \
    libmpfr-dev \
    autoconf \
    libtool-bin \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/emp-toolkit/emp-tool.git --branch master && \
    cd emp-tool && \
    cmake . && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf emp-tool  # Optional cleanup

RUN git clone https://github.com/flintlib/flint.git && cd flint && \
./bootstrap.sh && \
./configure && \
make && \
make install && \
make examples && \
cd .. && \ 
rm -rf flint  # Optional cleanup

# Set working directory
WORKDIR /picogram

# Copy repo and submodules
COPY . .

# Initialize and update git submodules
RUN git submodule update --init --recursive

# Set default command (optional)
CMD ["bash"]
