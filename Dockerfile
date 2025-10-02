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
    wget curl gnupg software-properties-common pkg-config nasm \
    && rm -rf /var/lib/apt/lists/*

# Install Intel oneAPI compilers
RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
    | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
    | tee /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update && apt-get install -y \
    intel-oneapi-compiler-dpcpp-cpp \
    intel-oneapi-compiler-fortran \
    && rm -rf /var/lib/apt/lists/*

# Environment setup
ENV INTEL_ONEAPI_ROOT=/opt/intel/oneapi
RUN echo "source /opt/intel/oneapi/setvars.sh" >> ~/.bashrc

# Entry point script
RUN echo '#!/bin/bash\n\
    source /opt/intel/oneapi/setvars.sh\n\
    export IPCL_DIR=/opt/intel/ipcl\n\
    export CMAKE_PREFIX_PATH=/opt/intel/ipcl/lib/cmake/ipcl-2.0.0:$CMAKE_PREFIX_PATH\n\
    export LD_LIBRARY_PATH=/opt/intel/ipcl/lib:$LD_LIBRARY_PATH\n\
    exec "$@"' > /entrypoint.sh && chmod +x /entrypoint.sh

RUN git clone https://github.com/emp-toolkit/emp-tool.git --branch master && \
    cd emp-tool && \
    cmake . && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf emp-tool  # Optional cleanup

RUN git clone https://github.com/flintlib/flint.git && cd flint && \
    ./bootstrap.sh && \
    ./configure --with-gmp=/usr/local --with-mpfr=/usr/local && \
    make -j$(nproc) && \
    make check && \
    make install && \
    cd .. && \
    rm -rf flint  # Optional cleanup

RUN GIT_TERMINAL_PROMPT=0 git clone https://github.com/intel/pailliercryptolib.git ipcl_src && \
    cd ipcl_src && \
    git checkout v2.0.0 && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/intel/ipcl -DIPCL_DISABLE_TEST=ON && \
    make -j$(nproc) install && \
    cd ../.. && \
    rm -rf ipcl_src

# Set working directory
WORKDIR /zebragram

# Copy the rest of the repo
COPY . .

# Set default command (optional)
ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]
