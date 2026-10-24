# Native arm64 image of LLNL Spheral at develop 116c71f (PR #520 merge, 2026-09-28).
#
# Follows LLNL's own root Dockerfile at that commit (same apt packages, same tpl-manager /
# host-config-build / make sequence) with these differences:
#   - built for linux/arm64 (LLNL's ghcr.io image is x86_64_v4 and dies under Rosetta);
#   - the source is cloned by full SHA with submodules, not COPY'd from a local tree;
#   - Spheral is installed at /opt/spheral (a real install prefix, not a symlink);
#   - make -j $JCXX Spheral_CXX / make -j $JPY with JCXX=16, JPY=4;
#   - the upstream ATS suite is not run in the image (Task 3 runs the parts M0 needs);
#   - TPLs and Spheral are separate RUN layers, so a late failure does not rebuild the TPLs.
# No arm64-specific fix was needed (first build, 2026-10-05: 1 h 51 min). Build with
# spheral_frag/container/build.sh; pin.json records what the image contains.

FROM ubuntu:24.04

ARG SPEC=gcc
ARG HOST_CONFIG=docker-$SPEC
ARG SPHERAL_SHA=116c71f2ea559d569456ba694194e5a6de3c75f7
ARG JCXX=16
ARG JPY=4

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
        ca-certificates \
        netbase \
        iproute2 \
        build-essential \
        git \
        gfortran \
        autotools-dev \
        autoconf \
        sqlite3 \
        pkg-config \
        uuid \
        gettext \
        cmake \
        openmpi-bin \
        libopenmpi-dev \
        libncurses-dev \
        libgdbm-dev \
        libffi-dev \
        libssl-dev \
        libexpat-dev \
        libreadline-dev \
        libbz2-dev \
        locales \
        python3 \
        python3-dev \
        python3-venv \
        python3-pip \
        python3-tk \
        unzip \
        libtool \
        wget \
        curl \
        libcurl4-openssl-dev \
        tk-dev \
        iputils-ping && \
    rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US.UTF-8

# Source by full SHA, with submodules (blt, PYB11Generator, PolyClipper, ATS).
WORKDIR /home/spheral/workspace
RUN git clone https://github.com/LLNL/spheral.git . && \
    git checkout $SPHERAL_SHA && \
    git submodule update --init --recursive && \
    test "$(git rev-parse HEAD)" = "$SPHERAL_SHA"

# Third-party libraries through Spack (the multi-hour layer). Spack lives in /home/spack,
# installs in /home/installs. Spack detects the VM's Apple CPU (CPU part 0x000) as generic aarch64.
RUN python3 scripts/devtools/tpl-manager.py --spec spheral%$SPEC --tpl-dir /home

# Configure Spheral against those TPLs, installing straight into /opt/spheral.
RUN mv *.cmake $HOST_CONFIG.cmake && \
    python3 scripts/devtools/host-config-build.py --host-config $HOST_CONFIG.cmake --install-dir /opt/spheral

WORKDIR /home/spheral/workspace/build_$HOST_CONFIG/build
RUN make python_build_env
RUN make python_runtime_env
RUN make -j $JCXX Spheral_CXX
RUN make -j $JPY
RUN make install

ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR /opt/spheral
