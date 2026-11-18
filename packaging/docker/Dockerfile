FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl gcc-14 g++-14 make binutils tar gzip && \
    rm -rf /var/lib/apt/lists/*

# Install Eclipse Temurin JDK 25 LTS & Rust 1.85+ cleanly in builder
RUN curl -fsSL "https://api.adoptium.net/v3/binary/latest/25/ga/linux/x64/jdk/hotspot/normal/eclipse" -o /tmp/jdk25.tar.gz && \
    mkdir -p /opt/jdk25 && tar -xzf /tmp/jdk25.tar.gz -C /opt/jdk25 --strip-components=1 && rm /tmp/jdk25.tar.gz
RUN curl -fsSL https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain 1.85.0

ENV PATH="/root/.cargo/bin:/opt/jdk25/bin:${PATH}"
ENV CC="gcc-14" CXX="g++-14" AR="ar" PREFIX="/opt/eimdl" SRC_DIR="/src"

WORKDIR /src
COPY . .
RUN bash packaging/bioconda/eimdl/build.sh

FROM ubuntu:24.04 AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends libstdc++6 ca-certificates && \
    rm -rf /var/lib/apt/lists/*
COPY --from=builder /opt/jdk25 /opt/jdk25
COPY --from=builder /opt/eimdl /usr/local
ENV PATH="/opt/jdk25/bin:/usr/local/bin:${PATH}"

LABEL org.opencontainers.image.title="EIMDL / HyperGraSS v3" \
      org.opencontainers.image.description="Evolutionary lossless genomic compressor (.hgs2/.hgs3) + compressed-domain analytics" \
      org.opencontainers.image.licenses="Apache-2.0"

ENTRYPOINT ["eimdl"]
CMD ["--version"]
