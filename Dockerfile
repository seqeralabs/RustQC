# ---- Build stage ----
FROM rust:1-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    pkg-config \
    clang \
    && rm -rf /var/lib/apt/lists/*

ARG GIT_SHORT_HASH=unknown
ARG CPU_TARGET=""

WORKDIR /build
COPY Cargo.toml Cargo.lock build.rs ./
COPY cpp/ cpp/
COPY src/ src/

RUN RUSTFLAGS="${CPU_TARGET:+-C target-cpu=$CPU_TARGET}" \
    GIT_SHORT_HASH="${GIT_SHORT_HASH}" \
    cargo build --release && strip target/release/rustqc

# ---- Runtime stage ----
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    fontconfig \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/rustqc /usr/local/bin/rustqc

CMD ["rustqc"]
