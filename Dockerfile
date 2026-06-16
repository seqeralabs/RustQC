# ---- Build stage ----
FROM rust:1-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    libfontconfig1-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

ENV CXX=g++

ARG GIT_SHORT_HASH=unknown
ARG CPU_TARGET=""

WORKDIR /build
COPY Cargo.toml Cargo.lock build.rs ./
COPY cpp/ cpp/
COPY src/ src/

RUN HOST_TRIPLE=$(rustc -vV | awk '/^host:/ {print $2}') && \
    GIT_SHORT_HASH="${GIT_SHORT_HASH}" \
    cargo build --release --target "$HOST_TRIPLE" \
        ${CPU_TARGET:+--config "target.'$HOST_TRIPLE'.rustflags=['-C', 'target-cpu=$CPU_TARGET']"} \
    && strip "target/$HOST_TRIPLE/release/rustqc" \
    && cp "target/$HOST_TRIPLE/release/rustqc" /rustqc

# ---- Runtime stage ----
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    fontconfig \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /rustqc /usr/local/bin/rustqc

CMD ["rustqc"]
