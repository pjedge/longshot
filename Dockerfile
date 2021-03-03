FROM ubuntu:18.04 AS install-longshot

## Install dependencies
#### you can use either wget or git+ca-certificates, no need to get both
#### musl-tools used to install x86_64-unknown-linux-musl so your binary has no external dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    cargo \
    zlib1g-dev \
    xz-utils \
    libclang-dev \
    clang \
    build-essential \
    curl \
    wget \
    ca-certificates \
    git \
    musl-tools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## Install RUST
#### x86_64-unknown-linux-musl is used so the binary won't have external dependencies
RUN curl https://sh.rustup.rs -o install-rustup.sh
RUN chmod +x install-rustup.sh
RUN ./install-rustup.sh -y --target x86_64-unknown-linux-musl
ENV PATH=/root/.cargo/bin:$PATH

## Download longshot
### Fixed version to 0.4.3
RUN wget https://github.com/pjedge/longshot/archive/v0.4.3.tar.gz
RUN tar -xzvf v0.4.3.tar.gz
WORKDIR /longshot-0.4.3

## Clone longshot repository
#### This is an alternative way to the download option. If you use this, comment the download section out.
#### This has no fixed version, and will just get the latest on master
# RUN git clone https://github.com/pjedge/longshot
# WORKDIR /longshot

## Install longshot
RUN cargo build --target x86_64-unknown-linux-musl --release
RUN cargo install --target x86_64-unknown-linux-musl --path .

## Copy the binary to a minimal docker image
FROM alpine:3.12.4
# RUN apk add --no-cache bash # My workflow required a bash entrypoint. Commented it out for a more general use.
COPY --from=install-longshot /root/.cargo/bin/longshot /bin/longshot
