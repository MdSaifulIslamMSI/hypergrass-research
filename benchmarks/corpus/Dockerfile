FROM --platform=linux/amd64 eclipse-temurin@sha256:56f10ddc224d7b900698422ffe82432ed0e178f606bd66859ed1d179825dfbe3

RUN apk add --no-cache python3 build-base gzip xz bash

ENV JARVIS3_CFLAGS="-Wall -Ofast -DPROGRESS" \
    GECO3_CFLAGS="-Wall -O3"

COPY build_image.py sources.lock.json /opt/benchmark/
COPY vendor/jarvis3.tar.gz vendor/geco3.tar.gz /opt/benchmark/vendor/

# Keep complete source trees, including upstream licenses, alongside the lock.
# Only package installation above needs network access.
RUN --network=none python3 /opt/benchmark/build_image.py --verify-sources /opt/benchmark \
    && mkdir -p /opt/sources/jarvis3 /opt/sources/geco3 \
    && tar -xzf /opt/benchmark/vendor/jarvis3.tar.gz --strip-components=1 -C /opt/sources/jarvis3 \
    && tar -xzf /opt/benchmark/vendor/geco3.tar.gz --strip-components=1 -C /opt/sources/geco3 \
    && make -C /opt/sources/jarvis3/src CC=gcc CFLAGS="$JARVIS3_CFLAGS" \
    && make -C /opt/sources/geco3/src CC=gcc CFLAGS="$GECO3_CFLAGS" \
    && install -m 0755 /opt/sources/jarvis3/src/JARVIS3 /usr/local/bin/JARVIS3 \
    && install -m 0755 /opt/sources/geco3/src/GeCo3 /usr/local/bin/GeCo3 \
    && install -m 0755 /opt/sources/geco3/src/GeDe3 /usr/local/bin/GeDe3 \
    && python3 /opt/benchmark/build_image.py --record-environment /opt/benchmark-environment.json

WORKDIR /work
CMD ["bash"]
