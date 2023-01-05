FROM ubuntu:22.04

RUN apt-get update && \
    DEBAIN_FRONTEND=noninteractive apt-get install -qy xinetd python3 python3-pip && \
    pip3 install pycryptodome

RUN useradd -m chal && \
    chown -R root:root /home/chal && \
    chmod -R 755 /home/chal


CMD ["/usr/sbin/xinetd", "-dontfork"]
