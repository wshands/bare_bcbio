FROM bcbio/bcbio:latest

MAINTAINER Walt Shands jshands@ucsc.edu

ENV TEMPDIR=/work

WORKDIR $TEMPDIR

USER root

COPY UCSCbcbioTool.py /usr/local/bin
RUN chmod a+x /usr/local/bin/UCSCbcbioTool.py

# switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 ubuntu

# bcbio_nextgen creates tmp directories in /home/ubuntu so we need to create 
# /home/ubuntu in the root as owned by ubuntu before switching to user ubuntu
RUN mkdir /home/ubuntu
RUN chown ubuntu:ubuntu /home/ubuntu

RUN chmod a+rw /mnt/biodata

USER ubuntu

VOLUME /mnt/biodata/

