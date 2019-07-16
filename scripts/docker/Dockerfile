FROM chgarth/ttk-base:latest

MAINTAINER Christoph Garth <garth@cs.uni-kl.de>

USER root

# install TTK
ARG ttk=v0.9.7
ENV TTK_VERSION=${ttk}

COPY ttk.sh /root/
RUN  /root/build /root/ttk.sh

USER paraview

