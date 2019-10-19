ARG paraview=5.6.1
FROM topologytoolkit/ttk-base:${paraview}

MAINTAINER Christoph Garth <garth@cs.uni-kl.de>

USER root

# install TTK
ARG ttk=v0.9.8
ENV TTK_VERSION=${ttk}

COPY ttk.sh /root/
RUN  /root/build /root/ttk.sh

USER paraview

