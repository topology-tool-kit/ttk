ARG paraview=5.8.1

FROM topologytoolkit/ttk-base:${paraview}
LABEL maintainer="Christoph Garth <garth@cs.uni-kl.de>"

# install TTK
ARG ttk=dev
ENV TTK_VERSION=${ttk}

COPY ttk.sh /root/
RUN  /root/helper.sh /root/ttk.sh
