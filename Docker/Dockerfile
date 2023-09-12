# nf-perturbseq image

################################################################################

FROM romanhaa/nf-pgp-perturbseq:0.1

# STAR
RUN rm -r tools && \
    wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
    tar -xzf 2.7.9a.tar.gz && \
    cd STAR-2.7.9a/source && \
    make
ENV PATH="/STAR-2.7.9a/source/:$PATH"

RUN conda install pigz

################################################################################