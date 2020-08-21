FROM ubuntu:16.04

LABEL base.image="variant_QC:latest"
LABEL tags="Variant QC"

# Replace sh with bash
RUN cd /bin && rm sh && ln -s bash sh

ARG DEBIAN_FRONTEND=noninteractive

# Configure locale and timezone
RUN apt-get update && apt-get install -y locales && \
    echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure -f noninteractive locales && \
    /usr/sbin/update-locale LANG=en_US.UTF-8; \
    apt-get clean && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

# Install miniconda3 at /opt/conda
RUN apt-get update \
    && apt-get -y install --no-install-recommends curl bzip2 ca-certificates \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /opt/conda \
    && rm -rf /tmp/miniconda.sh \
    && apt-get -y remove curl bzip2 \
    && apt-get -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && /opt/conda/bin/conda clean --all --yes

ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8" \
    PATH="/opt/conda/bin:${PATH}"                               \
    SHELL="/bin/bash"

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

## Install git
RUN apt-get update \
    && apt-get install -y --no-install-recommends git apt-transport-https gnupg2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN conda install -y python=3.6 pindel varscan gatk4 samtools snakemake pandas bcftools jupyterlab && \
    conda clean -y --all

USER root
WORKDIR /home/

ADD ./docker-entrypoint.sh /home/
ADD ./.bashrc /home/

#RUN source /home/.bashrc
#RUN jupyter lab --allow-root --no-browser --ip 0.0.0.0 --port 8888 &
RUN chmod +x ./docker-entrypoint.sh
#CMD ./docker-entrypoint.sh
#ENTRYPOINT ["/bin/bash"]
ENTRYPOINT ["/home/docker-entrypoint.sh"]
#RUN git clone https://github.com/ksanao/variant_QC

