# http://www.science.smith.edu/dftwiki/index.php/Tutorial:_Docker_Anaconda_Python_--_4#Exercise_5:_Creating_a_Container_with_Anaconda
# https://github.com/ContinuumIO/docker-images/issues/89 for activation of conda env in docker
# adapted from https://hub.docker.com/r/continuumio/anaconda3/dockerfile
# adapted from https://gist.githubusercontent.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a/raw/b0a841c4be8873ba73faba217e798d3d9e207823/Dockerfile_v5.5.sh
FROM debian:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH
ENV CONFIGFLAGS $CONFIGFLAGS" -lz"

# Identify the maintainer of an image
LABEL maintainer="dome.simone@gmail.com"

# upgrade OS
RUN apt-get update --fix-missing
RUN apt-get install -y locales
RUN apt-get install -y --no-install-recommends \
    dpkg build-essential ca-certificates \
    wget bzip2 \
    git sudo make zlib1g-dev && \
	# install miniconda and setup env
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    conda update conda && \
	# install MToolBox dependencies
    git clone --branch devel_docker https://github.com/mitoNGS/MToolBox_snakemake.git && cd MToolBox_snakemake && \
    conda env create -n mtoolbox -f envs/mtoolbox.yaml && conda clean -a -y
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate mtoolbox" >> ~/.bashrc
ENV PATH="/MToolBox_snakemake:/MToolBox_snakemake/scripts:$PATH"
#RUN echo "export PATH=/opt/conda/envs/mtoolbox/bin/:/MToolBox_snakemake:/MToolBox_snakemake/scripts:$PATH" >> ~/.bashrc
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: MToolBox :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> ~/.bashrc

CMD ["/bin/bash", "-l"]
