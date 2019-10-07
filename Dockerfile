# http://www.science.smith.edu/dftwiki/index.php/Tutorial:_Docker_Anaconda_Python_--_4#Exercise_5:_Creating_a_Container_with_Anaconda
# https://github.com/ContinuumIO/docker-images/issues/89 for activation of conda env in docker
# adapted from https://hub.docker.com/r/continuumio/anaconda3/dockerfile
# adapted from https://gist.githubusercontent.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a/raw/b0a841c4be8873ba73faba217e798d3d9e207823/Dockerfile_v5.5.sh
FROM debian:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV PATH=/opt/conda/bin:$PATH

### Install MToolBox
# update conda
RUN conda update conda

# fetch repo
RUN git clone https://github.com/mitoNGS/MToolBox_snakemake.git

# install mtoolbox conda env
WORKDIR MToolBox_snakemake
RUN ["conda", "env", "create", "-n", "mtoolbox", "-f", "envs/mtoolbox.yaml"]

RUN apt-get install -y curl grep sed dpkg sudo
RUN sudo apt-get -y install make gcc g++
RUN sudo apt-get -y install zlib1g-dev
ENV CONFIGFLAGS $CONFIGFLAGS" -lz"
# install bamUtils
#RUN ["/bin/bash", "-c", "conda", "activate", "mtoolbox"]
RUN git clone https://github.com/statgen/bamUtil.git
WORKDIR bamUtil
RUN make cloneLib
RUN make
RUN make install INSTALLDIR=$(dirname $(which python))
#RUN conda deactivate
WORKDIR ..
WORKDIR ..

RUN echo "export PATH=/opt/conda/envs/mtoolbox/bin/:/MToolBox_snakemake:/MToolBox_snakemake/scripts:$PATH" >> /root/.bashrc
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: MToolBox :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /root/.bashrc

# Run Jupytewr notebook as Docker main process
#CMD ["jupyter", "notebook", "--allow-root", "--notebook-dir=/home/ubuntu/notebooks", "--ip='*'", "--port=8888", "--no-browser"]
# Activate mtoolbox conda environment as Docker main process
# SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash", "-l"]

