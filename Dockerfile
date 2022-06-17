FROM ubuntu:16.04

# Make bash the default login shell, i.e. supported
# by the conda init command
SHELL [ "/bin/bash", "--login", "-c" ]

# Install wget
RUN apt update && apt upgrade && apt install wget -y

# Create a non-root user
ARG username=zocalo
ARG uid=1000
ARG gid=100
ENV USER $username
ENV UID $uid
ENV GID $gid
ENV HOME /home/$USER
RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --gid $GID \
    --home $HOME \
    $USER

COPY environment.yaml requirements.txt /tmp/
RUN chown $UID:$GID /tmp/environment.yaml /tmp/requirements.txt
COPY docker-entrypoint.sh /usr/local/bin/
RUN chown $UID:$GID /usr/local/bin/docker-entrypoint.sh && \
    chmod u+x /usr/local/bin/docker-entrypoint.sh

USER $USER
# install miniconda as the non-root user
ENV MINICONDA_VERSION py39_4.12.0
ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh
# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH
# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

# create a project directory inside user home
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build the conda environment
ENV ENV_PREFIX $PROJECT_DIR/env
RUN conda update --name base --channel defaults conda --yes
RUN conda env create --prefix $ENV_PREFIX --file /tmp/environment.yaml --force
RUN conda clean --all --yes

ENTRYPOINT [ "/usr/local/bin/docker-entrypoint.sh" ]

CMD zocalo.service -s DLSPyMcaFitter
