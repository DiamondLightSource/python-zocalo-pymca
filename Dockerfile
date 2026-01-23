FROM debian:12-slim

RUN apt-get update && apt-get install -y curl bzip2
# install miniconda
RUN curl -LO "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" && \
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniforge3-Linux-x86_64.sh
# make non-activate conda commands available
COPY requirements.txt .
RUN /opt/miniconda/bin/conda create -yp /opt/pymca -c conda-forge --file requirements.txt python
COPY . /tmp/pymca
RUN /opt/pymca/bin/python -mpip install /tmp/pymca --no-cache-dir --no-dependencies
ENV HOME /home

#FROM gcr.io/distroless/python3-debian12:debug-nonroot
#COPY --from=0 /opt/pymca /opt/pymca
CMD ["/opt/pymca/bin/zocalo.service", "-s", "DLSPyMcaFitter"]  