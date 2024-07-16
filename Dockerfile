FROM debian:12-slim

RUN apt-get update && apt-get install -y curl bzip2 gnuplot
# install miniconda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
# make non-activate conda commands available
COPY requirements.txt .
RUN /opt/miniconda/bin/conda create -yp /opt/pymca -c conda-forge --file requirements.txt python
COPY . /tmp/pymca
RUN /opt/pymca/bin/python -mpip install /tmp/pymca --no-cache-dir --no-dependencies

#FROM gcr.io/distroless/python3-debian12:debug-nonroot
#COPY --from=0 /opt/pymca /opt/pymca
CMD ["/opt/pymca/bin/zocalo.service", "-s", "DLSPyMcaFitter"]    