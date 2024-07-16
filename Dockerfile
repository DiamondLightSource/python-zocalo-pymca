FROM debian:12-slim

RUN apt-get update && apt-get install -y curl bzip2
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
COPY requirements.txt .
RUN /bin/micromamba create -yp /opt/pymca -c conda-forge gnuplot python --file requirements.txt
COPY . /tmp/pymca
RUN /opt/pymca/bin/python -mpip install /tmp/pymca

#FROM gcr.io/distroless/python3-debian12:debug-nonroot
#COPY --from=0 /opt/pymca /opt/pymca
CMD ["/opt/pymca/bin/zocalo.service", "-s", "DLSPyMcaFitter"]  