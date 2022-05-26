FROM python:3.10

RUN useradd --create-home --shell /bin/bash zocalo
USER zocalo

WORKDIR /home/zocalo/pymca

RUN python -m pip install --user --upgrade pip
RUN python -m pip install --user --upgrade --no-cache-dir cython numpy setuptools wheel
COPY requirements.txt .
RUN python -m pip install --user --upgrade --no-cache-dir -r requirements.txt

COPY . .

RUN python -m pip install --user --no-cache-dir .
