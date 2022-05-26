FROM python:3.10

WORKDIR /home/zocalo/pymca

RUN python -m pip install --upgrade pip
RUN python -m pip install --upgrade --no-cache-dir cython numpy setuptools wheel
COPY requirements.txt .
RUN python -m pip install --upgrade --no-cache-dir -r requirements.txt

COPY . .

RUN python -m pip install --no-cache-dir .
