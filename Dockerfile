FROM python:3.9-slim-bullseye 
LABEL maintainer="jsolislemus <jose.solislemus@kcl.ac.uk>"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends git libsm6 libxext6 libxrender-dev libgl1-mesa-glx ffmpeg && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN git clone git@bitbucket.org:aneic/meshtool.git /meshtool && \
    cd /meshtool && make && \
    ln -s /meshtool/meshtool /usr/local/bin/meshtool 


ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV && \
    /opt/venv/bin/python3 -m pip install --upgrade pip && \
    mkdir -p /code && mkdir -p /data

COPY . /code/
RUN pip install -r /code/requirements.txt 

WORKDIR /code/

# CMD [-h]
# ENTRYPOINT ["/opt/venv/bin/python3", "-u", "/code/docker/entrypoint.py"]

