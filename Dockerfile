FROM python:3.9-slim-bullseye 
LABEL maintainer="jsolislemus <jose.solislemus@kcl.ac.uk>"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends git libsm6 libxext6 libxrender-dev libgl1-mesa-glx ffmpeg build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://bitbucket.org/aneic/meshtool.git && \
    cd /meshtool && make && \
    ln -s /meshtool/meshtool /usr/local/bin/meshtool 

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV && \
    /opt/venv/bin/python3 -m pip install --upgrade pip && \
    mkdir -p /code && mkdir -p /data
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
ENV PYTHONPATH="${PYTHONPATH}:/code"

COPY ./requirements.txt /requirements.txt
RUN /opt/venv/bin/python3 -m pip install -r /requirements.txt 

COPY . /code/
WORKDIR /code/

CMD ["-h"]
ENTRYPOINT ["/opt/venv/bin/python3", "-u", "/code/docker/entrypoint_script.py"]

