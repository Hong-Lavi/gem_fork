# Install ubuntu
FROM ubuntu:14.04
MAINTAINER Hyun Uk Kim, Jae Yong Ryu

ENV GMSM_PATH="/gmsm/"
ENV PATH="${PATH}:${GMSM_PATH}"

RUN apt-get -y update
RUN apt-get install -y python2.7
RUN apt-get install -y python2.7-dev
RUN apt-get install -y python-pip
RUN apt-get install -y python-tox
RUN apt-get install -y ncbi-blast+
RUN apt-get install -y sed

# Install major dependencies
COPY . /gmsm
WORKDIR /gmsm/
RUN pip install pip --upgrade
RUN pip install -r requirements.txt

VOLUME ["/input", "/output"]
WORKDIR /gmsm/
RUN chmod +x run_gmsm.py

ENTRYPOINT ["/bin/bash"]
