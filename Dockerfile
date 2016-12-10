FROM ubuntu:16.04

RUN apt-get update -y 
RUN apt-get install -y build-essential
RUN apt-get install -y gdb 
RUN apt-get install -y git 
RUN apt-get install -y vim 
RUN apt-get install -y libcr-dev mpich
RUN apt-get install -y libblas-dev liblapack-dev 
