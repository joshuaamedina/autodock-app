FROM tacc/tacc-ubuntu18-mvapich2.3-ib:latest

# Installing mpi4py
RUN apt-get update && apt-get upgrade -y && apt-get install -y python3-pip
RUN apt-get install wget -y
RUN curl -O https://bootstrap.pypa.io/pip/3.6/get-pip.py && python3 get-pip.py
RUN pip3 install mpi4py && pip3 install numpy==1.19.5
RUN pip install -U numpy vina

RUN wget -O ADFR.tar  https://ccsb.scripps.edu/adfr/download/1038/
RUN tar -xf ADFR.tar
RUN cd ADFRsuite_x86_64Linux_1.0 && ./install.sh -c 0
RUN cd
ENV PATH="$PATH:/ADFRsuite_x86_64Linux_1.0/bin"
