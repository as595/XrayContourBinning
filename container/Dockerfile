FROM carloferrigno/heasoft:latest

USER root
#DS9
COPY ds9.ubuntu20.8.2.1.tar.gz /opt
#    wget -q https://ds9.si.edu/download/ubuntu20/ds9.ubuntu20.8.2.1.tar.gz && \
RUN cd /opt && \
    tar xvfz ds9.ubuntu20.8.2.1.tar.gz && \
    chmod a+x ds9 && \
    mv ds9 /usr/local/bin && \
    rm -f ds9.ubuntu20.8.2.1.tar.gz

#default shell to bash
SHELL ["/bin/bash", "-c"]

# SAS-compatible Conda environemnt creation

#gets SAS
ENV SAS_PERL=/usr/bin/perl
COPY sas_19.1.0-Ubuntu18.04.tgz /opt/
#    wget -q http://sasdev-xmm.esac.esa.int/pub/sas/19.0.0/Linux/Ubuntu18.04/sas_19.0.0-Ubuntu18.04.tgz && \
RUN cd /opt && \
    tar xfz sas_19.1.0-Ubuntu18.04.tgz && \
    ./install.sh && \
    rm sas_19.1.0-Ubuntu18.04.tgz

RUN ln -s /opt/xmmsas_20210317_1624  /opt/xmmsas

RUN apt-get update \
 && apt-get -y upgrade \
 && apt-get -y install \
        rsync

# CCF
RUN rsync -v -a --delete --delete-after --force --include=*.CCF --exclude=*/ sasdev-xmm.esac.esa.int::XMM_VALID_CCF /opt/ccf/

#Adds init file
ADD init_sas.sh /usr/local/bin/init_sas.sh

RUN pip3 install pip --upgrade
ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt --upgrade && \
    rm /requirements.txt

RUN apt-get update && \
    apt-get -y install software-properties-common

RUN add-apt-repository ppa:linuxuprising/libpng12 && \
    apt-get update && \
    apt-get -y install libpng12-0

USER heasoft
WORKDIR /home/heasoft

