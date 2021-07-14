# Environment

The Dockerfile here is a copy of the [XMM-Newton SAS Dockerfile](https://gitlab.astro.unige.ch/ferrigno/sas-docker/-/tree/master), authored by Carlo Ferrigno. If you make use of it please cite Carlo's work.

The [requirements.txt](./requirements.txt) file in this repository is slightly updated in order to enable the PyContourBin functionality. 

## Building the image

To build the docker image you should first get the source tar balls for ds9 and xmm-sas

```
wget -q https://ds9.si.edu/download/ubuntu20/ds9.ubuntu20.8.2.1.tar.gz
wget -q http://sasdev-xmm.esac.esa.int/pub/sas/19.1.0/Linux/Ubuntu18.04/sas_19.1.0-Ubuntu18.04.tgz
```

and then build the image using

```
docker build -t contour .
```

(or a tag of your choice).
