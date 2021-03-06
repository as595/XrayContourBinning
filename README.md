# XrayContourBinning

Code to perform approximate contour binning similar to the method described in [Sanders et al. 2006](https://arxiv.org/abs/astro-ph/0606528). 

This code is designed to work with the combined background-subtracted and exposure-corrected adaptively-smoothed image produced by the XMMSAS task `adapt`. 

Unlike the original contour binning method described in [Sanders et al. 2006](https://arxiv.org/abs/astro-ph/0606528), which defines regions based on a signal-to-noise constraint applied to the smoothed count image, we instead define regions based on an approximate total count constraint applied to the adaptively-smoothed combined count rate image.

## Run the code

The input parameters for each run are contained in the configuration files located in the [configs](./configs) directory. To run a particular experiment use:

```python
python main.py --config configs/config_a1314.cfg
```
An overview of the configuration file format can be found [here](./configs/README.md).

## Example notebook

A [jupyter notebook](./ContourBins.ipynb) is also provided to illustrate the steps in the contour binning process in a more accessible way.

## Environment

To (optionally) create region files for spectral fitting in XMM-SAS the code requires the XMM-SAS software and its dependencies to be available. The easiest way to achieve this is to use [this Docker image for XMM-Newton SAS](https://gitlab.astro.unige.ch/ferrigno/sas-docker/-/tree/master). If you build the docker image directly from the [Dockerfile](https://gitlab.astro.unige.ch/ferrigno/sas-docker/-/blob/master/Dockerfile) using the updated [requirements.txt](./container/requirements.txt) file from this repository then there is no need to install additional Python libraries for running the PyContourBin functions.
