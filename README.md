## XrayContourBinning

Code to perform approximate contour binning similar to the method described in [Sanders et al. 2006](). 

This code is designed to work with the combined background-subtracted and exposure-corrected adaptively-smoothed image produced by the XMMSAS task `adapt`. 

Unlike the original contour binning method described in [Sanders et al. 2006](), which defines regions based on a signal-to-noise constraint applied to the smoothed count image, we instead define regions based on an approximate total count constraint applied to the adaptively-smoothed combined count rate image.