.. _guide:

User Guide
============


The `covseisnet` package provides tools for array signal processing, with a focus on data from seismic networks. The central analyzed mathematical construction is the array covariance matrix (sometimes called network covariance matrix). The core signal detection algorithms are based on the analysis of the eigenvalues of this matrix. Eigenvector decomposition provides a basis for a blind source separation. In addition, the covariance matrix can be used as input for classical array processing tools such as beamforming and inter-station cross-correlations. Covseisnet objects are mostly inherited from obspy (seismology-oriented python package) and `numpy`.

The code repository is hosted on GitHub at https://github.com/covseisnet/covseisnet

Guides
------

.. toctree::
    :maxdepth: 2

    roadmap
    guide_stream
    guide_covariance
    guide_traveltimes
    guide_correlation
    guide_beamforming