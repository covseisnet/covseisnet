.. covnet documentation master file, created by
   sphinx-quickstart on Mon Mar 26 10:49:33 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Home
====

This package is an open-source project dedicated to the analysis of array seismic data. The project includes several classical tools such as beamforming, cross-correlation and spectral analysis. Seismic data reading and processing is mostly based on the ObsPy library (https://github.com/obspy/obspy/wiki).

Contents
--------

.. toctree::
   :maxdepth: 2

   api
   license

Manual installation
-------------------

At the moment, the only way to get this package is to clone it from GitHub. In a terminal, you can create a directory of your choice (for instance a directory where you can store all packages manually installed), and execute the following lines:

.. code-block:: bash

    git clone https://github.com/leonard-seydoux/covnet.git .
    cd covnet
    pip install .

Once you installed the package, you can verify your install with executing the following command in a python shell (outside the package repository).

.. code-block:: python

    import covnet as cn

References
----------

The following publications have made use of this package, or of a development version of it. Please refer to the description provided therein in order to get in-depth descriptions of the tools defined in the package as well as an overview of possible applications.

- Soubestre, Jean, et al. Network‐Based Detection and Classification of Seismovolcanic Tremors: Example From the Klyuchevskoy Volcanic Group in Kamchatka. *Journal of Geophysical Research: Solid Earth* 123.1 (2018): 564-582.

- Seydoux, Léonard, Julien de Rosny, and Nikolai M. Shapiro. "Pre-processing ambient noise cross-correlations with equalizing the covariance matrix eigenspectrum." Geophysical Journal International 210.3 (2017): 1432-1449.

- Seydoux, L., et al. "Spatial coherence of the seismic wavefield continuously recorded by the USArray." Geophysical Research Letters 43.18 (2016): 9644-9652.

- Seydoux, Léonard, et al. "Detecting seismic activity with a covariance matrix analysis of data recorded on seismic arrays." Geophysical Journal International 204.3 (2016): 1430-1442.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Licence
-------

This package is released under the GNU General Public License v3.0 (see :ref:`license`)

Todo
----

- Add the installation process from ``conda``.
- Add recent papers and doi links.
- Add a "How to cite this package" section.
- Discuss the copyright and licence formatting.
- Review the ``conf.py`` file meta data.

