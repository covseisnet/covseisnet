# CovSeisNet

Analaysing seismic data with array processing tools.

[![PyPI Version](https://img.shields.io/pypi/v/covseisnet.svg)](https://test.pypi.org/project/covseisnet)
[![Conda Version](https://img.shields.io/conda/v/francistong/covseisnet)](https://anaconda.org/francistong/covseisnet)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/conda/l/francistong/covseisnet)](https://www.gnu.org/licenses/lgpl.html)
[![Python Version](https://img.shields.io/pypi/pyversions/covseisnet)](https://test.pypi.org/project/covseisnet)
[![Platform](https://img.shields.io/conda/pn/francistong/covseisnet)](https://anaconda.org/francistong/covseisnet)
[![Pipeline Build](https://gricad-gitlab.univ-grenoble-alpes.fr/covseisnet/covseisnet/badges/develop/pipeline.svg)]()
[![Coverage](https://img.shields.io/codecov/c/gh/covseisnet/covseisnet?token=17eccecf8a2846e497a9c35f10e396de)]()

## Get and pre-process seismic data

### Read

The package provides a top-level function `read` (module: `data.py`) based on the Obspy's read function. It returns a object of class `stream` with all the useful tools provided by Obspy. It can read almost all file formats (`.sac`, `.mseed`, ...). 

```
import arrayprocessing as ap

stream = ap.read('/path/to/traces/*.sac')
```

### Pre-processing methods

In addition to the original methods attached to Obspy's `stream` class, some new methods are defined. 

#### Trim the data without using `UTCDateTime` in the main code

This method simply allows not to import `obspy.UTCDateTime` in the main script in order to trim the data with the `trim` method. Additional kwargs to `starttime` and `endtime` are passed to the `trim` method of `obspy.core.stream.Stream` class.

```
import arrayprocessing as ap

stream = ap.read('/path/to/traces/*.sac')
stream.cut(starttime='2017-01-01 01:34', endtime='2017-01-01 01:34')
```

#### Binarize the traces in the temporal domain

Divides each trace by its envelope. In order not to divide by 0, a regularization kwarg `epsilon` is define by default to `1e-10`.

```
import arrayprocessing as ap

stream = ap.read('/path/to/traces/*.sac')
stream.binarize()
```

#### Stationarize in the temporal domain

Divides each trace by their smooth envelopes. The smoothing is calculated with the Savitzky-Golay algorithm. The `window` kwargs controls the length of smoothing window, and the `order` controls the order of the method. Defaults are `window=11, order=1`. As with the binarization methods, a default `epsilon` is defined to `1e-10`.

```
import arrayprocessing as ap

stream = ap.read('/path/to/traces/*.sac')
stream.binarize()
```
#### Whiten each trace in the frequency domain

Divides each trace spectrum with its amplitude. if the `method` kwarg is set to `onebit`, no smoothing is applied to the amplitude. Otherwise, the amplitude is smoothed with the Savitzky-Golay algorithm, with a window of length `smooth` and a order 1.

```
import arrayprocessing as ap

stream = ap.read('/path/to/traces/*.sac')
stream.whiten(method='smooth, smooth=11)
```

### Show data

You can use the methods provided by obspy in order to show the data. Also, a `show` method is defined, currently under tests.

## Calculate spectrograms

The spectrograms are calculated with the short-time Fourier transform function provided by Scipy (`scipy.stft`). This function is computationally efficient. The parameters are `segment_duration_sec` which defines the duration of the window in seconds onto which the spectra are calculated. The `bandwidth` kwarg is usedful to cut the spectra within a frequency band (in Hz), and therefore to reduce the memory usage.

Other keyword arguments are passed to the `scipy.stft` function. Default are: 
- Sampling rate `fs` (obtained from the stream directly)
- Number of dots per segments `nperseg` (obtained from the stream directly)
- Overlapping between segments `noverlap`: by default 50%
- Number of frequency points for calculating the Fourier transform `nfft` (by default next power of two of the segments length)
- Window shape `window` (by default Hanning)
- Calculate one side of the spectra `return_onesided` (by default `true` because signals are real, we do not need both parts because of hermitian symmetry)

The function returns three variables: the spectrograms (`shape = [n_traces, n_frequencies, n_times]`), the frequncies in Hz, and the times in absolute `matplotlib.dates format`.
```
import arrayprocessing as ap

stream = ap.read('/path/to/traces/*.sac')
# some preprocessing

spectrogram, frequencies, times = stream.stft(segment_duration_sec=16, bandwidth=[3, 10])
```


### Show speectrograms

You can show the spectrograms with the `spectrogram` method. Doc will come soon.

![](https://github.com/leonard-seydoux/arrayprocessing/blob/master/arrayprocessing_mindmap_data.png)

## Show speectrograms

You can show the spectrograms with the `spectrogram` method. Doc will come soon.



![](https://github.com/leonard-seydoux/arrayprocessing/blob/master/arrayprocessing_mindmap_covariance.png)

## correlation.py

![](https://github.com/leonard-seydoux/arrayprocessing/blob/master/arrayprocessing_mindmap_correlation.png)

## antenna.py

## synthetic.py

# To-do list

- Create a `calculate_coherence()` method with a kwarg `shanon_index` or `spectral_width`
- Enrich the `get_eigenvectors()` method with datetime selection and frequency range selection
- Create `self.z` set to 0 when antenna is 2D, otherwise set to th values available in txt or csv file.
- Add a warning message in synthetic.spherical if none of the `xyz` or `llz` kwargs is set. Clarify `spheri al` function.
