# Package structure of `covseisnet`

# Public objects

|Class name| Inheritance|Description|
|-|-|-|
|[`SynchroStream`](#synchrostream)|`obspy.core.stream.Stream`| Regular stream object with adapted methods.|
|`MultiStationSingleChannelCovariance`|`CovarianceMatrix`| Regular stream object with adapted methods.|
|`SingleStationMultiChannelCovariance`|`CovarianceMatrix`| Regular stream object with adapted methods.|
|`MultiStationMultiChannelCovariance`|`CovarianceMatrix`| Regular stream object with adapted methods.|

# Private objects

|Class name| Inheritance|Description|
|-|-|-|
|[`CovarianceMatrix`](#covariancematrix)|`np.ndarray`| Covariance matrix generic object.|
|[`CorrelationMatrix`](#correlationmatrix)|`np.ndarray`| Covariance matrix generic object.|


## `SynchroStream`
 > An `obspy.core.stream.Stream` object (array of `obspy.core.stream.Trace`) with the additional methods listed below.
- `synchonize()`: synchronize all traces on the same time sampling
- `times()`: returns only a single time vector
- `cut(start, end)`: select the traces within two dates
- `preprocess(**options)`: preprocess each trace in time or frequency domain


# Core objects

## `CovarianceMatrix`
 > A `numpy.ndarray` object with flexible physical dimension. The matrix dimension is always the two last dimensions. Other physical dimension (time, frequency) are vectorized, and the methods operate on all covariances at once. Selection of subsets of matrices in a `CovarianceMatrix` object is made via usual `numpy` slicing.
- `eigenvalues()`: get eigenvalues
- `eigenvectors(rank)`: get eigenvectors of given `rank`
- `coherence(kind)`: calculate coherence from spectral width or entropy.
- `_flat()`: (private) array vectorization

## `CorrelationMatrix`
 > Likewise `CovarianceMatrix`, it is a `numpy.ndarray` object with flexible physical dimension. The matrix dimension is always the two last dimensions. Other physical dimension (time, frequency) are vectorized, and the methods operate on all covariances at once. Selection of subsets of matrices in a `CovarianceMatrix` object is made via usual `numpy` slicing.
- `bandpass()`: bandpass filter
- `envelope()`: return correlation envelopes
- `smooth()`: smooth the envelopes
