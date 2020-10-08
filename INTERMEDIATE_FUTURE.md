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


## SynchroStream
> inherits from `obspy.core.stream.Stream` with additional methods listed below. 
- `synchonize()` synchronize all traces on the same time sampling 
- `times()` returns only a single time vector 
- `cut(start, end)` select the traces within two dates 
- `preprocess()` preprocess each trace in time or frequency domain
- `set_header()` attach metadata (if header exist, format it our way, else fetch DC)
- `header` our needed metadata

## Coordinates
> A set of coordinate points in a given coordinate system, with the following methods
- `argsort_distances()` sorting index as a function of inter-points distances
- `argsort_distances_from(point)` sorting index as a function of the specified point
- `get_distances()` returns all point to point distances
- `get_distances_from(point)` returns all grid point distance to specified point
- `project(projection)` from its own geographical coordinate system to the new coordinate system (based on `cartopy`)
- `vectorized()` returns a vectorized version of the object (looping acceleration)
- `coordinate_system` actual coordinate system of the object (e.g. `'geographical'`, `'UTM'`, `'cartesian'`)
- `keys`

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
