#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Covariance matrix in spectral domain.

Todo
----
- Pass the function :func:`covseisnet.covariance.calculate` as a top level
  function?
- Should we consider `privatizing` the :func:`covseisnet.covariance.xcov`
  routine?
"""

import numpy as np
import obspy

from scipy.linalg import eigvalsh, eigh


class CovarianceMatrix(np.ndarray):
    r"""Covariance matrix.

    This class is a subclass of :class:`numpy.ndarray`, meaning that
    all the methods available with regular numpy arrays are also available
    here, plus additional arrayprocessing-oriented methods. Note that
    any numpy method or function applied to a
    :class:`~arrayprocessing.covariance.CovarianceMatrix` instance returns a
    :class:`~arrayprocessing.covariance.CovarianceMatrix` instance.

    The shape of a covariance matrix calculate from :math:`N` traces is at
    least :math:`N \times N`. Depending on the averaging size and frequency
    content, the covariance matrix can be of shape

    - ``(n_sta, n_sta)`` if a single frequency and time sample is obtained.

    - ``(n_freq, n_sta, n_sta)`` for a single time sample and ``n_freq``
      frequency points

    - ``(n_times, n_freq, n_sta, n_sta)`` for ``n_times`` and ``n_freq``
      dimensions.

    All the methods defined in the the
    :class:`~arrayprocessing.covariance.CovarianceMatrix` class are performed
    on the flattened array with the private method
    :class:`~arrayprocessing.covariance.CovarianceMatrix._flat`, which allow
    to obtain as many :math:`N \times N` covariance matrices as time and
    frequency samples. Given a method that outputs a shape ``shape_out``
    given a covariance matrix of shape ``(n_sta, n_sta)``, the output of this
    method on a covariance matrix of shape ``(n_times, n_freq, n_sta, n_sta)``
    will be of shape ``(n_times, n_freq, *shape_out)``.

    Tip
    ---
    Any :class:`numpy.ndarray` can be turned into a
    :class:`~arrayprocessing.covariance.CovarianceMatrix` object with

    >>> import covseisnet as cn
    >>> import numpy as np
    >>> c = np.zeros((4, 4)).view(cn.CovarianceMatrix)
    >>> c
    CovarianceMatrix([[ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.]])
    """

    def __new__(cls, input_array):
        """Subclassing."""
        obj = np.asarray(input_array, dtype=complex).view(cls)
        return obj

    def coherence(self, kind="spectral_width", epsilon=1e-10):
        r"""Covariance-based coherence estimation.

        The measured is performed onto all the covariance matrices from
        the eigenvalues obtained with the method
        :meth:`~arrayprocessing.covariance.CovarianceMatrix.eigenvalues`.
        For a given matrix :math:`N \times N` matrix :math:`M` with
        eigenvalues :math:`\mathbf{\lambda} = \lambda_i` where
        :math:`i=1\ldots n`. The coherence is obtained as :math:`F(\lambda)`,
        with :math:`F` being defined by the `kind` parameter.

        - The spectral width is obtained with setting ``kind='spectral_width'``
          , and returns the width :math:`\sigma` of the
          eigenvalue distribution such as

          .. math::

              \sigma = \frac{\sum_{i=0}^n i \lambda_i}{\sum_{i=0}^n \lambda_i}


        - The entropy is obtained with setting ``kind='entropy'``, and returns
          the entropy :math:`h` of the eigenvalue distribution such as

          .. math::

              h = - \sum_{i=0}^n i \lambda_i \log(\lambda_i + \epsilon)


        Keyword arguments
        -----------------
        kind: str, optional
            The type of coherence, may be "spectral_width" (default) or
            "entropy".

        epsilon: float, optional
            The regularization parameter for log-entropy calculation. Default
            to ``1e-10``.

        Returns
        -------

        :class:`np.ndarray`
            The spectral width of maximal shape ``(n_times, n_frequencies)``.

        """
        if kind == "spectral_width":
            eigenvalues = self.eigenvalues(norm=sum)
            indices = np.arange(self.shape[-1])
            return np.multiply(eigenvalues, indices).sum(axis=-1)
        elif kind == "coherence":
            eigenvalues = self.eigenvalues(norm=sum)
            log_eigenvalues = np.log(eigenvalues + epsilon)
            return -np.sum(eigenvalues * log_eigenvalues, axis=-1)
        else:
            message = "{} is not an available option for kind."
            raise ValueError(message.format(kind))

    def eigenvalues(self, norm=max):
        """Eigenvalue decomposition.

        The eigenvalue decomposition is performed onto the two last dimensions
        of the :class:`~arrayprocessing.covariance.CovarianceMatrix` object.
        The function used for eigenvalue decomposition is
        :func:`scipy.linalg.eigvalsh`. It assumes that the input matrix is 2D
        and hermitian. The decomposition is performed onto the lower triangular
        part in order to save time.

        Keyword arguments
        -----------------
        norm : function, optional
            The function used to normalize the eigenvalues. Can be :func:`max`,
            (default), any other functions.

        Returns
        -------
        :class:`np.ndarray`
            The eigenvalues of maximal shape ``(n_times, n_freq, n_sta)``.

        """
        matrices = self._flat()
        eigenvalues = np.zeros((matrices.shape[0], matrices.shape[-1]))
        for i, matrix in enumerate(matrices):
            eigenvalues[i] = np.abs(eigvalsh(matrix)[::-1])
            eigenvalues[i] /= norm(eigenvalues[i])
        return eigenvalues.reshape(self.shape[:-1])

    def eigenvectors(self, rank=0, covariance=False):
        """Extract eigenvectors of given rank.

        The function used for extracting eigenvectors is
        :func:`scipy.linalg.eigh`. It assumes that the input matrix is 2D
        and hermitian. The decomposition is performed onto the lower triangular
        part.

        Keyword arguments
        -----------------
        rank : int, optional
            Eigenvector rank, 0 by default (first eigenvector).

        covariance: int, optional
            Outer-product of eigenvectors of rank ``rank``.

        Returns
        -------
        :class:`np.ndarray`
            The complex-valued eigenvector array of shape
            ``(n_times, n_freq, n_sta)`` if the parameter ``covariance`` is
            ``False``, else ``(n_times, n_freq, n_sta, n_sta)``.


        Todo
        ----
        Implement a new option on the ``rank`` in order to accept a list, so
        the filtered covariance can be obtained from multiple eigenvector
        ranks. This should be defined together with a ``normalize`` boolean
        keyword argument in order to take into account the eigenvalues or not,
        and therefore the isotropization of the covariance matrix would be
        here defined fully (so far, the user has to define a loop in the
        main script).

        """
        # Initialization
        matrices = self._flat()
        eigenvectors = np.zeros((matrices.shape[0], matrices.shape[-1]), dtype=complex)

        # Calculation over submatrices
        for i, m in enumerate(matrices):
            eigenvectors[i] = eigh(m)[1][:, -1 - rank]

        if covariance:
            ec = np.zeros(self.shape, dtype=complex)
            ec = ec.view(CovarianceMatrix)
            ec = ec._flat()
            for i in range(eigenvectors.shape[0]):
                ec[i] = eigenvectors[i, :, None] * np.conj(eigenvectors[i])
            ec = ec.reshape(self.shape)
            return ec.view(CovarianceMatrix)
        else:
            return eigenvectors.reshape(self.shape[:-1])

    def _flat(self):
        """Covariance matrices with flatten first dimensions.

        Returns
        -------

        :class:`np.ndarray`
            The covariance matrices in a shape ``(a * b, n, n)``.
        """
        return self.reshape(-1, *self.shape[-2:])

    def triu(self, **kwargs):
        """Extract upper triangular on flatten array.

        Keyword arguments
        -----------------
        **kwargs: dict, optional
            The keyword arguments passed to the :func:`numpy.triu` function.


        Example
        -------

        >>> import covseisnet as cn
        >>> import numpy as np
        >>> c = np.arange(8).reshape((2, 2, 2)).view(cn.CovarianceMatrix)
        >>> c
            CovarianceMatrix([[[0, 1],
                              [2, 3]],
                             [[4, 5],
                              [6, 7]]])
        >>> c.triu()
            CovarianceMatrix([[0, 1, 3],
                              [4, 5, 7]])

        """
        trii, trij = np.triu_indices(self.shape[-1], **kwargs)
        return self[..., trii, trij]


def calculate(stream, window_duration_sec, average, average_step=None, **kwargs):
    """Calculate covariance matrix from the streams.

    Arguments
    ---------
    stream: :class:`ArrayStream` or :class:`obspy.core.stream.Stream`
        The input data stream. If an :class:`obspy.core.stream.Stream` given,
        the calculation assumes a synchronized stream, and raises and error
        if not.

    window_duration_sec: float
        The Fourier calculation window in seconds.

    average: int
        The number of window used to estimate the sample covariance.

    Keyword arguments
    -----------------
    average_step: int, optional
        The sliding window step for covariance matrix calculation (in number
        of windows).

    **kwargs: dict, optional
        Additional keyword arguments passed to the
        :func:`~covseisnet.covariance.stft` function.

    Returns
    -------
    :class:`numpy.ndarray`
        The time vector of the beginning of each covariance window.

    :class:`numpy.ndarray`
        The frequency vector.

    :class:`covseisnet.covariance.CovarianceMatrix`
        The complex covariance matrix in a maximal shape
        ``(n_time, n_freq, n_sta, n_sta)``.


    Example
    -------
    Calculate the covariance matrix of the example stream with 1 second windows
    averaged over 5 windows:

    >>> import covseisnet as cn
    >>> stream = cn.data.read()
    >>> t, f, c = cn.covariance.calculate(stream, 1., 5)
    >>> print(c.shape)  # (n_times, n_freq, n_cha, n_cha)
        (28, 199, 3, 3)
    >>>> print(c[0, 1])  # first window, second frequency
    CovarianceMatrix([[  4.84823507e+08       +0.j        ,
                         8.36715570e+07+55825525.33551149j,
                        -2.47008776e+08-60728089.46938748j],
                      [  8.36715570e+07-55825525.33551149j,
                         1.67037017e+08       +0.j        ,
                        -1.82827585e+08+27408763.38879033j],
                      [ -2.47008776e+08+60728089.46938748j,
                        -1.82827585e+08-27408763.38879033j,
                         2.66448583e+08       +0.j        ]])


    """
    times, frequencies, spectra = stft(stream, window_duration_sec, **kwargs)

    # Parametrization
    step = average // 2 if average_step is None else average * average_step
    n_traces, n_windows, n_frequencies = spectra.shape

    # Times
    t_end = times[-1]
    times = times[:-1]
    times = times[: 1 - average : step]
    n_average = len(times)
    times = np.hstack((times, t_end))

    # Initialization
    cov_shape = (n_average, n_traces, n_traces, n_frequencies)
    covariance = np.zeros(cov_shape, dtype=complex)

    # Compute
    for t in range(n_average):
        covariance[t] = xcov(t, spectra, step, average)

    return (
        times,
        frequencies,
        covariance.view(CovarianceMatrix).transpose([0, -1, 1, 2]),
    )


def stft(
    stream,
    window_duration_sec,
    bandwidth=None,
    window_step_sec=None,
    window=np.hanning,
    times_kw=dict(),
    **kwargs
):
    """Short-time fourier transform.

    Arguments
    ---------
    stream: :class:`ArrayStream` or :class:`obspy.core.stream.Stream`
        The input data stream. If an :class:`obspy.core.stream.Stream` given,
        the calculation assumes a synchronized stream, and raises and error
        if not.

    window_duration_sec: float
        The Fourier calculation window in seconds.

    Keyword arguments
    -----------------
    window_step_sec: float, optional
        The step between two consecutive Fourier windows in seconds. Default
        (None) considers half the ``window_duration_sec``.

    bandwidth: list, optional
        The restricted frequencies onto which the Fourier transform is
        calculated. This improve performances when Fourier windows are large.

    times_kw: string, optional
        The keyword arguments passed to the
        :class:`covseisnet.data.ArrayStream.times` or
        :meth:`obspy.core.trace.Trace.times` method depending on the input
        stream class.

    window: function, optional
        The window function, by default :func:`numpy.hanning`. A list of
        available window functions is available at
        https://numpy.org/doc/stable/reference/routines.window.html. Any
        function taking as argument a integer is accepted.

    **kwargs: dict, optional
        Additional keyword arguments passed to the :func:`numpy.fft.fft`
        function. By default, ``kwargs['n']`` is defined as
        ``2 * npts - 1`` in order to calculate the cross-correlation.

    Returns
    -------
    :class:`numpy.ndarray`
        The time vector of the beginning of each Fourier window.

    :class:`numpy.ndarray`
        The frequency vector.

    :class:`numpy.ndarray`
        The complex spectrograms of each timeseries in a
        ``(n_sta, n_time, n_freq)`` array.
    """
    # Time vector
    fs = stream[0].stats.sampling_rate
    npts = int(window_duration_sec * fs)
    step = npts // 2 if window_step_sec is None else int(window_step_sec * fs)
    times_kw.setdefault("type", "relative")
    times_kw.setdefault("reftime", None)
    if type(stream) is obspy.core.stream.Stream:
        times = stream[0].times(**times_kw)[: 1 - npts : step]
    else:
        times = stream.times(**times_kw)[: 1 - npts : step]
    n_times = len(times)

    # Frequency vector
    kwargs.setdefault("n", 2 * npts - 1)
    frequencies = np.linspace(0, fs, kwargs["n"])

    if bandwidth is not None:
        fin = (frequencies >= bandwidth[0]) & (frequencies <= bandwidth[1])
    else:
        fin = np.ones_like(frequencies, dtype=bool)
    frequencies = frequencies[fin]

    # Calculate spectra
    spectra_shape = len(stream), n_times, sum(fin)
    spectra = np.zeros(spectra_shape, dtype=complex)
    for trace_id, trace in enumerate(stream):
        tr = trace.data
        for time_id in range(n_times):
            start = time_id * step
            end = start + npts
            segment = tr[start:end] * window(npts)
            spectra[trace_id, time_id] = np.fft.fft(segment, **kwargs)[fin]

    # Times are extended with last time of traces
    t_end = stream[0].times(**times_kw)[-1]
    times = np.hstack((times, t_end))
    return times, frequencies, spectra


def xcov(wid, spectra_full, overlap, average):
    """Calculation of the array covariance matrix from the array data vectors.

    Warning
    -------
    This function is not fully documented yet, as may be reformulated.

    To do
    -----
    Allow to possibly reformulate this part with Einstein convention for
    faster computation, clarity and TensorFlow GPU transparency.

    Arguments
    ---------
    spectra_full: :class:`numpy.ndarray`
        The stream's spectra.

    overlap: int
        The average step.

    average: int
        The number of averaging windows.

    Returns
    -------
    :class:`numpy.ndarray`
        The covariance matrix.
    """
    n_traces, n_windows, n_frequencies = spectra_full.shape
    beg = overlap * wid
    end = beg + average
    spectra = spectra_full[:, beg:end, :].copy()
    x = spectra[:, None, 0, :] * np.conj(spectra[:, 0, :])
    for swid in range(1, average):
        x += spectra[:, None, swid, :] * np.conj(spectra[:, swid, :])
    return x
