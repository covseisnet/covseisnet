import covseisnet as csn


def test_covariance_spec_width():

    import obspy

    # read ObsPy's example stream
    stream = obspy.read()
    window_duration_sec = 1.0
    average = 5
    times, frequencies, covariances = csn.covariancematrix.calculate(
        stream, window_duration_sec, average
    )

    assert covariances.all() != None

    spectral_width = covariances.coherence(kind="spectral_width")

    assert spectral_width.all() != None
