import covseisnet as csn
import pytest
import obspy

def test_arraystream():

    stream = csn.arraystream.read()

    # test non-existent method for preprocessing
    with pytest.raises(ValueError):
        stream.preprocess(domain='spectral', method='fake')
        
    with pytest.raises(ValueError):        
        stream.preprocess(domain='temporal', method='fake')

    # test non-existent domain type for preprocessing
    with pytest.raises(ValueError):        
        stream.preprocess(domain='fake')