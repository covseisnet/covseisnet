#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Import externally created travel time grids for use in covseisnet.

"""
import numpy as np


class TravelTime:
    r"""List-like object of multiple travel time grids.
    
    This class is an object created by loading in travel time grids from numpy .npy files.
    """
    

    def __init__(self, stream, T_path, identifiers="S"):
        """Load a travel time grid for each trace in a given stream.
               
        Parameters
        ----------
        
        stream: :class:`ArrayStream` or :class:`obspy.core.stream.Stream`
            The input data stream. It will attempt to load a traveltime grid
            for each trace in the stream.
        
        T_path: str
            Path to the travel time grids.
            
        identifiers: str, default
            Format of the travel time grid filenames. Defaults to "S" for using
            the station name as the filename. For situations where multiple
            traces have a common station name, use "NSL" which requires the 
            network code, station name, and location code to all be present in 
            the filename.
 

        Example
        -------

        >>> import covseisnet as csn
        >>> import os
        >>> path = './tt_grids'
        >>> stream = csn.arraystream.read()
        >>> print(stream)
        3 Trace(s) in Stream:
        BW.RJOB..EHZ | 2009-08-24T00:20:05.000000Z... | 100.0 Hz, 701 samples
        BW.RJOB..EHN | 2009-08-24T00:20:05.000000Z... | 100.0 Hz, 701 samples
        BW.RJOB..EHE | 2009-08-24T00:20:05.000000Z... | 100.0 Hz, 701 samples
        >>> os.listdir(path)
        ['BW.RJOB..npy']
        >>> traveltimes = csn.traveltime.TravelTime(stream, T_path=path, identifiers="NSL")
        
        """
        
        self.nsta = len(stream)

        list_nsl = [
            (stream[i].stats.network, stream[i].stats.station, stream[i].stats.location)
            for i in range(self.nsta)
        ]

        T = []

        for nsl in list_nsl:
            network = nsl[0]
            station = nsl[1]
            location = nsl[2]
            
            try:
                # file names contain only the station
                if identifiers == "S":
                    T_sta = np.load(T_path + "/" + station + ".npy")
    
                # file names contain network, station, and location codes
                elif identifiers == "NSL":
                    T_sta = np.load(
                        T_path + "/" + network + "." + station + "." + location + ".npy"
                    )
                
                T.append(T_sta)

            except (FileNotFoundError):
                raise SystemExit("Travel time grid file not found for trace", nsl)


        T = np.array((T))
        self.grid = T
        self.nx = T.shape[1]
        self.ny = T.shape[2]
        self.nz = T.shape[3]
