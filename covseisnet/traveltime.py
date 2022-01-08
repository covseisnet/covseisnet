#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Travel time grids.

"""
import numpy as np

class TravelTime:
    """List-like object of multiple travel time grids.



    """
    def __init__(self, stream, T_path, identifiers="S",):
        self.nsta = len(stream)
        
        list_nsl = [(stream[i].stats.network, stream[i].stats.station, stream[i].stats.location) for i in range(self.nsta)]

        T = []
        
        for nsl in list_nsl:
            network = nsl[0]
            station = nsl[1]
            location = nsl[2]
            
            #file names contain only the station
            if identifiers == "S": 
                T_sta = np.load(T_path+'/'+station+'.npy')

            #file names contain network, station, and location codes
            elif identifiers == "NSL":
                T_sta = np.load(T_path+'/'+network+'.'+station+'.'+location+'.npy') 

            T.append(T_sta)
            
        T = np.array((T))
        self.grid = T
        self.nx = T.shape[1]
        self.ny = T.shape[2]
        self.nz = T.shape[3]