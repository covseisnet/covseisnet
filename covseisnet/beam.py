#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Beamforming.

"""
import numpy as np
import sys

# class Beam(np.ndarray):
class Beam:

    # def __new__(cls, nwin, traveltime_grid):
    #     stockbeam = [np.zeros((traveltime_grid.nx, traveltime_grid.ny, traveltime_grid.nz)) for i in range(0,nwin)]
    #     return stockbeam
        # return super(Beam, cls).__new__(cls)
        # return np.zeros((nwin, traveltime_grid.nx, traveltime_grid.ny, traveltime_grid.nz))

    def __init__(self, shape):
    # def __init__(self):
        # print("inited")
        # self = np.zeros((traveltime.nx, traveltime.ny, traveltime.nz))
        # self = np.zeros(shape)
        self.grid = np.zeros(shape)
        # self.nrf =  np.zeros(shape)
        # self.likelihood = np.zeros((nwin, traveltime_grid.nx, traveltime_grid.ny, traveltime_grid.nz))
        # self.nrf = np.zeros(nwin)


    def set_extent(self, west, east, south, north, depth_top, depth_max):
        """ Limits of the beamforming domain.

        Args
        ----
            west, east, south, north, depth_top and depth_max (float):
                Extent of the 3D map in degrees for the azimuths and km for
                the depths. The depth is given in km and should be negative for
                elevation.
        """

        # self.extent = west, east, south, north, depth_top, depth_max
        self.xmin = west
        self.xmax = east
        self.ymin = south
        self.ymax = north
        self.zmin = depth_top
        self.zmax = depth_max
        self.lon = np.linspace(west, east, self.grid.shape[0])
        self.lat = np.linspace(south, north, self.grid.shape[1])
        self.dep = np.linspace(depth_top, depth_max, self.grid.shape[2])
        self.meshgrid = np.meshgrid(self.lon, self.lat, self.dep)
        self.meshgrid_size = len(self.grid[0].ravel())


    def max_likelihood(self):
        beam_max_index = np.nanargmax(self.grid) #1D index of max likelihood
        beam_max_indices = np.unravel_index(np.ravel_multi_index([beam_max_index], self.grid.flatten().shape), self.grid.shape) #get 3D index of max likelihood
        beam_max = self.grid[beam_max_indices] #return max likelihood value
        
        beam_max_lon = beam_max_indices[0]/self.grid.shape[0]*(self.xmax-self.xmin)+self.xmin
        beam_max_lat = beam_max_indices[1]/self.grid.shape[1]*(self.ymax-self.ymin)+self.ymin
        beam_max_depth = beam_max_indices[2]/self.grid.shape[2]*(self.zmax-self.zmin)+self.zmin
        
        return(beam_max_lon, beam_max_lat, beam_max_depth, beam_max)


    def calculate_nrf(self):
        return np.max(self.grid)*np.size(self.grid)/np.sum(self.grid) 

    def calculate_likelihood(self, xcorr, sampling_rate, traveltimes, close=None):
        """ Shift cross-correlation for each source in grid.
        xcorr.shape = (stat_combinations, lags)
        """

        # beam = np.zeros((T.nx, T.ny, T.nz)).view(csn.beam.Beam)

        # Initialization
        beam_max = 0
        trii, trij = np.triu_indices(traveltimes.nsta, k=1)
        # xcorr_shifted_best = 0
        n_lon = traveltimes.nx
        n_lat = traveltimes.ny
        n_dep = traveltimes.nz
        center = (xcorr.shape[1] - 1) // 2 + 1

        for k in range(n_lon * n_lat * n_dep):
            # Differential travel times
            # Complicated way of making a loop over lon, lat, dep
            i, j, k = np.unravel_index(k, (n_lon, n_lat, n_dep))
            # Extract the travel times of all stations from specific source at i, j ,k
            tt = traveltimes.grid[:, i, j, k]
            # Increase the dimension and find the substraction between all the stations in a NxN matrix
            tt = tt[:, None] - tt
            # Extract the upper triangle values (and invert sign)
            tt = -tt[trii, trij]

            if np.any(np.isnan(tt)):
                continue

            if close is not None:
                tt = tt[close]

            # Shift = center of the CC + arrival time
            dt_int = -(sampling_rate * tt).astype(int)
            dt_int = center + dt_int

            max_time_diff = np.max(np.abs(dt_int))
            
            if max_time_diff >= xcorr.shape[1]:
                sys.exit("ERROR: the max time difference " + str(max_time_diff) + " is bigger than the correlation duration " + str(xcorr.shape[1]))

            # beam is sum of the CCs after being shifted the arrival time difference
            # extract for each stat comb, the sum of the CCs with the delay
            beam = xcorr[range(xcorr.shape[0]), dt_int].sum()
            self.grid[i, j, k] = beam

            if beam_max < beam:
                # Keep that into memory
                beam_max = beam
                dt_int_abs = -(dt_int - center)

        # Move
        rows, column_indices = np.ogrid[:xcorr.shape[0], :xcorr.shape[1]]
        # Find where the time delay is bigger than the possible correlation time; this condition will never be fulfilled
        dt_int_abs[np.abs(dt_int_abs) > xcorr.shape[1]] = xcorr.shape[1] - 1
        # Move the negatives (almost all of them)
        dt_int_abs[dt_int_abs < 0] += xcorr.shape[1]
        column_indices = column_indices - dt_int_abs[:, np.newaxis]
        xcorr_best = xcorr[rows, column_indices]

        # return xcorr_best.T #don't return anything because xcorr not widely used