#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Beamforming.

"""
import numpy as np
import sys

class Beam:


    def __init__(self, nwin, traveltimes):
        self.nwin = nwin
        self.nx = traveltimes.nx
        self.ny = traveltimes.ny
        self.nz = traveltimes.nz
        self.traveltimes = traveltimes
        self.likelihood = np.zeros((self.nwin, self.nx, self.ny, self.nz))
        self.nrf = np.zeros(self.nwin)


    def set_extent(self, west, east, south, north, depth_top, depth_max):
        """ Limits of the beamforming domain.

        Args
        ----
            west, east, south, north, depth_top and depth_max (float):
                Extent of the 3D map in degrees for the azimuths and km for
                the depths. The depth is given in km and should be negative for
                elevation.
        """

        self.xmin = west
        self.xmax = east
        self.ymin = south
        self.ymax = north
        self.zmin = depth_top
        self.zmax = depth_max
        self.lon = np.linspace(west, east, self.nx)
        self.lat = np.linspace(south, north, self.ny)
        self.dep = np.linspace(depth_top, depth_max, self.nz)
        self.meshgrid = np.meshgrid(self.lon, self.lat, self.dep)
        self.meshgrid_size = len(self.likelihood[0][0].ravel())


    def max_likelihood(self, window_index):
        beam_max_index = np.nanargmax(self.likelihood[window_index]) #1D index of max likelihood
        beam_max_indices = np.unravel_index(np.ravel_multi_index([beam_max_index], self.likelihood[window_index].flatten().shape), self.likelihood[window_index].shape) #get 3D index of max likelihood
        beam_max = self.likelihood[window_index][beam_max_indices] #return max likelihood value
        
        beam_max_lon = beam_max_indices[0]/self.nx*(self.xmax-self.xmin)+self.xmin
        beam_max_lat = beam_max_indices[1]/self.ny*(self.ymax-self.ymin)+self.ymin
        beam_max_depth = beam_max_indices[2]/self.nz*(self.zmax-self.zmin)+self.zmin
        
        return(beam_max_lon, beam_max_lat, beam_max_depth, beam_max)


    def calculate_nrf(self, window_index):
        self.nrf[window_index] = np.max(self.likelihood[window_index])*np.size(self.likelihood[window_index])/np.sum(self.likelihood[window_index]) 
        return self.nrf[window_index]

    def calculate_likelihood(self, cross_correlation, sampling_rate, window_index, close=None):
        """ Shift cross-correlation for each source in grid.
        cross_correlation.shape = (stat_combinations, lags)
        """

        # Initialization
        cross_correlation = cross_correlation.T
        beam_max = 0
        trii, trij = np.triu_indices(self.traveltimes.nsta, k=1)
        # cross_correlation_shifted_best = 0
        n_lon = self.nx
        n_lat = self.ny
        n_dep = self.nz
        center = (cross_correlation.shape[1] - 1) // 2 + 1

        for k in range(n_lon * n_lat * n_dep):
            # Differential travel times
            # Complicated way of making a loop over lon, lat, dep
            i, j, k = np.unravel_index(k, (n_lon, n_lat, n_dep))
            # Extract the travel times of all stations from specific source at i, j ,k
            tt = self.traveltimes.grid[:, i, j, k]
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
            
            if max_time_diff >= cross_correlation.shape[1]:
                sys.exit("ERROR: the max time difference " + str(max_time_diff) + " is bigger than the correlation duration " + str(cross_correlation.shape[1]))

            # beam is sum of the CCs after being shifted the arrival time difference
            # extract for each stat comb, the sum of the CCs with the delay
            beam = cross_correlation[range(cross_correlation.shape[0]), dt_int].sum()
            self.likelihood[window_index, i, j, k] = beam

            if beam_max < beam:
                # Keep that into memory
                beam_max = beam
                dt_int_abs = -(dt_int - center)

        # Move
        rows, column_indices = np.ogrid[:cross_correlation.shape[0], :cross_correlation.shape[1]]
        # Find where the time delay is bigger than the possible correlation time; this condition will never be fulfilled
        dt_int_abs[np.abs(dt_int_abs) > cross_correlation.shape[1]] = cross_correlation.shape[1] - 1
        # Move the negatives (almost all of them)
        dt_int_abs[dt_int_abs < 0] += cross_correlation.shape[1]
        column_indices = column_indices - dt_int_abs[:, np.newaxis]
        cross_correlation_best = cross_correlation[rows, column_indices]

        return self.likelihood[window_index]
        # return cross_correlation_best.T #don't because cross_correlation not widely used