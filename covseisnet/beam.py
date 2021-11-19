#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Beamforming.

"""
import numpy as np
import sys

class Beam(np.ndarray):

    def set_extent(self, west, east, south, north, depth_top, depth_max):
        """ Limits of the beamforming domain.

        Args
        ----
            west, east, south, north, depth_top and depth_max (float):
                Extent of the 3D map in degrees for the azimuths and km for
                the depths. The depth is given in km and should be negative for
                elevation.
        """

        self.extent = west, east, south, north
        self.lon = np.linspace(west, east, self.shape[0])
        self.lat = np.linspace(south, north, self.shape[1])
        self.dep = np.linspace(depth_top, depth_max, self.shape[2])
        self.grid = np.meshgrid(self.lon, self.lat, self.dep)
        self.grid_size = len(self.grid[0].ravel())


    def calculate_nrf(self):
        return np.max(self)*np.size(self)/np.sum(self)

    def calculate_heterogeneous(self, xcorr, fs, num_stn, ttimes, close=None):
        """ Shift cross-correlation for each source in grid.
        xcorr.shape = (stat_combinations, lags)
        """
        # t_a = time.time()
        # print(inventory.dim) #number of stations
        # Initialization
        beam_max = 0
        trii, trij = np.triu_indices(num_stn, k=1)
        # xcorr_shifted_best = 0
        n_lon = self.shape[0]
        n_lat = self.shape[1]
        n_dep = self.shape[2]
        center = (xcorr.shape[1] - 1) // 2 + 1

        # jaja = 0
        for k in range(n_lon * n_lat * n_dep):
            # Differential travel times
            # Complicated way of making a loop over lon, lat, dep
            i, j, k = np.unravel_index(k, (n_lon, n_lat, n_dep))
            # src_distance = ttimes['distances'][:, i, j]
            # sdid = [np.abs(ttimes['epicentral_distances'] - d).argmin()
            #         for d in src_distance]
            # tt = np.array([ttimes['ttimes'][s, sdid[s], k]
            #                for s in range(stations.dim)])
            # Extract the travel times of all stations from specific source at
            # i, j ,k
            tt = ttimes[:, i, j, k]
            # Increase the dimension and find the substraction between all the
            # stations in a NxN matrix
            tt = tt[:, None] - tt
            # Extract the upper triangle values (and invert sign)
            tt = -tt[trii, trij]

            if np.any(np.isnan(tt)):
                continue

            if close is not None:
                tt = tt[close]

            # Shift = center of the CC + arrival time
            dt_int = -(fs * tt).astype(int)
            dt_int = center + dt_int

            # if np.max(np.abs(dt_int)) > jaja:
            #     print(np.max(np.abs(dt_int)))
            #     jaja = np.max(np.abs(dt_int))


            max_time_diff = np.max(np.abs(dt_int))
            # print(dt_int)
            # print(xcorr.shape[1])
            if max_time_diff >= xcorr.shape[1]:
                sys.exit("ERROR: the max time difference " + str(max_time_diff) + " is bigger than the correlation duration " + str(xcorr.shape[1]))

            # beam is sum of the CCs after being shifted the arrival time
            # difference
            # extract for each stat comb, the sum of the CCs with the delay
            beam = xcorr[range(xcorr.shape[0]), dt_int].sum()
            self[i, j, k] = beam

            if beam_max < beam:
                # Keep that into memory
                beam_max = beam
                dt_int_abs = -(dt_int - center)

        # t_b = time.time()
        # print('Bbeam', t_b-t_a)
        # t_a = time.time()
        # Move
        rows, column_indices = np.ogrid[:xcorr.shape[0], :xcorr.shape[1]]
        # Find where the time delay is bigger than the possible correlation
        # time; this condition will never be fulfilled
        dt_int_abs[np.abs(dt_int_abs) > xcorr.shape[1]] = xcorr.shape[1] - 1
        # Move the negatives (almost all of them)
        dt_int_abs[dt_int_abs < 0] += xcorr.shape[1]
        column_indices = column_indices - dt_int_abs[:, np.newaxis]
        xcorr_best = xcorr[rows, column_indices]

        return xcorr_best.T