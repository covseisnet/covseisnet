#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Beamforming to calculate the likelihood function and network response function.

"""
import numpy as np
import sys


class Beam:
    r""" An object containing the likelihood function and associated products.
    
    A beam object is initiated by providing a :class:`~covseisnet.traveltime.TravelTime` object and the number of time windows used 
    """

    def __init__(self, nwin, traveltimes):
        """Create a Beam object by providing a traveltime object and specifying
        the number of time windows.
               
        Parameters
        ----------
        nwin: int
            The number of windows in the correlation matrix.
            
        traveltimes: :class:`~covseisnet.traveltime.TravelTime`
            An object containing all relevant travel time grids.
        

        """
        self.nwin = nwin
        self.nx = traveltimes.nx
        self.ny = traveltimes.ny
        self.nz = traveltimes.nz
        self.traveltimes = traveltimes
        self.likelihood = np.zeros((self.nwin, self.nx, self.ny, self.nz))
        self.nrf = np.zeros(self.nwin)
        self.correlation_shifted = []
        self.correlation_unshifted = []

    def set_extent(self, x_min, x_max, y_min, y_max, z_min, z_max):
        """Sets the geographical limits of the beamforming domain.

        Specify the extent of the 3D grid. Units need to match those of the 
        travel time grids. Horizontally the units may be in degrees of 
        longitude and latitude or in km. The depth is normally in km and should
        decrease for increasing elevation.

        Parameters
        ----------
        x_min: float 
            Minimum value in the X direction
            
        x_max: float
            Maximum value in the X direction
        
        y_min: float
            Minimum value in the Y direction
            
        y_max: float
            Maximum value in the Y direction
            
        z_min: float
            Minimum value in the Z direction
            
        z_max: float
            Maximum value in the Z direction
                
        """

        self.xmin = x_min
        self.xmax = x_max
        self.ymin = y_min
        self.ymax = y_max
        self.zmin = z_min
        self.zmax = z_max
        self.lon = np.linspace(x_min, x_max, self.nx)
        self.lat = np.linspace(y_min, y_max, self.ny)
        self.dep = np.linspace(z_min, z_max, self.nz)
        self.meshgrid = np.meshgrid(self.lon, self.lat, self.dep)
        self.meshgrid_size = len(self.likelihood[0][0].ravel())

    def max_likelihood(self, window_index):
        """ Extract maximum likelihood location
        
        Returns a tuple containing the maximum likelihood location and value in the window 
        specified by the window index.
               
        Parameters
        ----------
        window_index : int
            The window index of the desired likelihood function.
       
        Returns
        -------
        :class:`~tuple`              
            :class:`float`: The x coordinate of the maximum of the likelihood function. 
            
            :class:`float`: The y coordinate of the maximum of the likelihood function.
            
            :class:`float`: The z coordinate or depth of the maximum of the likelihood function.
            
            :class:`float`: The maximum value (unnormalized) of the likelihood function.
        

        """

        beam_max_index = np.nanargmax(
            self.likelihood[window_index]
        )  # 1D index of max likelihood
        beam_max_indices = np.unravel_index(
            np.ravel_multi_index(
                [beam_max_index], self.likelihood[window_index].flatten().shape
            ),
            self.likelihood[window_index].shape,
        )  # get 3D index of max likelihood
        beam_max = self.likelihood[window_index][
            beam_max_indices
        ]  # return max likelihood value

        beam_max_x = beam_max_indices[0] / self.nx * (self.xmax - self.xmin) + self.xmin
        beam_max_y = beam_max_indices[1] / self.ny * (self.ymax - self.ymin) + self.ymin
        beam_max_z = beam_max_indices[2] / self.nz * (self.zmax - self.zmin) + self.zmin

        return (beam_max_x, beam_max_y, beam_max_z, beam_max)

    def calculate_nrf(self, window_index):
        """ Calculates and returns the network response function for the window
        specified by the window index.
               
        Parameters
        ----------
        window_index: int
            The window index of the desired likelihood function upon which the
            network response function is to be calculated.
       
        Returns
        -------
        :class:`~numpy.ndarray`
            The value of the network response function in that window. 
        

        """
        self.nrf[window_index] = (
            np.max(self.likelihood[window_index])
            * np.size(self.likelihood[window_index])
            / np.sum(self.likelihood[window_index])
        )
        return self.nrf[window_index]

    def calculate_likelihood(
        self, cross_correlation, sampling_rate, window_index, close=None
    ):
        """Shift the cross-correlation for each source in the grid.
               
        Parameters
        ----------
        cross_correlation: :class:`~covseisnet.correlationmatrix.CorrelationMatrix`
            The correlation matrix from which the 3D likelihood function is 
            calculated.
            
        sampling_rate: float
            Sampling rate in Hz.
            
        window_index: int
            The index of the window upon which to calculate the 3D likelihood 
            function.
            
        close: 
        """

        # Initialization
        cross_correlation = cross_correlation.T
        self.correlation_unshifted.append(cross_correlation)

        beam_max = 0
        trii, trij = np.triu_indices(self.traveltimes.nsta, k=1)
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
                sys.exit(
                    "ERROR: the max time difference "
                    + str(max_time_diff)
                    + " is bigger than the correlation duration "
                    + str(cross_correlation.shape[1])
                )

            # beam is sum of the CCs after being shifted the arrival time difference
            # extract for each stat comb, the sum of the CCs with the delay
            beam = cross_correlation[range(cross_correlation.shape[0]), dt_int].sum()
            self.likelihood[window_index, i, j, k] = beam

            if beam_max < beam:
                # Keep that into memory
                beam_max = beam
                dt_int_abs = -(dt_int - center)

        # Move
        rows, column_indices = np.ogrid[
            : cross_correlation.shape[0], : cross_correlation.shape[1]
        ]

        # Find where the time delay is bigger than the possible correlation time; this condition will never be fulfilled
        dt_int_abs[np.abs(dt_int_abs) > cross_correlation.shape[1]] = (
            cross_correlation.shape[1] - 1
        )
        # Move the negatives (almost all of them)
        dt_int_abs[dt_int_abs < 0] += cross_correlation.shape[1]
        column_indices = column_indices - dt_int_abs[:, np.newaxis]
        self.correlation_shifted.append(cross_correlation[rows, column_indices])
