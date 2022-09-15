#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
CovSeisNet Real-Time
================================================================================


"""
import covseisnet as csn
import numpy as np
import csn_rt as csn_rt
import argparse
import sys

parser = argparse.ArgumentParser(description="Run covseisnet continuously and save results to a hdf5 file. Run in display mode to view plots in a web browser.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('archive_path', help="path to the hdf5 archive file")

parser.add_argument('--nrf',  action="store_true", help="enable computation of Network Response Function")
parser.add_argument('--data_source', help="where to get the data from, specified as URL or FDSN web service name or file path")
parser.add_argument('-p', '--protocol', default='FDSN', help="how to get waveform data: 'FDSN' (default), 'SeedLink', or 'File'")
parser.add_argument('-g', '--grid_path', help="path to the directory containing the traveltime grids")
parser.add_argument('-m', '--metadata_path', help="path to the station xml file containing station metadata")
parser.add_argument('-r', '--sampling_rate', type=float, default=20, help="Data will be resampled to this value (Hz) before processing")
parser.add_argument('-d', '--duration', type=int, default=3600, help="Length of data (seconds) to process")
parser.add_argument('-w', '--window_duration_sec', type=float, default=60, help="length of Fourier calculation window (seconds)")
parser.add_argument('-a', '--average', type=int, default=30, help="number of window used to estimate the sample covariance")
parser.add_argument('-o', '--overlap', type=float, default=0.5, help="overlap of windows, default is 50 percent")
parser.add_argument('--sigma', type=float, default=20, help="standard deviation in correlation smoothing")
parser.add_argument('--delay', type=int, default=0, help="time (seconds) to wait for arrival of data; to compensate for expected data latency")
parser.add_argument('-f', '--filter_corner_frequencies', default='0.5,10', help="bandpass filter corner frequencies")
parser.add_argument('-n', '--network', default='*', help="Network code to retrieve data for")
parser.add_argument('-s', '--station', default='*', help="Station code to retrieve data for")
parser.add_argument('-l', '--location', default='*', help="Location code to retrieve data for")
parser.add_argument('-c', '--channel', default='*', help="Channel code to retrieve data for")

#arguments for plotting
parser.add_argument('--display',  action="store_true", help="run in display mode to plot data in web browser")
parser.add_argument('--refresh', type=int, default=30, help="refresh interval (seconds) of display")
parser.add_argument('--color_scale', default='jet_r', help="Color scale to use for spectral width")
parser.add_argument('--ymin', type=float, default=0.5, help="y-axis lower limit, default is 0.5Hz")
parser.add_argument('--ymax', type=float, default=None, help="y-axis upper limit, default is Nyquist frequency")
parser.add_argument('--tmax', type=int, default=604800, help="maximum length of data (in seconds) to display, default is 1 week")


args = parser.parse_args()

data_source = args.data_source
archive_path = args.archive_path
data_protocol = args.protocol
sampling_rate = args.sampling_rate
network = args.network
channel = args.channel
location = args.location
list_stations = [s.strip() for s in args.station.split(",")] #User may provide multiple stations in a comma delimited list
# traveltime_filepath = args.grid_path
# metadata_filepath = args.metadata_path
duration_sec = args.duration
low_pass, high_pass = [float(s.strip()) for s in args.filter_corner_frequencies.split(",")]
window_duration_sec = args.window_duration_sec
average = args.average
overlap= args.overlap
sigma = args.sigma
delay = args.delay
refresh_secs = args.refresh
display_secs = args.tmax
ymax = args.ymax

if data_source == None and args.display == False:
    print("No data source provided, therefore running in Display Mode...")

if args.display or data_source==None:
    print("Display Mode activated")
    
    if ymax==None:
        ymax = sampling_rate/2    #if user does not specifiy a ymax for plot, set to half of sampling rate 
    
    app = csn_rt.plot_archive(archive_path, ymax, color_scale=args.color_scale, ymin=args.ymin, refresh_secs=refresh_secs, display_secs=display_secs)
    
    if __name__ == '__main__':
        app.run_server(debug=True, host='0.0.0.0')
        # app.run_server(host='0.0.0.0', debug=False, use_reloader=True, dev_tools_hot_reload=True) #detects changes in code and auto-reloads and replots
    
    sys.exit()
    
# CHECK TT path provided if NRF option checked
if args.nrf:
    print("Network Response Function option enabled")
    if args.grid_path==None:
        print("Path to traveltime grids, using --grid_path , is missing. Exiting...")
        sys.exit()
        
else: #set beam to null since we will not be computing it
    beam = None

preproc_spectral_secs = window_duration_sec * average * overlap
spectral_width_length = sampling_rate * window_duration_sec * 2 - 1

t = csn_rt.prep_archive(archive_path, spectral_width_length, duration_sec, sampling_rate)
print("Data starting at", t, "will be processed.")
    
# Retrieve waveforms and metadata
stream, inv  = csn_rt.fetch_data(data_source, data_protocol, t, duration_sec, delay=delay, nrf=args.nrf , station_list=list_stations, network_code=network, location_code=location, channel_code=channel)

# resample data
stream.resample(sampling_rate)

# Pre-processing
print('Pre-processing waveforms')   
               
## merge traces to have one trace per station
stream.merge(method=1,fill_value='interpolate',interpolation_samples=-1)

n_sta = len(stream)
      
## synchronize traces in the stream
stream = stream.synchronize(t, duration_sec, method="linear")

## filtering
stream.detrend(type='demean')
stream.detrend(type='linear')
stream.filter(type='bandpass', freqmin=low_pass, freqmax=high_pass)  

## Preprocess
stream.preprocess(domain='spectral', method='onebit', window_duration_sec=preproc_spectral_secs)
stream.trim(t, t+duration_sec) #preprocessing can add artifacts to ends

# Calculate Covariance matrix and wavefield coherence
# ---------------------------------------------------
    
print('Computing Covariance Matrix ')   


times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)


# Spectral width

print('Computing Spectral Width  ')   
spectral_width = covariances.coherence(kind="spectral_width")

print('Extracting Cross-Correlations ')  
   
# Eigenvector decomposition - covariance matrix filtered by the 1st eigenvector to show the dominant source
covariance_1st = covariances.eigenvectors(covariance=True, rank=0)

# Extract cross-correlations
lags_1st, correlation_1st =  csn.correlationmatrix.cross_correlation(covariance_1st, sampling_rate)

lags_full, correlation_full =  csn.correlationmatrix.cross_correlation(covariances, sampling_rate)
  
nwin = correlation_full.nwin() # number of time windows

# nrf flag selected, therefore do beamforming and compute nrf
if args.nrf:
    print('Computing the Network Response Function and Likelihood Function')  
    
    # load traveltime grids from file
    traveltimes = csn.traveltime.TravelTime(stream, args.grid_path)
    
    # compute NRF and likelihood function
    beam = csn.beam.Beam(nwin, traveltimes)
    
    nech = correlation_full.shape[1]
    npairs = correlation_full.shape[2]
    ech_center = (nech-1) // 2 + 1
    
    R = np.zeros([nwin,traveltimes.nx, traveltimes.ny, traveltimes.nz])
    delta_t = np.zeros([npairs,traveltimes.nx, traveltimes.ny, traveltimes.nz])
    
    n=0
    while n <npairs:
        for a in range(0,n_sta):
            for b in range(a+1,n_sta):
                delta_t[n] = traveltimes.grid[a] - traveltimes.grid[b]
                n=n+1
    
    delta_t2 = np.where(delta_t<0,np.floor(delta_t),np.ceil(delta_t))
    delta_int = delta_t2.astype(np.int64)
    
    NRF=np.zeros((nwin))
    
    # Loop through all windows and calculate likelihood and nrf
    for i in range(0, nwin):
        print('Processing window', i+1, 'of', nwin)
    
        correl = correlation_full[i]
    
        # Filter correlation
        correl = correl.bandpass(low_pass, high_pass, sampling_rate)
        
        # Smooth correlation
        correl = correl.hilbert_envelope()
        correl = correl.smooth(sigma=sigma) #default sigma is 5
        
        beam.calculate_likelihood(correl, sampling_rate, i)
        beam.calculate_nrf(i) 
        
        # beam_max = beam.max_likelihood(i)

#write new data to archive file
csn_rt.write_archive(archive_path, times, nwin, spectral_width, t, beam)

print('End of computation')