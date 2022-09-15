# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 10:33:15 2022

@author: ftong
"""

import sys
import covseisnet as csn
import os
from obspy import UTCDateTime, read_inventory
import h5py
import numpy as np
import time
import xarray as xr
import dash
import plotly.graph_objects as go
import plotly.express as px
from dash import dcc, html
from dash.dependencies import Input, Output


def prep_archive(archive_path, spectral_width_length, duration_sec, sampling_rate, minute=0):


    # check for presence of archive file
    if os.path.exists(archive_path): #archive exists
        print("Loading from archive file...")
        
        # check archive file is compatible with run parameters
        with h5py.File(archive_path, "r") as hf:
        
            if spectral_width_length != hf['spectral width'].attrs['length']:
                print("Spectral width length of chosen parameters does not match that of the archive file. Please adjust the parameters to match or discard the archive file.")
                sys.exit()
                
            if sampling_rate != hf['spectral width'].attrs['sampling_rate']:
                print("Target sampling rate does not match that of the archive file. Please adjust the parameters to match or discard the archive file.")
                sys.exit()
    
            dset_time = hf['time']
            t_interval = dset_time.attrs['time interval']
            
            if t_interval == 0:
                print("Archive file is invalid. Please recreate!")
                sys.exit()
            
            t = UTCDateTime(dset_time[-1]) + t_interval #get start time of next window to process
            
        
    else: # Create archive file
        print("Specified archive file does not exist. Creating...")
    
        t = UTCDateTime.now() 
        t = t - duration_sec
        t = t.replace(minute=minute, second=0, microsecond=0) #set start time of archive to be rounded specified minute of the hour

        
        
        with h5py.File(archive_path, "w") as hf:
            dset_sw = hf.create_dataset("spectral width", (0, spectral_width_length), maxshape=(None, spectral_width_length), dtype='float')
            dset_sw.attrs['length'] = spectral_width_length
            dset_sw.attrs['sampling_rate'] = sampling_rate
            hf.create_dataset("network response function", (0, ), maxshape=(None,), dtype='float')
            
            utf8_type = h5py.string_dtype('utf-8', 30)
            dset_time = hf.create_dataset("time", (0, ), maxshape=(None,), dtype=utf8_type)
            dset_time.attrs['time interval'] = 0 #initialize with 0
            
    return t  


def fetch_data(data_source, data_protocol, t, duration_sec, delay=0, nrf=False, station_list=None, network_code=None, location_code=None, channel_code=None, metadata_filepath=None):
        
    if data_protocol == 'FDSN' or data_protocol == 'SeedLink':
        
        if data_protocol == 'FDSN':
        
            from obspy.clients.fdsn import Client
            
        else:
            
            from obspy.clients.seedlink.basic_client import Client
        
        # PRODUCTION    
        stream = csn.arraystream.ArrayStream()
        print('Connecting to', data_source, 'by', data_protocol)    
        # print(t)
        print("Requesting data from", t, "to", t + duration_sec)
        
        while (t + duration_sec > UTCDateTime.now() + delay):
            print("Data not yet available, waiting 60 seconds")
            time.sleep(60)
        
        client = Client(data_source)
        for sta in station_list:
            try:
                print('Requesting', sta)
                st = client.get_waveforms(network_code, sta, location_code, channel_code, t, t + duration_sec)
                stream.append(st[0])
            except:
                print('...Failed')
                continue
    
    
    elif data_protocol == 'File':         # DEVELOPMENT - get data from file instead of from server

        file = data_source
        print ("Loading waveforms from", file)

        stream = csn.arraystream.read(file)
        stream.trim(t, t + duration_sec)

    
    
    else:
        print("Unsupported data retrieval method:", data_protocol)
        sys.exit()
    
    
    print(stream.__str__(extended=True))
    
    if len(stream)==0:
        print("Failed to fetch waveforms")
        sys.exit()
         
        
    if metadata_filepath == None and data_protocol == 'FDSN': # get metadata from FDSN server
        client = Client(data_source)
        inv = client.get_stations(network=network_code, location=location_code, channel=channel_code, starttime=t, endtime=t + duration_sec)
    
    elif metadata_filepath == None:
        inv = None #ok for now since we have not implemented plotting of locations
    
    else:
        inv = read_inventory(metadata_filepath) #read metadata from file
            


    return stream, inv


def write_archive(archive_path, times, nwin, spectral_width, t, beam):

    print("Writing to archive file")
    
    # save data to h5f file
    with h5py.File(archive_path, "a") as hf:
            
        dset_sw = hf['spectral width']        
        dset_sw.resize(dset_sw.shape[0]+nwin, axis=0) #increase number of windows in archive
    
        #save timestamp (datetime) and window length (seconds)
        dset_time = hf['time']
        dset_time.resize(dset_time.shape[0]+nwin, axis=0) #increase number of windows in archive
        dset_time.attrs['time interval'] = times[1]
                
        dset_nrf = hf['network response function']
        dset_nrf.resize(dset_nrf.shape[0]+nwin, axis=0) #increase number of windows in archive
    
        #append new values of sw, nrf, and time stamps of windows. inefficient, improve later
        for i in range(0,nwin):
            dset_sw[dset_sw.shape[0]-nwin+i]=spectral_width[i]
            dset_time[dset_time.shape[0]-nwin+i]= t + times[i] 

            if beam != None:
                dset_nrf[dset_nrf.shape[0]-nwin+i]=beam.nrf[i]
            else:
                dset_nrf[dset_nrf.shape[0]-nwin+i]=0 #fill with zero

    

        

def read_archive(archive_path):
    
    if os.path.exists(archive_path):
                    
        with h5py.File(archive_path, "r") as hf:
            
            dset_sw = hf['spectral width']        
            dset_nrf = hf['network response function']
            dset_time = hf['time']
            sw_complete = dset_sw[:] #load entire available spectral width         
            nrf_complete = dset_nrf[:] #load entire available nrf
            time_complete = [UTCDateTime(dset_time[i]) for i in range(0,len(dset_time))]
            sampling_rate = dset_sw.attrs['sampling_rate']
            time_interval = dset_time.attrs['time interval']
            
        return sw_complete, nrf_complete, time_complete, sampling_rate, time_interval

    else:
        print("Archive file does not exist at specified path. Exiting...")
        sys.exit()




def plot_archive(archive_path, ymax, color_scale='jet_r', ymin=0.5,  refresh_secs=30, display_secs=604800): #tmax 604800 is 1 week 
    
    

    
        
    
    # display 
    es = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    
    app = dash.Dash(__name__, external_stylesheets=es)
    
    app.layout = html.Div(children=[
            
        html.Div([
            html.H4('Network Response Function'),
            html.Div(id='live-update-text-nrf'),
            dcc.Graph(id='live-update-graph-nrf'),
            dcc.Interval(
                id='interval-component-nrf',
                interval=refresh_secs*1000, # in milliseconds
                n_intervals=0
            )
        ]),
        
        html.Div([
            html.H4('Spectral Width'),
            html.Div(id='live-update-text'),
            dcc.Graph(id='live-update-graph'),
            dcc.Interval(
                id='interval-component',
                interval=refresh_secs*1000, # in milliseconds
                n_intervals=0
            )
        ])
    ])
    
    @app.callback(Output('live-update-graph', 'figure'),
                Input('interval-component', 'n_intervals'))
    def update_graph_live(n):
    
        sw_complete, nrf_complete, time_complete, sampling_rate, time_interval = read_archive(archive_path)
        

            
        #trim data set to desired time interval ver 2
        t_index = int(display_secs/time_interval)+1
        time_complete = time_complete[-1*t_index:]
        sw_complete = sw_complete[-1*t_index:]




        # Place data into an xarray
        sw = xr.DataArray(
            data=sw_complete,
            dims=["time", "frequency"],
            coords=dict(
                time=time_complete,
                frequency=np.linspace(0, sampling_rate, sw_complete.shape[1])
            ),
            attrs=dict(
                description="Covariance matrix spectral width",
            ),
        )
        fig = px.imshow(sw.T, color_continuous_scale=color_scale, origin='lower', aspect='auto')
        fig.update_layout(xaxis_title='Time', yaxis_title='Frequency (Hz)')
        fig.update_xaxes(tickson='boundaries', rangeselector_xanchor='right')
        fig.update_yaxes(type='log', range=[np.log10(ymin),np.log10(ymax)], autorange=False)
    
        return fig 
    
    @app.callback(Output('live-update-graph-nrf', 'figure'),
                Input('interval-component-nrf', 'n_intervals'))
    def update_graph_live_nrf(n):
    
        sw_complete, nrf_complete, time_complete, sampling_rate, time_interval = read_archive(archive_path)
                
        #trim data set to desired time interval
        t_index = int(display_secs/time_interval)+1
        time_complete = time_complete[-1*t_index:]
        nrf_complete = nrf_complete[-1*t_index:]

        xs = time_complete
        ys = nrf_complete
        fig = go.Figure(data=go.Scatter(x=xs, y=ys))
        fig.update_layout(xaxis_title='Time', yaxis_title='NRF, unnormalized')
    
        return fig 


    
    
    return app