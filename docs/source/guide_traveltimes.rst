:: _guide_travel_times:

Importing travel time grids
===========================

Travel time grids saved as numpy .npy files can be imported into `covseisnet`. There are many ways to generate travel time grids, for example by using `TauP <https://docs.obspy.org/packages/obspy.taup.html>`_.

The :class:`~covseisnet.traveltime.TravelTime` class provides a list-like structure in which the traveltime grids can be loaded. Given a stream object, it will attempt to import a traveltime grid for each trace in the stream.