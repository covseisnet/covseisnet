#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa

__all__ = [
    "covariancematrix",
    "arraystream",
    "beam",
    "correlationmatrix",
    "traveltime",
]

from . import *
from .covariancematrix import CovarianceMatrix
from .correlationmatrix import CorrelationMatrix
from .traveltime import TravelTime
from .beam import Beam
from .arraystream import read
from pkg_resources import get_distribution

__version__ = get_distribution('covseisnet').version