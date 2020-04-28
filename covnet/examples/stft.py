#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Covariance matrix in spectral domain."""

import covnet as cn

stream = cn.data.read()
t, f, c = cn.covariance.calculate(stream, 1., 5)
print(c.shape)
print(c[0, 1])
print(c.shape)
