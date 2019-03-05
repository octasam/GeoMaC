#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 10:14:00 2019

@author: luiggi
"""

import maya

a = float(input('dame el extremo a:'))
b = float(input('dame el extremo b:'))
c = 0
d = 10

Nx = int(input('Nx:'))
Ny = int(input('Ny:'))

x, y = maya.mesh(a,b,c,d,Nx,Ny)

print(x,y)