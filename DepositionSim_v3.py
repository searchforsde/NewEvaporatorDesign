# -*- coding: utf-8 -*-
"""
This script models a nozzle at (x,y) = (width/2,0) emitting with a cos(theta)**n 
angular density and depositing on a rectangular substrate of (width, length) 
and at height away from the nozzle.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# total flux is the integral of c*cos(theta)**n over the hemisphere of emission
c = 1
n = 2

# dimensions of the rectangular substrate
width = 100 # x dimension
length = 1000 # y dimension

# height above the nozzle tip
height = 100

# resolve down to "pixels" of this dimension on the substrate
resolution = 0.1

# [x, y] position of the substrate centroid
rectangle_centroid = np.array([width/2, length/2], dtype='f')

# nozzle starting position
nozzle_start = np.array([width/2, 0], dtype='f')

# nozzle final position
nozzle_final = np.array([width/2, length], dtype='f')

# length of time for depostion
t_span = 10 ** -3

# time step for integration purposes
dt = 10 ** -5


# All variables have been set. Now process data

pixelsnumber_width = int(width / resolution) # pixel count along x
pixelsnumber_length = int(length / resolution) # pixel count along y
sub_deposition = np.zeros([pixelsnumber_length, pixelsnumber_width]) # deposition map

# sub_x is the x-position of every substrate pixel
sub_x = np.arange(pixelsnumber_width) * resolution # starts out as a vector
sub_x = sub_x + rectangle_centroid[0] - np.mean(sub_x) # correct the mean value
sub_x = np.repeat(np.expand_dims(sub_x, 0), pixelsnumber_length, axis = 0) # expand into an array

# sub_y is the y-position of every substrate pixel
sub_y = np.arange(pixelsnumber_length) * resolution # starts out as a vector
sub_y = np.flip(sub_y) + rectangle_centroid[1] - np.mean(sub_y) # correct the orientation and mean value
sub_y = np.repeat(np.expand_dims(sub_y, 1), pixelsnumber_width, axis = 1) # expand into an array

# every pixel's distance from the nozzle
dist_array = np.sqrt((sub_x - nozzle_start[0]) ** 2 + (sub_y - nozzle_start[1]) ** 2 + height ** 2)

# cos(theta) from the nozzle midline to each pixel
cos__array = np.divide(height, dist_array)

dxy = (nozzle_final - nozzle_start) * dt/t_span

t = 0 # set starting time
while t <= t_span:
    # add deposition to each pixel
    sub_deposition = sub_deposition + np.divide(dt * c * cos__array ** n, dist_array ** 2) 
    t += dt # increment the time
    nozzle_start += dxy # change in nozzle position
    dist_array = np.sqrt((sub_x - nozzle_start[0])** 2 + (sub_y - nozzle_start[1])** 2 + height ** 2) # update distances from nozzle
    cos__array = np.divide(height, dist_array) # update cos(theta) from nozzle midline
    
sub_deposition = sub_deposition / np.nanmax(sub_deposition) * 100

x_min = rectangle_centroid[0] - width/2
x_max = rectangle_centroid[0] + width/2
y_min = rectangle_centroid[1] - length/2
y_max = rectangle_centroid[1] + length/2

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1 = plt.imshow(sub_deposition, cmap='magma', extent = [x_min, x_max, y_min, y_max]) # starter plot. Needs to be scaled and labeled.
cbar = fig.colorbar(ax1)