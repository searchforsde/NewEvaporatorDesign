# -*- coding: utf-8 -*-
"""
This script models a nozzle at (x,y,z) = (0,0,0) emitting with a cos(c*theta) 
angular density (to pi/2 rad) and depositing on a rectangular substrate of 
w = delta_x, l = delta_y, and at height z above the nozzle.
"""
# total flux is the integral of c*cos^n(theta) over the hemisphere of emission
c = 1  # constant in flux
n = 1  # power of cosine

# dimensions of the rectangular substrate
w = 104.9 # x dimension
l = 10000 # y dimension

# height above the nozzle tip
height = 50

# resolve down to "pixels" of this dimension on the substrate
resolution = 0.1

# starting [x, y] position of the substrate centroid
rectangle_centroid = [-50,0]

# change in (x, y) position of the substrate at the end of the deposition time
position_change = [100,0]

# length of time for depostion
t_span = 1

# time step for integration purposes
dt = 10 ** -3

# All variables have been set. Now process data
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

pixels_w = int(w / resolution) # pixel count along x
pixels_l = int(l / resolution) # pixel count along y
sub_deposition = np.zeros([pixels_w, pixels_l]) # deposition map

# sub_x is the x-position of every substrate pixel
sub_x = np.arange(pixels_w) * resolution # starts out as a vector
sub_x = sub_x + rectangle_centroid[0] - np.mean(sub_x) # correct the mean value
sub_x = np.repeat(np.expand_dims(sub_x, 0), pixels_l, axis = 0) # expand into an array

# sub_y is the y-position of every substrate pixel
sub_y = np.arange(pixels_l) * resolution # starts out as a vector
sub_y = np.flip(sub_y) + rectangle_centroid[1] - np.mean(sub_y) # correct the orientation and mean value
sub_y = np.repeat(np.expand_dims(sub_y, 1), pixels_w, axis = 1) # expand into an array

dist_array = sub_x ** 2 + sub_y ** 2 + height ** 2 # every pixel's distance from the nozzle

cos__array = np.divide(height, dist_array) # cos(theta) from the nozzle midline to each pixel

t = 0 # set starting time
while t < t_span:
    sub_deposition = sub_deposition + np.divide(dt * c * cos__array ** (n + 1), dist_array ** 2) # add deposition to each pixel
    t += dt # increment the time
    sub_x = sub_x + position_change[0] * dt / t_span # update x-positions
    sub_y = sub_y + position_change[1] * dt / t_span # update y-positions
    dist_array = sub_x ** 2 + sub_y ** 2 + height ** 2 # update distances from nozzle
    cos__array = np.divide(height, dist_array) # update cos(theta) from nozzle midline
    
imgplot = plt.imshow(sub_deposition) # starter plot. Needs to be scaled and labeled.