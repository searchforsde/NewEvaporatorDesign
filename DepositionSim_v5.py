# -*- coding: utf-8 -*-
"""
This script models a nozzle at (x,y) = (0,0) emitting with a cos(theta)**n 
angular density and depositing on a rectangular substrate of (width, length) 
and at height away from the nozzle.
"""

import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
# import matplotlib.image as mpimg

# total flux is the integral of c*cos(theta)**n over the hemisphere of emission
c = 1
n = 2

# boundaries of the rectangular substrate
x_low = -50
x_high = 50
y_low = 0
y_high = 1000

# resolve down to "pixels" of this dimension on the substrate
resolution = 2

# height above the nozzle tip
height = 100

# define walls attached to the nozzles with surfaces parallel to the x-axis
left_wall_x_pos = -7
right_wall_x_pos = 7
wall_height = 20

# nozzle starting positions
nozzles_pos = np.array([[-7, 7],
                        [-7,-7],
                        [0,0],
                        [7,7],
                        [7,-7],
                        [-7,0],
                        [7,0]], dtype='f')
#nozzles_pos = np.array([[0, 0]], dtype='f')

# total nozzle array position change
nozzle_pos_change = np.array([0, 1000], dtype='f')

# length of time for depostion
t_span = 10 ** -3

# time step for integration purposes
dt = 10 ** -5

# %% All variables have been set. Now process data

nozzle_count = nozzles_pos.shape[0] # number of nozzles

# define x and y position of every substrate pixel
x = np.arange(x_low, x_high + resolution, resolution)
y = np.arange(y_low, y_high + resolution, resolution)
sub_x, sub_y = np.meshgrid(x, y)

sub_deposition = np.zeros(sub_x.shape,dtype='f') # deposition map

dxy = nozzle_pos_change * dt/t_span # [x,y] position change of the nozzle each time step

t = 0 # set starting time

# %%
while t <= t_span:
    
    for nozzle_idx in range(nozzle_count):
        nozzle_pos = nozzles_pos[nozzle_idx] # position of the current nozzle
        
        # effective x-cutoff of the walls for the current nozzle
        left_cutoff_pos = nozzle_pos[0] + (left_wall_x_pos - nozzle_pos[0])*height/wall_height
        right_cutoff_pos = nozzle_pos[0] + (right_wall_x_pos - nozzle_pos[0])*height/wall_height
        
        # l = left cutoff index
        try:
            l = np.where(sub_x[0] <= left_cutoff_pos)[0][-1] 
        except:
            l = 0
        # r = right cutoff index
        try:
            r = np.where(sub_x[0] >= right_cutoff_pos)[0][0]+1
        except:
            r = sub_deposition.shape[1]
        
        # define every pixel's distance from the nozzle
        dist_array = np.sqrt((sub_x[:,l:r] - nozzle_pos[0]) ** 2 + (sub_y[:,l:r] - nozzle_pos[1]) ** 2 + height ** 2) 
        # angle from nozzle axis to each pixel
        cos__array = np.divide(height, dist_array) 
        # add deposition to each pixel
        sub_deposition[:,l:r] = sub_deposition[:,l:r] + dt * c * cos__array ** n / dist_array ** 2
    t += dt # increment the time
    nozzles_pos += dxy # change in nozzle position
    
sub_deposition = sub_deposition / np.nanmax(sub_deposition) * 100

plt.figure(figsize=(40, 20))
plt.subplot(211)
ax1 = plt.imshow(sub_deposition, cmap='magma', extent = [x_low, x_high, y_low, y_high]) # starter plot. Needs to be scaled and labeled.
cbar = plt.colorbar(ax1)   
#cbar.ColorbarBase(ax1, cmap=cm, norm=plt.colors.Normalize(vmin=0, vmax=100))
cbar.ax.set_yticklabels(['20%','30%','40%','50%','60%','70%','80%','90%','100%'])

plt.figure(figsize=(1, 10))
levels = np.array([90,95,99])
ax2 = plt.contour(sub_x, sub_y, sub_deposition, levels, colors='k', extent = [x_low, x_high, y_low, y_high]) # can switch between color and black using this: cmap='magma'
fmt = {}
strs = [str(level) + '%' for level in levels]
for l, s in zip(ax2.levels, strs):
    fmt[l] = s
plt.clabel(ax2, inline=1, fontsize=15, fmt=fmt)


