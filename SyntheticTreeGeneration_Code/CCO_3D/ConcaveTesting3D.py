
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:46:48 2016

@author: cjaquet
"""

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from skimage.segmentation import random_walker
from scipy.interpolate import RegularGridInterpolator
import mpl_toolkits.mplot3d as a3
import matplotlib.cm as cm
import matplotlib as mpl
import time
#from pyamg import *



#intialization
center= np.array([14,14,14])
ext_radius_f = 10.
int_radius_f = 4.
slicing = 12

cx,cy,cz = int(center[0]), int(center[1]),int(center[2])
ext_radius = int(ext_radius_f)
int_radius =int(int_radius_f)

#initialization
im = np.zeros((cx+ext_radius*3, cy+ext_radius*3,cz+ext_radius*3)).astype('uint8')
markers = np.zeros((cx+ext_radius*3, cy+ext_radius*3, cz+ext_radius*3)).astype('uint8')
margin = int(np.ceil(ext_radius/10.))
z,y,x = np.ogrid[-ext_radius-margin: ext_radius+margin, -ext_radius-margin: ext_radius+margin, -ext_radius-margin:ext_radius+margin]  

#setting mask
index = x**2 + y**2 + z**2 <= ext_radius**2
im[cz-ext_radius-margin:cz+ext_radius+margin, cy-ext_radius-margin:cy+ext_radius+margin, cx-ext_radius-margin:cx+ext_radius+margin][index] = 1
z_i, y_i, x_i = np.ogrid[-int_radius: int_radius, -int_radius: int_radius,-int_radius: int_radius]
index_int = x_i**2 + y_i**2 + z_i**2<= int_radius**2
im[cz-int_radius:cz+int_radius, cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 0
#setting markers
#external
index = x**2 + y**2 +z**2 > ext_radius**2
markers[cz-ext_radius-margin:cz+ext_radius+margin, cy-ext_radius-margin:cy+ext_radius+margin, cx-ext_radius-margin:cx+ext_radius+margin][index] = 2
#internal
index_int = x_i**2 + y_i**2 +z_i**2 < int_radius**2
markers[cz-int_radius:cz+int_radius, cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 3
#random walker
deb = time.time()
result = random_walker(im, markers, copy =True, return_full_prob = True)  #, mode = 'cg_mg'
end = time.time()
print end - deb, "seconds for a ", im.shape[0]*im.shape[1]*im.shape[2] ,"pixel image" 
#print "rdm walker shape",result.shape
#print markers.shape

#make it a half sphere
result[1][center[2]:im.shape[2],:,:] = 0

wid=8
hei=8
fig = plt.figure(figsize=(wid*2, hei))
ax = fig.add_subplot(122,projection='3d')#a3.Axes3D()#pl.figure(figsize=(wid, hei))
#ax.imshow(im[20],20)
xx, yy = pl.ogrid[0:im.shape[1], 0:im.shape[2]]
#xx, yy = np.meshgrid(np.linspace(0,1,12), np.linspace(0,1,13))
# create vertices for a rotated mesh (3D rotation matrix)
X =  xx 
Y =  yy
Z =  slicing*np.ones(X.shape)
#fig = plt.figure(projection='3d')


#m = cm.ScalarMappable(cmap=cm.jet)
#m.set_array(markers[:,:,slicing])

#cmap for markers
colors = [(1.0,1.0,1.0)]
colors.extend([(0.0,0.0,1.)])
colors.extend([(1.0,0.72,0.06)])
colors.extend([(1.0,0.0,0.0)])
cmap = mpl.colors.ListedColormap(colors)

#cmap for random walker
N=50
colors = [(1.0,1.0,1.0)]
colors.extend(plt.cm.jet(np.linspace(0., 1., N)))
colors.extend([(1.0,1.0,1.0)])
cmap =mpl.colors.ListedColormap(colors) #plt.cm.jet

maping = ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(result[1][slicing,:,:]),shade=False)
#ax.colorbar()#maping,ax=ax#maping, shrink=0.5, aspect=5
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

ax1 = fig.add_subplot(121)
ax1.imshow(markers[:,:,slicing])

#cset = ax.contourf(X, Y, result[0][:,:,slicing], 100, zdir='z', offset=0.5, cmap=plt.cm.BrBG)
#plt.colorbar(cset)
plt.show()