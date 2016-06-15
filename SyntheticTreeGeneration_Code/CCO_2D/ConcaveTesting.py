# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:46:48 2016

@author: cjaquet
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.segmentation import random_walker

#### testing random walker to create potential between two surfaces #####

##image definition
if False:
    ext_radius = 200
    int_radius = 50
    a = np.zeros((512, 512)).astype('uint8')
    cx, cy = 256, 256 # The center of circle
    y, x = np.ogrid[-ext_radius: ext_radius, -ext_radius: ext_radius]
    index = x**2 + y**2 <= ext_radius**2
    a[cy-ext_radius:cy+ext_radius, cx-ext_radius:cx+ext_radius][index] = 1
    
    
    y_i, x_i = np.ogrid[-int_radius: int_radius, -int_radius: int_radius]
    index_int = x_i**2 + y_i**2 <= int_radius**2
    a[cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 0
      
    #markers definition
    markers = np.zeros((512, 512)).astype('uint8')
    index = np.logical_and(x**2 + y**2 >= ext_radius**2 * 0.98, x**2 + y**2 <= ext_radius**2)
    markers[cy-ext_radius:cy+ext_radius, cx-ext_radius:cx+ext_radius][index] = 2
    index_int = np.logical_and(x_i**2 + y_i**2 >= int_radius**2 * 0.90, x_i**2 + y_i**2 <= int_radius**2)
    markers[cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 3

    result = random_walker(a, markers, copy =True, return_full_prob = True)
    
    plt.subplot(1,4, 1)
    plt.imshow(a)
    plt.subplot(1,4,2)
    plt.imshow(markers)
    plt.subplot(1,4,3)
    plt.imshow(result[0])
    plt.subplot(1,4,4)
    plt.imshow(result[1])
    plt.show()
    
