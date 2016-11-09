# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:13:38 2016

@author: jaquetc
"""
import numpy as np
import NodeClass3D as nclass
import CCO_3DFunctions as cco_3df
import copy
import random
import os
import os.path
import sys  # Need to have acces to sys.stdout
import cPickle as pickle
import pickle
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d
from pylab import figure, gca, Line2D
import pylab as pl
############# Visualisation tools ####################

def plot_tree(tree, name, half):
    wid=16 # Must be proportional to x and y limits below
    hei=16
    ax = a3.Axes3D(pl.figure(figsize=(wid, hei)))
    
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    center = tree.center
    r_ext = tree.real_final_radius
    if half:
        r_ext = tree.real_final_radius#*2. if root at border no need to multiply by 2
    r_int = tree.max_curv_rad
    x = center[0]+ r_ext * np.outer(np.cos(u), np.sin(v))
    y = center[1]+ r_ext * np.outer(np.sin(u), np.sin(v))
    z = center[2]+ r_ext * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1) #,
    
    x = center[0]+ r_int * np.outer(np.cos(u), np.sin(v))
    y = center[1]+ r_int * np.outer(np.sin(u), np.sin(v))
    z = center[2]+ r_int * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1)

    p = Circle((center[0:2]), r_ext, facecolor = 'b', alpha=0.1) #Add a circle in the yz plane
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z = center[1], zdir = 'x')
    
    im = tree.w_pot
    xx, yy = pl.ogrid[0:im.shape[1], 0:im.shape[2]]
    #xx, yy = np.meshgrid(np.linspace(0,1,12), np.linspace(0,1,13))
    # create vertices for a rotated mesh (3D rotation matrix)
    X =  xx 
    Y =  yy
    Z =  (center[2]-1)*np.ones(X.shape)

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
    
    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]-1,:,:].transpose()),shade=False,alpha=0.2)
    Z =  (center[2]+1)*np.ones(X.shape)
    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]+1,:,:].transpose()),shade=False,alpha=0.2)
#    Z =  (center[2]+10)*np.ones(X.shape)
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]+10,:,:]),shade=False,alpha=0.2) 
#    Z =  (center[2]+20)*np.ones(X.shape)
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]+20,:,:]),shade=False,alpha=0.2)
#    Z =  (center[2]+20)*np.ones(X.shape)
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]+20,:,:]),shade=False,alpha=0.2)
#    
    
    
    #setting figure so that we get linewidth in data unit
    bound = r_ext*1.5
    ax.set_xlim(center[0] - bound,center[0]+ bound)#using same bound diff for all axis 
    ax.set_ylim(center[1] - bound,center[1]+ bound)#--> proportional to the figure size
    ax.set_zlim(center[2] - bound,center[2]+ bound)
    ax.axes.set_aspect('equal')
    x1,x2,y1,y2 = ax.axis()
    print "axis",x1,x2,y1,y2
    yrange =   y2 - y1 #don't know if necessary to calculate pointlinewid_factor
    # line is in points: 72 points per inch
    point_hei=hei*72 
    # For the calculation below, you have to adjust width by 0.8
    # because the top and bottom 10% of the figure are labels & axis
    pointlinewid_factor = point_hei * 0.8 /(0.5*yrange) # corresponding width in pts ( /yrange ?)
    
    inv_length_fac = 1. / tree.length_factor
    for sgmt in tree.nodes:
        if (sgmt.parent() >= 0):
            distal = sgmt.coord
            proximal = tree.get_node(sgmt.parent()).coord
            radius = tree.get_radius(sgmt.index) *inv_length_fac
            print "radius", radius
            verts = [zip([distal[0], proximal[0]],[distal[1], proximal[1]],[distal[2], proximal[2]])]
            tri = a3.art3d.Poly3DCollection(verts)
            tri.set_color('k')
           

            tri.set_linewidth(radius*2.*pointlinewid_factor/(r_ext*3.) )
            ax.add_collection3d(tri)
            

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.grid(False)
    #plt.savefig(name+".png")
    #plt.savefig(name+".pdf")
    print "y_range", yrange
    plt.show()

nterm = 2000
seed=42
kterm=290

#tree = pickle.load( open( "./Results/InterTree_Nt%i_kt%i_s%i_half_nr.p"% (nterm,kterm, seed), "rb" ) )
tree = pickle.load( open( "./Results/InterTree_Nt250_kt20_s42_realcutof.p", "rb" ) )


plot_tree(tree,"figname", True)