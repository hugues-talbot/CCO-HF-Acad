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
from matplotlib.patches import Circle, Ellipse, PathPatch
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
    
    #ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]-1,:,:].transpose()),shade=False,alpha=0.2)
    Z =  (center[2]+1)*np.ones(X.shape)
    #ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[center[2]+1,:,:].transpose()),shade=False,alpha=0.2)
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
    
def plot_tree_ellips(tree,vol_descptr, name,cutof_val):
    wid=16 # Must be proportional to x and y limits below
    hei=16
    ax = a3.Axes3D(pl.figure(figsize=(wid, hei)))
    
    rx,ry,rz =vol_descptr[1][0], vol_descptr[1][1], vol_descptr[1][2]
    rix,riy,riz =vol_descptr[2][0], vol_descptr[2][1], vol_descptr[2][2]

    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = vol_descptr[0][0] + rx * np.outer(np.cos(u), np.sin(v))
    y = vol_descptr[0][1] + ry * np.outer(np.sin(u), np.sin(v))
    z = vol_descptr[0][2] + rz * np.outer(np.ones_like(u), np.cos(v))
    
    xi = vol_descptr[0][0] + rix * np.outer(np.cos(u), np.sin(v))
    yi = vol_descptr[0][1] + riy * np.outer(np.sin(u), np.sin(v))
    zi = vol_descptr[0][2] + riz * np.outer(np.ones_like(u), np.cos(v))

   
    ax.plot_surface(x, y, z,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1) #,
    ax.plot_surface(xi, yi, zi,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1)
    
    #Add an ellips in the yz plane
    e = Ellipse(xy=(vol_descptr[0][1] ,vol_descptr[0][2] ), width= ry*2, height=rz*2,facecolor='b',alpha = 0.1)
    ax.add_patch(e)
    art3d.pathpatch_2d_to_3d(e, z = vol_descptr[0][0], zdir = 'x')
    
    #Add an ellips in the xy plane
    z_cut = vol_descptr[0][2] + rz - cutof_val
    x_cut = cco_3df.length(np.sqrt(1.- ((z_cut-vol_descptr[0][2]) / float(rz))**2 ) * rx) 
    y_cut = cco_3df.length(np.sqrt(1.- ((z_cut-vol_descptr[0][2]) / float(rz))**2 ) * ry) 
    e = Ellipse(xy=(vol_descptr[0][0] ,vol_descptr[0][1] ), width= x_cut*2, height=y_cut*2,facecolor='b',alpha = 0.1)
    ax.add_patch(e)
    art3d.pathpatch_2d_to_3d(e, z = vol_descptr[0][2] + rz -20, zdir = 'z')
    
    
#    im = tree.w_pot
#    xx, yy = pl.ogrid[0:im.shape[1], 0:im.shape[2]]
#    #xx, yy = np.meshgrid(np.linspace(0,1,12), np.linspace(0,1,13))
#    # create vertices for a rotated mesh (3D rotation matrix)
#    X =  xx 
#    Y =  yy
#    Z =  (vol_descptr[0][2]-1)*np.ones(X.shape)
#
#    #cmap for markers
#    colors = [(1.0,1.0,1.0)]
#    colors.extend([(0.0,0.0,1.)])
#    colors.extend([(1.0,0.72,0.06)])
#    colors.extend([(1.0,0.0,0.0)])
#    cmap = mpl.colors.ListedColormap(colors)
#    
#    #cmap for random walker
#    N=50
#    colors = [(1.0,1.0,1.0)]
#    colors.extend(plt.cm.jet(np.linspace(0., 1., N)))
#    colors.extend([(1.0,1.0,1.0)])
#    cmap =mpl.colors.ListedColormap(colors) #plt.cm.jet
#    
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[vol_descptr[0][2]-1,:,:].transpose()),shade=False,alpha=0.2)
#    Z =  (vol_descptr[0][2]+1)*np.ones(X.shape)
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(im[vol_descptr[0][2]+1,:,:].transpose()),shade=False,alpha=0.2)
##
#    
    #setting figure so that we get linewidth in data unit
    bound = vol_descptr[1][2]*1.5
    ax.set_xlim(vol_descptr[0][0] - bound,vol_descptr[0][0]+ bound)#using same bound diff for all axis 
    ax.set_ylim(vol_descptr[0][1] - bound,vol_descptr[0][1]+ bound)#--> proportional to the figure size
    ax.set_zlim(vol_descptr[0][2] - bound,vol_descptr[0][2]+ bound)
    ax.axes.set_aspect('equal')
    x1,x2,y1,y2 = ax.axis()
    print "axis",x1,x2,y1,y2
    yrange =   y2 - y1 #don't know if necessary to calculate pointlinewid_factor
    # line is in points: 72 points per inch
    point_hei=hei*72 
    # For the calculation below, you have to adjust width by 0.8
    # because the top and bottom 10% of the figure are labels & axis
    pointlinewid_factor = point_hei * 0.8 #/yrange # corresponding width in pts ( /yrange ?)
    
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
          
            tri.set_linewidth(radius*2.*pointlinewid_factor /(vol_descptr[1][2]*3.))
            ax.add_collection3d(tri)
            

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.grid(False)
    #plt.savefig(name+".png")
    #plt.savefig(name+".pdf")
    plt.show()

nterm = 2000
seed=42
kterm=290

#tree = pickle.load( open( "./Results/InterTree_Nt%i_kt%i_s%i_half_nr.p"% (nterm,kterm, seed), "rb" ) )
tree = pickle.load( open( "./Results/tree_Nt250_kt15_s42_ellip.p", "rb" ) )
v_center = np.array([90,90,90])#np.array([14.,14., 14.])#np.array([80.,80.,80.])#np.array([80.,80.,80.])#
#v_ext_radius =45#10.#50.#50#
#v_int_radius =35#4.#15.#15#        

r_ext = np.array([30,40,70])
r_int = np.array([20,30,60])
cutof = 20
#in schreiner non convex cco: the total ellispoid volume is 48cm3
#to use a similar volume in sphere we should take: r_ext = 45mm and r_int =35mm (so center = 55,55,55) 
v_descptr = [v_center, r_ext, r_int]
plot_tree_ellips(tree,v_descptr, "fig_elli", 20)
#plot_tree(tree,"figname", True)