# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 13:48:17 2016

@author: cjaquet
"""
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
import Kamiya_3D as kami
import NodeClass3D as nclass
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from pylab import figure, gca, Line2D
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt

##################################Visualisation tool###########################
def plot_trees(trees, area_descptr):
    #fig = plt.figure(figsize=(8,8))
    #ax = fig.add_subplot(111, projection='3d')#fig.gca(projection='3d')
    #ax = Axes3D(fig)    
    wid=8
    hei=8
    ax = a3.Axes3D(pl.figure(figsize=(wid, hei)))
    colors = ['b', 'r']
    ind = 0
    for tree in trees:      
        for sgmt in tree.nodes:
            if (sgmt.parent() >= 0):
                distal = sgmt.coord
                print "distal", distal
                proximal = tree.get_node(sgmt.parent()).coord
                radius = tree.get_radius(sgmt.index) 
                #print radius
                #ax.add_line(Line2D([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = colors[ind], alpha = 0.4))
                
                #ax.plot3D([distal[0], proximal[0]],[distal[1], proximal[1]], [distal[2], proximal[2]], linewidth = radius, zdir='z', color = colors[ind], alpha = 0.4)
                
                verts = [zip([distal[0], proximal[0]],[distal[1], proximal[1]],[distal[2], proximal[2]])]
                tri = a3.art3d.Poly3DCollection(verts)
                tri.set_color('b')
                tri.set_linewidth(radius*72*hei/(area_descptr[1]*3.))#*0.8
                tri.set_alpha(0.5)
                #tri.set_norm
                ax.add_collection3d(tri)
                
                #plot3D(x,y,z, zdir = 'z')
                #ax.contour(x,y,z,stride = 0.96, zdir = 'z')
        ind = ind +1   
    ax.set_xlim(-25,125)#area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5
    ax.set_ylim(-25,125)#area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5
    ax.set_zlim(-25,125)#area_descptr[0][2] - area_descptr[1]*1.5,area_descptr[0][2]+ area_descptr[1]*1.5
    plt.show()
    #Axes3D.plot()
    #return fig


###############################################################################


Q_perf = 8.33e3
N_term = 250.		
Q_term = Q_perf / N_term
P_drop = 1.33e7 - 8.38e6
viscosity = 3.6 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 (check works with radius and length in mm)

######## Testing Kamiya's algorithm #################

##initialization

# tree 
tree = nclass.Tree([], 24, Q_perf, P_drop, viscosity)

z_val = 0.
#source point 		
area_center = np.array([100.,50.,z_val])
area_radius = 50
area_descptr = [area_center, area_radius]
area = np.pi*area_radius**2
root_position = np.array([100,100,z_val])#100,100
root_node = nclass.Node(0,root_position, Q_perf, -1)
tree.add_node(root_node)


#first segment end
#generate random position > dist_crit

first_node_position = np.array([120,20,z_val+10])#100,20#first_segmt_end(area)# question: should it be over the d_tresh?
print "first segment end point", first_node_position
first_node = nclass.Node(1, first_node_position, Q_perf,0)
#first_node.set_betas(np.array([0.6,0.8]))
tree.add_node(first_node)


test = np.array([80,20,z_val-10])#75,40# np.array([75,40])#get_new_location(tree, area_descptr, N_term)
print "new location", test
#second_node = nclass.Node(2, test, Q_perf, 1)
#tree.add_node(second_node)
radius_ori = tree.get_radius(1)
r = np.ones(3)*radius_ori
print "initial radii", r
f = np.array([Q_perf, Q_perf*1./2., Q_perf*1./2.])#np.array([5666,2833,2833])#Q_perf, Q_perf*1./4., Q_perf*1./4.
print "initial flow", f
c0 = root_position
c1 = first_node_position
c2 = test
print "initial positions", c0,c1,c2
iter_max = 100
tolerance = 0.01
convergence, res = kami.kamiya_single_bif_opt(r, f, c0, c1, c2, iter_max, tolerance, True)
gamma = 3.
ind = 0
trees = []
for i in res:
    if ind == 0 : #or ind == len(res)-1
        tre = nclass.Tree([], 24, Q_perf, P_drop, viscosity)
        
        tre.add_node(root_node)
        
        bif_node = nclass.Node(1, np.array([i[0], i[1], i[2]]), Q_perf,0)
        betas = nclass.calculate_betas(i[3][1]/i[3][2],gamma)
        bif_node.set_betas(betas)
        bif_node.set_child_0_index(2)
        bif_node.set_child_1_index(3)
        tre.add_node(bif_node)
        
        old_child = nclass.Node(2, first_node_position, Q_perf/2.,1)
        tre.add_node(old_child)
        
        new_child = nclass.Node(3, test, Q_perf/2.,1)
        tre.add_node(new_child)
        trees.append(tre)
        #cco.plot_tree(tre, area_descptr, str(ind))
    ind = ind + 1
    
plot_trees(trees, area_descptr)