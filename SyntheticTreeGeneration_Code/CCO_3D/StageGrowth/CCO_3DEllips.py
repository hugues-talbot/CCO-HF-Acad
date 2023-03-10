# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:35:19 2016

@author: jaquetc
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import Kamiya3D as kami
import NodeClass3D as nclass
import CCO_3DFunctions as cco_3df
import copy
from pylab import figure, gca, Line2D
import pylab as pl
import random
import os
import os.path
import cPickle as pickle
import sys  # Need to have acces to sys.stdout
import time
from multiprocessing import Pool, TimeoutError
import pickle
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, Ellipse, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d

############# Visualisation tools ####################

def plot_tree(tree, vol_descptr, name):
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
    plt.savefig(name+".png")
    plt.savefig(name+".pdf")
    plt.show()


    
##################################################################################
############# Karch algo : CCO ####################

timing = True
store_data = True
parallelized = False
half = True
cutof = True
cutof_val = 20

if timing:
    debut = time.time()
    print debut
if store_data:
    fd = open('./Results/CCO3D_newton_negpot_ellipsqd.txt','w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it    
    sys.stdout = fd # Now your file is used by print as destination 
    

if True:

    NTerm = 500
    if half:
        NTerm = 250
    seed = 42

    np.random.seed(seed)
    process_nb = 16
 
    #### Parameters to define: ##
    ## About tree
      
    Q_perf = 8.33e3
    N_term = NTerm
    InterTerm = 50
    Q_term = Q_perf / N_term

    P_drop = 1.33e4 - 8.38e3 # when Nterm = 4 000, the P_drop is 1.33e7 -7.98e6 #when =Nterm=250 :1.33e7 - 8.38e6
    viscosity = 3.6e-3 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 =3.6e-3 Pa.s = 3.6e-9 MPa.s 
    gamma = 3.0    
    # need to be carefull about unit: relativ between P_drop and viscosity
    # ex: 1.33e7  - 7.98e6 (kg mm-1 s-2)   and 3.6    (kg mm-1 s-2)
    # ex: 1.33e-2 - 7.98e-3 (MPa = N mm-2) and 3.6e-9 (MPa)
    # ex: 1.33e4  - 7.98e3 (Pa)            and 3.6e-3 (Pa.s)  
    N_con = 20

    # About  convexe perfusion surface : defines a disc surface 
    v_center = np.array([90,90,90])#np.array([14.,14., 14.])#np.array([80.,80.,80.])#np.array([80.,80.,80.])#
    #v_ext_radius =45#10.#50.#50#
    #v_int_radius =35#4.#15.#15#        

    r_ext = np.array([30,40,70])
    r_int = np.array([20,30,60])
    
    #in schreiner non convex cco: the total ellispoid volume is 48cm3
    #to use a similar volume in sphere we should take: r_ext = 45mm and r_int =35mm (so center = 55,55,55) 
    v_descptr = [v_center, r_ext, r_int]
    vperf = (4./3.) * np.pi*(r_ext[0]*r_ext[1]*r_ext[2] - r_int[0]*r_int[1]*r_int[2])
    filename = "potential_rext%i-%i-%i_rint%i-%i-%i_ellips" %(r_ext[0],r_ext[1],r_ext[2] ,r_int[0],r_int[1],r_int[2])
    if half:
        vperf = (4./3.) * np.pi*(r_ext[0]*r_ext[1]*r_ext[2] - r_int[0]*r_int[1]*r_int[2]) /2.
        filename = "potential_half_rext%i-%i-%i_rint%i-%i-%i_ellips" %(r_ext[0],r_ext[1],r_ext[2] ,r_int[0],r_int[1],r_int[2])
    filepath = "./Results/"+filename+".npy"
    if os.path.isfile(filepath):
        print "loading potential from numpy file %s" %filepath
        potential = np.load(filepath)
    else:
        potential = cco_3df.potential_ellipsoid(v_center, r_int,r_ext,half,cutof, cutof_val)
        np.save("./Results/"+filename, potential)
    
    
    #### initialization ##    
    store_cet = []
    tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity, vperf, np.power((r_ext[0]*r_ext[1]*r_ext[2] - r_int[0]*r_int[1]*r_int[2]),1/3.), potential, r_ext[2], v_center, r_ext[2], gamma)
    if half:      
        tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity, vperf, np.power((r_ext[0]*r_ext[1]*r_ext[2] - r_int[0]*r_int[1]*r_int[2])/2.,1/3.), potential, r_ext[2], v_center, r_ext[2], gamma) #v_ext_radius/2.?
     
    if half:
        root_position = np.array([v_center[0],v_center[1], v_center[2]+r_ext[2]])
    if cutof:
        #y_root = v_center[1]
        z_root = v_center[2] + r_ext[2] - cutof_val
        x_radius = cco_3df.length(np.sqrt(1.- ((z_root-v_center[2]) / float(r_ext[2]))**2 ) * r_ext[0]) #dont need it because y_root = v_center[1] - (y_root / r_ext[1])**2
        y_radius = cco_3df.length(np.sqrt(1.- ((z_root-v_center[2]) / float(r_ext[2]))**2 ) * r_ext[1])
        y_root = v_center[1] + 0.5 * y_radius
        x_root = v_center[0] + np.sqrt(1. - ((y_root - v_center[1]) /float(y_radius))**2) *x_radius
        root_position = np.array([x_root,y_root, z_root -1.6])
        
        
    print "root_position",root_position
    
    #when computing surface growth they accept 50mm segments diving below the surface as tolerance, with a 40mm external radius
    tol_schreiner = 50 #we use the same length
    schreiner_ext_radii = 40
    tol_dist_c = schreiner_ext_radii - np.sqrt(schreiner_ext_radii**2-(tol_schreiner / 2.)**2)  #and our external diameter  
    bellybutton = np.array([v_center[0]+r_ext[0],v_center[1], v_center[2]])
    bellybutton[0] = bellybutton[0]-tol_dist_c  
    w_tol = tree.get_w(bellybutton) 
    print "wtol when surface", w_tol
    
    if (tree.inside_perf_terr_exact(root_position)==True):
        print "root inside perf"
            
        root_node = nclass.Node(0,root_position, Q_perf, -1)
        root_node.set_child_0_index(1)
        root_node.set_label(0)
        tree.add_node(root_node)    
        
        #first segment end: randomly picked inside perfusion surface
        if InterTerm >0:          
                surface = True
                surface_tol = w_tol
        else: 
                surface = False
                surface_tol = 0.
        
        tree.update_length_factor()
        curvature_tolerance =0.05*tree.max_curv_rad
        first_node_position = tree.first_segmt_end(curvature_tolerance, True, surface_tol)
        first_node = nclass.Node(1, first_node_position, Q_perf,0)
        first_node.set_label(0)
        tree.add_node(first_node)
        tree.update_flow()
        tree.update_length_factor()
        tree.depthfirst_resistances(0)        
    
        counter = np.zeros(6)
        dead_end_counter = 0
        d_tresh_factor = 1
              

        while tree.get_k_term() < N_term:
            kterm = tree.get_k_term()
                
            if kterm < InterTerm:          
                surface = True
            else: 
                surface = False
                surface_tol = 0.
            
            success, new_child_location, d_tresh = tree.get_new_location(d_tresh_factor,surface)
            if (success == False):
                print "impossible to satisfy distance criteria", "d_tresh", d_tresh
                break       
            
            cet = [] 
            if parallelized == False:
                cet = [[] for i in range (N_con)]
            adding_location = False
            added_location = []
            
            dtype_r=[("convgce", int), ("volume", float), ("betas", float,(2)),("branching_location",float,(3)),("old_child_index", int)]
                    
            test_N_con_max = False
            
            # test closest neighbors
            neighbors = tree.find_neighbors(new_child_location, N_con)           
            args = [[tree, neighbors[i],new_child_location,curvature_tolerance,surface, surface_tol] for i in range (len(neighbors))]
    
            #process all neighbors connection test batch by batch
            if parallelized == True:
                while (len(cet) < N_con) and (len(cet) < len(neighbors)):           
                    end_lim = len(cet) + process_nb if (len(cet) + process_nb < len(neighbors)) else len(neighbors)            
                    pool = Pool(processes =  process_nb)               
                    res = pool.map(cco_3df.test_connection_list,args[len(cet): end_lim])    
                    cet = cet + [res[0][1:]]
                    #getting stat
                    pool.close()
            else:
                for n_index in range(len(neighbors)):
                    tree_copy = copy.deepcopy(tree)
                    res= tree_copy.test_connection(neighbors[n_index], new_child_location,curvature_tolerance,surface,surface_tol)                    
                    cet[n_index] = res[1:]
                    print "res[1:]", res[1:]                  
    
            
            cet_filtered = filter(None,cet)
            cet_values = np.array(cet_filtered, dtype_r)
            
            #if there is a candidate that converges
            if (np.sum(cet_values['convgce']) >= 1) or (np.sum(cet_values['convgce']) > 0 and tree.get_k_term() == 1):        
                cet_sel = cet_values[cet_values['convgce']>0]
                cet_sorted = np.sort(cet_sel, order = "volume")
                cet_final=cet_sorted[0]
                adding_location = True
                added_location.append(cet_final.tolist()[1:])
               
            if (adding_location): # optimal connection found!
                store_cet.append(filter(None,cet))
                opt = added_location[-1]
               
                if (tree.add_connection(opt[3], new_child_location, opt[2], opt[1])):
                    print "k termmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm is now ", tree.get_k_term()
                    kterm=tree.get_k_term()
                    d_tresh_factor = 1.
                    if kterm%10 == 0:
                        #plot_tree(tree, v_descptr, )
                        name ="./Results/InterTree_Nt%i_kt%i_s%i_ellip" %(NTerm,kterm,seed)
                        pickle.dump(tree, open(name + ".p", "wb"))
##
                    if kterm ==  30:
                        break
                else:
                    print "failed to add connection on tree"
            else:              
                print "ktemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm", tree.get_k_term()
                print "location doesn't provide an optimal connection, testing new location"
                dead_end_counter = dead_end_counter +1
                print "dead_end_counter = ", dead_end_counter
                if dead_end_counter == 10: #if too small doesn't produce homogeneous tree?
                    dead_end_counter = 0
                    d_tresh_factor = d_tresh_factor *0.85
                    print "dead end: decrease d_tresh of 20% to look for new location", d_tresh_factor

            #keep going until reach Nterm!
        
        
        name = "./Results/tree_Nt%i_kt%i_s%i_ellip" %(NTerm, tree.get_k_term(),seed)
        if half:
            name = "./Results/tree_Nt%i_kt%i_s%i_ellip" %(NTerm, tree.get_k_term(),seed)

        
        plot_tree(tree, v_descptr, name)
        pickle.dump(tree, open(name + ".p", "wb"))


if store_data:
    sys.stdout=old_stdout # here we restore the default behavior
    fd.close() # to not forget to close your file

if timing:
    fin =time.time()
    print "duration = ", fin-debut, "secondes"
    

