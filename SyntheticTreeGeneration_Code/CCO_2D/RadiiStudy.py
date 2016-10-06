# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:01:15 2016

@author: jaquetc
"""

import numpy as np
import matplotlib.pyplot as plt

if False:
    input_file = open("./Results/tree_Nt4000_s42_33all.txt",'r')
    lines = [float(line.rstrip('\n')) for line in input_file]
    #radii_list = input_file.readlines()
    radii = np.array(lines)
    print "len radii", radii.shape
    print "max diam ",np.max(radii)*2.
    plt.hist(radii*2., 100, color ='b',alpha=0.5)
    plt.title("Histogram of vessel diameter (mm)")
    plt.savefig("./Results/tree_Nt4000_s42_33all.png")
    plt.show()

if False:
    #one plot mean over 10 models
    all_means=np.zeros((10,40))
    overall_mean = np.zeros((40))
    overall_std = np.zeros((40))
    seed_start = 33
    seed_end = 43
    seed_index=np.arange(seed_start,seed_end,1)
    max_level= 0
    for i in range(10):      
        radii_along_bif = np.load("./Results/radii_Nt250_s%i_M255_PS.npy" %(seed_index[i]))
        means = np.zeros(len(radii_along_bif))
        #stds = np.zeros(len(radii_along_bif))
        level = (len(radii_along_bif))
        if max_level < level:
            max_level = level
        for j in range (len(radii_along_bif)):
            diam = np.array(radii_along_bif[j])*2.
            means[j] = np.mean(diam)
            #stds[i] = np.std(diam)
        all_means[i][0:level] = means
    overall_std = np.std(all_means, axis = 0) 
    overall_mean = np.mean(all_means, axis=0)
    
    #cleaning the last bifurcation levels (not all trees reach this level so need to count only some trees)
    nb_c = max_level - 20
    mean_over_0=[[] for k in range(nb_c)]
    concerned_levels  = np.arange(max_level-nb_c,max_level)
    for i in range(nb_c):      
        for j in range(10):
            if (all_means[j][concerned_levels[i]]> 0.0):
                mean_over_0[i].append(all_means[j][concerned_levels[i]])                                
            else:
                print "miss"              
    for i in range(nb_c):
        overall_mean[concerned_levels[i]] = np.mean(mean_over_0[i])
        overall_std[concerned_levels[i]] = np.std(mean_over_0[i])
            
    bif_levels = np.arange(max_level)    
    plt.errorbar(bif_levels,overall_mean[0:max_level], overall_std[0:max_level])
    plt.title("Diameter along bifurcation levels (mm) in 10 generated tree with 250Nterm")
    plt.xlabel("Bifurcation level")
    plt.ylabel("Diameter (mm)")
    plt.xlim(-1, max_level +1 )
    plt.ylim(0,5)
    plt.savefig("./Results/DiamAlongBif_Nt250_over10trees_M211_PS_seeds%ito%i.png" %(seed_start,seed_end))
    plt.show()
    
if False:
    #one plot of one model
    seed_index= 35
    radii_along_bif = np.load("./Results/radii_Nt250_s%i_M301_PK.npy" %(seed_index))
    means = np.zeros(len(radii_along_bif))
    stds = np.zeros(len(radii_along_bif))
    for j in range (len(radii_along_bif)):
        diam = np.array(radii_along_bif[j])*2.
        means[j] = np.mean(diam)
        stds[j] = np.std(diam)
            
    bif_levels = np.arange(len(radii_along_bif))    
    plt.errorbar(bif_levels,means, stds, fmt='o')
    plt.title("Diameter along bifurcation levels (mm) in a single generated tree with 250Nterm")
    plt.xlabel("Bifurcation level")
    plt.ylabel("Diameter (mm)")
    plt.xlim(-1, max_level +1 )
    plt.ylim(0,3)
    plt.savefig("./Results/DiamAlongBif_Nt250_over10trees_M255_PK_seeds%ionly.png" %(seed_index))
    plt.show()

if True:
    # several plots with mean over 10 models, gamma varying
    seed_start = 33
    seed_end = 43
    
    #gamma 211
    all_means=np.zeros((10,40))
    overall_mean = np.zeros((40))
    overall_std = np.zeros((40))
    seed_index=np.arange(seed_start,seed_end,1)
    max_level= 0
    for i in range(10):      
        radii_along_bif = np.load("./Results/radii_Nt250_s%i_M211.npy" %(seed_index[i]))
        means = np.zeros(len(radii_along_bif))
        #stds = np.zeros(len(radii_along_bif))
        level = (len(radii_along_bif))
        if max_level < level:
            max_level = level
        for j in range (len(radii_along_bif)):
            diam = np.array(radii_along_bif[j])*2.
            means[j] = np.mean(diam)
            #stds[i] = np.std(diam)
        all_means[i][0:level] = means
    overall_std = np.std(all_means, axis = 0) 
    overall_mean = np.mean(all_means, axis=0)
    
    #cleaning the last bifurcation levels (not all trees reach this level so need to count only some trees)
    nb_c = max_level - 20
    mean_over_0=[[] for k in range(nb_c)]
    concerned_levels  = np.arange(max_level-nb_c,max_level)
    for i in range(nb_c):      
        for j in range(10):
            if (all_means[j][concerned_levels[i]]> 0.0):
                mean_over_0[i].append(all_means[j][concerned_levels[i]])                                
            else:
                print "miss"              
    for i in range(nb_c):
        overall_mean[concerned_levels[i]] = np.mean(mean_over_0[i])
        overall_std[concerned_levels[i]] = np.std(mean_over_0[i])
            
    bif_levels = np.arange(max_level)    
    plt.errorbar(bif_levels,overall_mean[0:max_level], overall_std[0:max_level], color = 'b', label = "gamma 2.10")
    
    #gamma 255
    all_means=np.zeros((10,40))
    overall_mean = np.zeros((40))
    overall_std = np.zeros((40))
    seed_index=np.arange(seed_start,seed_end,1)
    max_level= 0
    for i in range(10):      
        radii_along_bif = np.load("./Results/radii_Nt250_s%i_M255.npy" %(seed_index[i]))
        means = np.zeros(len(radii_along_bif))
        #stds = np.zeros(len(radii_along_bif))
        level = (len(radii_along_bif))
        if max_level < level:
            max_level = level
        for j in range (len(radii_along_bif)):
            diam = np.array(radii_along_bif[j])*2.
            means[j] = np.mean(diam)
            #stds[i] = np.std(diam)
        all_means[i][0:level] = means
    overall_std = np.std(all_means, axis = 0) 
    overall_mean = np.mean(all_means, axis=0)
    
    #cleaning the last bifurcation levels (not all trees reach this level so need to count only some trees)
    nb_c = max_level - 20
    mean_over_0=[[] for k in range(nb_c)]
    concerned_levels  = np.arange(max_level-nb_c,max_level)
    for i in range(nb_c):      
        for j in range(10):
            if (all_means[j][concerned_levels[i]]> 0.0):
                mean_over_0[i].append(all_means[j][concerned_levels[i]])                                
            else:
                print "miss"              
    for i in range(nb_c):
        overall_mean[concerned_levels[i]] = np.mean(mean_over_0[i])
        overall_std[concerned_levels[i]] = np.std(mean_over_0[i])
            
    bif_levels = np.arange(max_level)    
    plt.errorbar(bif_levels,overall_mean[0:max_level], overall_std[0:max_level], color = 'r', label="gamma 2.55")
    
    #gamma 300   
    all_means=np.zeros((10,40))
    overall_mean = np.zeros((40))
    overall_std = np.zeros((40))
    seed_index=np.arange(seed_start,seed_end,1)
    max_level= 0
    for i in range(10):      
        radii_along_bif = np.load("./Results/radii_Nt250_s%i_33.npy" %(seed_index[i]))
        means = np.zeros(len(radii_along_bif))
        #stds = np.zeros(len(radii_along_bif))
        level = (len(radii_along_bif))
        if max_level < level:
            max_level = level
        for j in range (len(radii_along_bif)):
            diam = np.array(radii_along_bif[j])*2.
            means[j] = np.mean(diam)
            #stds[i] = np.std(diam)
        all_means[i][0:level] = means
    overall_std = np.std(all_means, axis = 0) 
    overall_mean = np.mean(all_means, axis=0)
    
    #cleaning the last bifurcation levels (not all trees reach this level so need to count only some trees)
    nb_c = max_level - 20
    mean_over_0=[[] for k in range(nb_c)]
    concerned_levels  = np.arange(max_level-nb_c,max_level)
    for i in range(nb_c):      
        for j in range(10):
            if (all_means[j][concerned_levels[i]]> 0.0):
                mean_over_0[i].append(all_means[j][concerned_levels[i]])                                
            else:
                print "miss"              
    for i in range(nb_c):
        overall_mean[concerned_levels[i]] = np.mean(mean_over_0[i])
        overall_std[concerned_levels[i]] = np.std(mean_over_0[i])
            
    bif_levels = np.arange(max_level)    
    plt.errorbar(bif_levels,overall_mean[0:max_level], overall_std[0:max_level], color = 'g', label = "gamma 3.00")
    
    
    
    plt.title("Diameter along bifurcation levels (mm) in 10 generated tree (250Nterm), \n with different Murray coefficient gamma")
    plt.xlabel("Bifurcation level")
    plt.ylabel("Diameter (mm)")
    plt.xlim(-1, max_level +1 )
    plt.ylim(0,5)
    plt.legend()
    plt.savefig("./Results/DiamAlongBif_Nt250_seeds%ito%i_gammacomparison.png" %(seed_start,seed_end))
    plt.show()