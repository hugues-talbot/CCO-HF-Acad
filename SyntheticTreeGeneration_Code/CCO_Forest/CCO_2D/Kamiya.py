# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:34:54 2016

@author: cjaquet

 Kamiya algo: single bifurcation optimization 
based on minimal intravascular volume and physiological constraints


"""


import numpy as np
from scipy import optimize

### loop process
# given a bifurcation starting_point (convexe average of the 3 points), and estimated radii (using tree existing segments) 
# calculate segments radii with Poiseuille and minimal vascular volume (non linear system)
# use them to calculate new bifurcation point, which serves as a new input to calculate again radii, and so on... 
# if position converges and total tree volume decreases: solved 


#ci contains coordinate xi and y1 of the point
#f contains f0,f1,f2
#returns position of starting point
def starting_point(c0,c1,c2,f):
    if (len(c0) != 2) or (len(c1) != 2) or (len(c2) != 2) or (len(f)!=3):
        print "error in starting point inputs"
        return 0.,0.
    x_coord = (c0[0]*f[0] + c1[0]*f[1] + c2[0]*f[2]) / (2* f[0])
    y_coord = (c0[1]*f[0] + c1[1]*f[1] + c2[1]*f[2]) / (2* f[0])
    return x_coord, y_coord
    
def calculate_segment_lengths(c0,c1,c2,x,y, length_factor):
    coords = np.array([x,y])
    l = np.ones(3)
    l[0] = rescaled_length(c0 - coords, length_factor)
    l[1] = rescaled_length(c1 - coords, length_factor)
    l[2] = rescaled_length(c2 - coords, length_factor)
    return l

#dp are the pressure loss: P1-P0 and P2-P0 
#the r_ori is an estimated one to serve the iteration start
def calculate_dp_from_Poiseuille(f,l,r_ori, visc):
    kappa = 8 * visc / np.pi    
    r0_ori4 = r_ori[0]*r_ori[0]*r_ori[0]*r_ori[0]
    r1_ori4 = r_ori[1]*r_ori[1]*r_ori[1]*r_ori[1]
    r2_ori4 = r_ori[2]*r_ori[2]*r_ori[2]*r_ori[2]
    dp1 = kappa * (f[0]*l[0]/r0_ori4 + f[1]*l[1]/r1_ori4)
    dp2 = kappa * (f[0]*l[0]/r0_ori4 + f[2]*l[2]/r2_ori4)
    return dp1,dp2 
       
#sqrd_r contains r0**2, r1**2, r2**2
def calculate_new_bif_coords(c0,c1,c2,f,l,sqrd_r):
    new_x = ( c0[0]*sqrd_r[0]/l[0] + c1[0]*sqrd_r[1]/l[1] + c2[0]*sqrd_r[2]/l[2] ) /  ( sqrd_r[0]/l[0] + sqrd_r[1]/l[1] + sqrd_r[2]/l[2] )
    new_y = ( c0[1]*sqrd_r[0]/l[0] + c1[1]*sqrd_r[1]/l[1] + c2[1]*sqrd_r[2]/l[2] ) /  ( sqrd_r[0]/l[0] + sqrd_r[1]/l[1] + sqrd_r[2]/l[2] )
    return new_x, new_y
    
    
def non_linear_solver(f,l,k,dp1,dp2,r_ori,gamma):
    def func(z):
        r1=z[0]
        r2=z[1]
        #print "z", z
        zeta= 3+gamma
        r14= r1*r1
        r16 = r1**(zeta/2.) #r14*r1
        r24= r2*r2
        r26 = r2**(zeta/2.)#r24*r2
        ans = ( dp1/k*r14*(np.power(f[0]*(r16/f[1]+r26/f[2]),4./zeta )) - f[0]*l[0]*r14 - f[1]*l[1]*np.power(f[0]*(r16/f[1] +r26/f[2]), 4./zeta),
                dp2/k*r24*(np.power(f[0]*(r16/f[1]+r26/f[2]),4./zeta )) - f[0]*l[0]*r24 - f[2]*l[2]*np.power(f[0]*(r16/f[1] +r26/f[2]), 4./zeta))         
        return(ans)
    #print"r_ori",  r_ori
    dict_f={}
    dict_f['solver']="root"
    dict_f['method']="hybr"
    dict_f['disp']= False
    #optimize.optimize.
    #optimize.show_options(solver=None, method = None, disp =  False)
    sol=optimize.root(func,r_ori, jac = False )# options = {'solver' : 'root', 'method' :'hybr','disp':False}
    #print "sol message", sol.message
    #print "Value =", sol.fun
    #print "sol", sol.x
    
    return sol.success, sol.x
    
#l contains l0, l1, l2
def calculate_squared_radii(f, l, dp1, dp2, r_ori, visc,gamma):
    k=8*visc/np.pi
    success, r_1_2 = non_linear_solver(f,l,k,dp1,dp2, r_ori,gamma)
    r1 = r_1_2[0]
    r2 = r_1_2[1]
    if (success == True):
        r0 = np.power(f[0]*(r1**gamma/f[1] + r2**gamma/f[2]) , 1/gamma)  
        return r0,r1,r2
    else:
        return 0.,0.,0.
        
def single_bif_volume(r, l):
    return np.sum(np.pi*(r**2)*l, axis = 0)
    
#input r : r[0] and r[1] are the radii before connection added (r0 = r1 = r1 = radius of original segment before connection added) 
#        this estimated r is used to calculate dp, then is updated in the calculate_radii function     
def kamiya_loop_r2(x_ini,y_ini,c0,c1,c2,f, r, length_factor, visc,gamma):
    l = calculate_segment_lengths(c0,c1,c2,x_ini,y_ini, length_factor)
    dp1, dp2 = calculate_dp_from_Poiseuille(f,l,r, visc)
    r[0], r[1], r[2] = calculate_squared_radii(f,l,dp1, dp2, np.power(r,2)[1:3], visc,gamma)
    if (r[0] == 0.):
        return False, x_ini, y_ini, np.sqrt(r), l
    else:        
        x_new,y_new = calculate_new_bif_coords(c0,c1,c2,f,l,r)
        return True, x_new, y_new, np.sqrt(r), l 
        
#when calculating segment length, need to consider a factor 
#because current coordinates are scaled from a bigger perfusion territory
def rescaled_length(vector, factor):
    return np.sqrt(np.sum(vector**2))*factor

def calculate_seg_lengths(c0,c1,c2,x,y, length_factor):
    coords = np.array([x,y])
    l = np.zeros(3)
    l[0] = rescaled_length(c0 - coords, length_factor)
    l[1] = rescaled_length(c1 - coords, length_factor)
    l[2] = rescaled_length(c2 - coords, length_factor)
    return l


# r =r0,r1,r2 in mm
# f = flow in mm3/s
# c are coordinates of each segment end
def kamiya_single_bif_opt(r, f, c0, c1, c2, iter_max, tolerance, store_whole_results, length_factor):
    x, y = starting_point(c0,c1,c2,f)
    print "starting point coord ", x, y
    print "starting radii ", r
    it = 0
    gdt = 1
    l_ori = calculate_segment_lengths(c0,c1,c2,x,y, length_factor)
    volume = single_bif_volume(r, l_ori)
    storage = []
    
    while it < iter_max and gdt > tolerance:     
        cvge, x_c, y_c, r_c, l = kamiya_loop_r2(x, y, c0, c1, c2, f, r, length_factor)
        if (cvge):
            vol_c = single_bif_volume(r_c,l) 
            gdt = np.sqrt((x_c - x)**2 + (y_c - y)**2) * length_factor      
            if (gdt < tolerance) and (vol_c < volume):              
                print "\n convergence found at iteration ", it
                #print "position gdt ", gdt
                #print "volume gdt ", volume - vol_c
                #print "optimal volume ", vol_c
                print  "optimal position ", x_c, y_c
                print "optimal radii ", np.sqrt(r_c)
                storage.append([x_c, y_c, r_c])
                return True, storage #x_c, y_c, r
            else:
                it = it + 1
                x , y , volume, r = x_c, y_c, vol_c, r_c
                if store_whole_results:
                    storage.append([x, y, r]) #### need it if used in Kamyia testing               
                #print " position gdt ", gdt
                #print "volume gdt ", volume - vol_c
                #print "volume ", vol_c
                #print "iterating with new position ", it
        else:
            print "\n no convergence"
            return False, storage # x, y, r
    
    print " \n max iteration reached"
    return False, storage # x, y, r
            
