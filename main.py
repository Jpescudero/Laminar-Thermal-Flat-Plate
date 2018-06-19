import numpy as np
import matplotlib.pyplot as plt


from Flat_Plate import *
from FreeStream import *
from Boundary_Layer import *
from Plotting import *

def main(plot=False):

	## Problem Specification
    ## =================================================================== 

    # Suction veocity (m/s)
    vw = 0 

    # Lenght of the plate (m)  
    L  = 1      

    # Plate Surface Temperature (K)  
    Tp  = 300 

	# Creation of a flat Plate: (Length, Wall Temperature,suction velocity)
    fp = flat_plate(L,Tp,vw)

    # Plot Vectors
    Re      = []
    Cd_l    = []
    Cd_t    = []
    Cd_l_t  = []
    xt      = []

    
    for ue in np.logspace(3,3,100):

        print(ue)

    	# Creation of a Free stream: (Free Stream Velocity, Characteristic Length)
        fs = freestream(ue,L)
    	
        Re.append(fs.Re)

        print(fs.Re)

    	# Estimation of boundary layer height
        delta_est = L/(fs.Re)**0.5

    	# Grid Flat Plate Generation
        n_x = 500
        n_y = 200

        # Creation of the flat plate grid
        fp.grid(n_x,n_y,delta_est,ue,fs.nu)

    	#Creation of a laminar Boundary Layer
        l_bl = laminar_boundary_layer(vw,fp.eta,fs.ue,fs.rho,fs.mu,fs.nu,fp.x,fs.Pr,fp.L,fp.w,fs.T,fp.Tp,fp.ETA,fp.X)
        Cd_l.append(l_bl.Cd)
        xt.append(l_bl.xt)
   

        #Creation of a turbulent Boundary Layer
        model =["log_wake","seventh_law"]
        model = model[0]
        t_bl = turbulent_boundary_layer(fs.ue,fs.rho,fs.mu,fs.nu,fp.x,fp.L,fp.w,fp.Y,model)
        Cd_t.append(t_bl.Cd)

        # Creation of a transiotional boundary layer
        if l_bl.xt !=0:
            
            tl_bl = transition_boundary_layer(fs.ue,fs.rho,fs.mu,fs.nu,fp.x,l_bl,fp.L,fp.w,fp.Y)
            Cd_l_t.append(t_bl.Cd)

        else:
            Cd_l_t.append(l_bl.Cd)


    ## Plot
    ## =================================================================================
    #bl = [l_bl,t_bl]

    #plotting.f_blasius(fp,fs,l_bl,delta_est)
    #plotting.boundary_layer(fp,fs,tl_bl)
    #plotting.cf(fp,fs,tl_bl)
    #plotting.tw(fp,fs,t_bl)

    plotting.drag_reynolds(Re,Cd_l,Cd_t,Cd_l_t)
    #plotting.xt_reynolds(Re,xt)

    if plot:
    	plt.show()

    if plot:
        plt.show()
    
if __name__ == "__main__":
    
    main(plot=True)
