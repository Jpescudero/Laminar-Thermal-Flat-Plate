import numpy as np
import matplotlib.pyplot as plt


from Flat_Plate import *
from FreeStream import *
from Boundary_Layer import *
from Plotting import *



def main(plot=False):

	## Problem Specification
    ## ==========================================================================
    L  = 1      # Lenght of the plate (m)  
    Tp  = 300   # Plate Surface Temperature (K)  

	# Creation of a flat Plate Object: (Length, Wall Temperature,suction velocity)
    fp = flat_plate(L,Tp)

    ## Laminar Case
    ## ===========================================================================
    # Freestream velocity (m/s)
    ue = 10    

    # Creation of a Free stream object
    fs = freestream(ue,L)

    # Estimation of boundary layer height
    delta_est = L/(fs.Re)**0.5

    # Print information
    case = "Laminar Case"
    fs.print_data(delta_est,case)

    # Creation of the flat plate grid
    fp.grid(delta_est,ue,fs.nu)
    
    ## Free Stream
    ## ************************************************************************* 

    # Suction Non Dimensional Velocity
    vw = [-2,-1,0,0.2,0.45]
    
    # Falkner-Skan wedge angles (rad)
    beta = [-0.15,-0.1,0,0.1,0.3,0.5]

    # Prand Number
    Pr   = [0.1,0.03,0.1,0.3,1,10]

    # Creation of a dictionary
    parameters = {}
    parameters["vw"] = vw
    parameters["beta"] = beta
    parameters["Pr"] = Pr

    ## Suction
    ## *************************************************************************
    # Creation of a dictionary of boundary layers
    laminars_bl = {}

    for idx,VW in enumerate(vw):

        # Creation of a laminar boundary layer for each vw=0
        l_bl = laminar_boundary_layer(fp,fs,vw[idx],0)

        # Updating dictionary
        laminars_bl[str(vw[idx])]=l_bl
    
    # Plot
    plotting.laminar(fp,fs,laminars_bl,delta_est,parameters["vw"],"vw",1)

    ## Pressure Gradient
    ## *************************************************************************
    # Creation of a dictionary of boundary layers
    laminars_bl = {}

    for idx,BETA in enumerate(beta):

        # Creation of a laminar boundary layer for each vw=0
        l_bl = laminar_boundary_layer(fp,fs,0,beta[idx])

        # Updating dictionary
        laminars_bl[str(beta[idx])]=l_bl
    
    # Plot
    plotting.laminar(fp,fs,laminars_bl,delta_est,parameters["beta"],"beta",2)

    ## Temperature
    ## ***************************************************************************
    # Creation of a dictionary of boundary layers
    laminars_bl = {}

    for idx,PR in enumerate(Pr):

        # Creation of a laminar boundary layer for each vw=0
        fs.Pr = Pr[idx]
        l_bl  = laminar_boundary_layer(fp,fs,0,0)

        # Updating dictionary
        laminars_bl[str(Pr[idx])]=l_bl
    
    # Plot
    plotting.laminar_thermal(fp,fs,laminars_bl,delta_est,parameters["Pr"],"Pr",3)
    

    ## Turbulent Case
    ## ===========================================================================
    

    # Reynolds Influence
    # ============================================================================

    # Plot Vectors
    Re      = []
    Cd_l    = []
    Cd_t    = []
    Cd_l_t  = []
    xt      = []
    
    for ue in np.logspace(2,4,100):

    	# Creation of a Free stream: (Free Stream Velocity, Characteristic Length)
        fs = freestream(ue,L)
    	
        Re.append(fs.Re)

        # Estimation of boundary layer height
        delta_est = L/(fs.Re)**0.5

        # Creation of the flat plate grid
        fp.grid(delta_est,ue,fs.nu)

    	#Creation of a laminar Boundary Layer
        l_bl = laminar_boundary_layer(fp,fs,0,0)
        Cd_l.append(l_bl.Cd)
        xt.append(l_bl.xt)
   

        #Creation of a turbulent Boundary Layer
        model =["log_wake","seventh_law"]
        model = model[0]
        t_bl = turbulent_boundary_layer(fp,fs,model)
        Cd_t.append(t_bl.Cd)

        # Creation of a transiotional boundary layer
        if l_bl.xt !=0:
            tl_bl = transition_boundary_layer(fp,fs,l_bl)
            Cd_l_t.append(t_bl.Cd)
        else:
            Cd_l_t.append(l_bl.Cd)


    plotting.drag_reynolds(Re,Cd_l,Cd_t,Cd_l_t,4)
    plotting.xt_reynolds(Re,xt,5)

    if plot:
    	plt.show()

    if plot:
        plt.show()
    
if __name__ == "__main__":
    
    main(plot=True)
