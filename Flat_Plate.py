import numpy as np

class flat_plate:

    def __init__(self,L,Tp,vw):
        
        self.L          = L     # Lenght of the plate (m)
        self.w          = 1     # Withd of the plate  (w)
        self.Tp         = Tp    # Surface Temperature (K)
        self.vw         = vw    # Suction velocity    (non-dimensonal)  
        self.x,self.y   = [],[] # x, y grid vectors
        self.X,self.Y   = [],[] # X, Y Grid            
        self.eta        = []    # eta
        self.ETA        = []    # ETA Grid
    
    def grid(self,n_x,n_y,delta_est,ue,nu):

        ## Grid Generation (By default grid is not generated)

        # Mesh X Grid
        self.x    = self.L*np.logspace(-5,0,n_x+1)
        #self.x    = np.linspace(0,self.L,1000) 
        #self.x[0] = 10e-6     

        # Mesh Y Grid boundary layer height
        self.y    = np.linspace(0,50*delta_est,n_y+1)  

        ## Lineal Vector of etas
        self.eta  = np.linspace(0, 30, n_x+1)

        # Meshgrid
        self.X, self.Y   = np.meshgrid(self.x, self.y)

        # Non dimensional distance Meshgrid              
        self.ETA    = (ue/(2*nu*self.X))**0.5 * self.Y         
