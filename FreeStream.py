##==================
# Freestream Class
# Inputs: 
# ue    : Free Stream velocity (m/s)
# L     : Characteristic Lenght (m)
#
# Outputs: 
# 
# Atributes



class freestream:

    def __init__(self,ue,L):
        
        self.ue         = ue     # Free stream veocity
        self.P          = 101325 # Pressure
        self.T          = 300    # Temperature
        self.k          = 0.02   # Air thermal conductivity 
        self.cp         = 0.3    # Air constant specific heat

        self.L          = L      # Charasteristic Lenght

        # Viscosity Sutherland´s Law
        self.mu   = (1.458*10e-6 * self.T**(3/2)) / (self.T + 110.4)

        # Ideal Gas
        self.rho  = self.P/(287*self.T)

        # Kinematic Viscosity
        self.nu   = self.mu/self.rho

        # Prant number
        self.Pr   = self.cp*self.mu/self.k

        # Reynolds
        self.Re   = self.ue*self.L/self.nu


    def print_data(self,delta_est,case):

        print(  case   , ":  ue =", self.ue, " m/s",
            " Re=" , round(self.Re,0)," Pr  = ", round(self.Pr,3),
            "delta ~ " , round(delta_est,5),"(m)")