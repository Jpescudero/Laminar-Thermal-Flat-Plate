import numpy as np
from scipy.integrate import odeint

class shooting_method():

    def __init__(self):

        self.nmax       = 200       # Maximum number of iterations
        self.eps        = 10e-4     # Error Allowed
        self.s          = [-1,1]    # Guessed values
        self.p          = -1        # position of the guess in vector df/dx
        self.q          = 0         # position of the boundary value in vector f(inf)
        self.s_guesses  = np.linspace(0.01,0.8)

    def plot_error_function(self,df,x,f0,finf,args_list):
        phi = []
        for s_guess in self.s_guesses:
            f0[-1] = s_guess
            if len(args_list)== 0: 
                f = odeint(df,f0,x)
            else:
                f = odeint(df,f0,x,args = args_list)

            phi.append(f[-1][self.q] - finf[self.q])

        fig1 = plt.figure(1)
        plt.plot(self.s_guesses,phi)
        plt.title('Error Phi-function for the Blasius equation shooting')
        plt.ylabel('phi')
        plt.xlabel('s')
        plt.grid(b=True, which='both')

    def error_function(self,phi0,phi1,s0,s1):
        if (abs(phi1-phi0)>0.0):   
            return    -phi1 *(s1 - s0)/float(phi1 - phi0)
        else:
            return 0.0

    def first_guess(self,df,x,f0,finf,args_list):

        f0[self.p]  = self.s[0]

        if len(args_list)== 0: 
            f = odeint(df,f0,x)
        else:
            f = odeint(df,f0,x,args = args_list)
        
        phi0 = f[-1][self.q] - finf[self.q]

        return phi0

    def execute(self,df,x,f0,finf,args_list):

        phi0 = self.first_guess(df,x,f0,finf,args_list)

        for n in range(self.nmax):

            # Updating next guess
            f0[self.p] = self.s[1]

            # Solving the differenctial equation
            if len (args_list)==0: 
                f = odeint(df,f0,x)
            else:
                f = odeint(df,f0,x,args = args_list)

            # Computing Phi
            phi1 = f[-1][self.q] - finf[self.q]

            # Computing Error
            ds = self.error_function(phi0,phi1,self.s[0],self.s[1])

            # Updating initial guesses
            self.s[0]  = self.s[1]
            self.s[1]  += ds
            phi0 = phi1

            # Convergence criteria
            if (abs(ds)<=self.eps):
                break

        return f