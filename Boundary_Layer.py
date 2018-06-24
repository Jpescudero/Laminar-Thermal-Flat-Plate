import numpy as np

from Shooting_Method import *

class boundary_layer:

	def __init__(self,*args):
		
		# X Values
		self.Rex,self.delta,self.delta_s,self.delta_ss,self.H,self.tw,self.cf = [],[],[],[],[],[],[]

	    # Drag
		self.Cd = 0

	    # Velocities
		self.u,self.v = [],[]

	def Drag(self,L,w,ue,rho,x):

		D  = 2 * w * np.trapz(self.tw, x)
		self.Cd = (2*D)/(rho*ue**2*w*L)
	    

## Turbulent Boundary Layer
## ===========================================================================

class turbulent_boundary_layer(boundary_layer):

	def __init__(self,*args):
		# call the parent class constructor with default values for potato
		# growth rate = 1: food need = 6: water need =3 
		super().__init__(*args)

		# Auto_solve
		self.solve(*args)

	def seventh_law(self,ue,rho,mu,nu,x):

        # Reynolds x
		self.Rex = ue*x/nu

        # Boundary layer thichkness with x
		self.delta = 0.38  *x/(self.Rex**(1/5))

        # Boundary layer displacement thichkness with x
        # eta- f when eta = inf
		self.delta_s = 0.048  *x/(self.Rex**(1/5))

        # Boundary layer momentum thichkness with x
		self.delta_ss = 0.037  *x/(self.Rex**(1/5))

        # Shape factor
		self.H = self.delta_s/self.delta_ss

        # Friction coefficient with x
		self.cf = 0.059 *(nu/(ue*self.delta))**(1/5)
		#self.cf = 0.0456 *(nu/(ue*self.delta))**(1/4)

        # Shear Stress
		self.tw = 0.5*rho*ue**2*self.cf

	def log_wake(self,ue,rho,mu,nu,x):


		B = 5.3
		A = 0.62
		k = 0.41

		def fturbulent (Red,Rex,k,B,A):
			return  1/((1/k)*np.log(Red)+B-A)**2


		def fu_s (Red,Rex,k,B,A):
			return  1/((1/k)*np.log(Red)+B-A)

		# Reynolds x
		self.Re_x = ue*x/nu

		## Initial condition approximation
		Re_delta0 = (x[0]*0.048/self.Re_x[0]**(1/5))*ue/nu

		## Differential equation for turbulent boundary layers
		Re_delta = odeint(fturbulent,Re_delta0,self.Re_x,args=(k,B,A,),tcrit=[x[0]])
		u_s   = ue * fu_s (Re_delta.flatten(),self.Re_x,k,B,A)

		# Boundary layer displacement thichkness with x
		self.delta_s  = nu*Re_delta.flatten()/ue
		self.delta_ss = self.delta_s

		# Boundary layer thichkness with x
		self.delta   = ue*self.delta_s/u_s	

		# Shape factor
		self.H = self.delta_s/self.delta_s

		# Shear Stress
		self.tw 	= u_s**2*rho

		# Cf
		self.cf = self.tw/ (0.5*rho*ue**2)


	def grid_values(self,ue,Y):

        # U for each point of the grid
		self.u = ue * (Y/self.delta)*(1/7)

		self.v = self.u*0

	def solve(self,fp,fs,model):

		ue,rho,mu,nu,x,L,w,Y = fs.ue,fs.rho,fs.mu,fs.nu,fp.x,fp.L,fp.w,fp.Y

		if model == "seventh_law":
			self.seventh_law(ue,rho,mu,nu,x)
		elif model =="log_wake":
			self.log_wake(ue,rho,mu,nu,x)

		self.Drag(L,w,ue,rho,x)
		self.grid_values(ue,Y)


## Laminar Boundary Layer
## ===========================================================================

class laminar_boundary_layer(boundary_layer):

	def __init__(self,*args):

		super(laminar_boundary_layer,self).__init__(*args)

		self.vw  = args[-2]   # Suction ve3locity (m/s)

        # Initial Conditions for blasius and thermal equations
		self.f0        = [-(2**0.5)*self.vw,0,0]
		self.finf      = [0,1,0]
		self.theta0    = [1,0]
		self.thetainf  = [0,0]

        # Blasius solution
		self.f,self.fd,self.fdd = [],[],[]

        # Thermal Solution
		self.theta,self.thetad  = [],[]
		self.Nu =[]

        # Transition Point xt
		self.xt = 0

		# Falkner-Skan velocity
		self.ue = []

		# Auto_solve
		self.solve(*args)

	def fblasius(self,f,eta,beta):
        # Blasius differential equation for solving BL velocities
		return [f[1],f[2], -f[0]*f[2]-beta*(1-f[1]**2)]

	def fthermal(self,theta,x,Pr,f,eta):
    # Thermal differential equation for solving Temperature
		return [theta[1], -Pr*np.interp(x,eta,f)*theta[1]]

	def x_values(self,eta,beta,ue,rho,mu,nu,x,Pr):

		m = -beta/(beta-2)

        # Reynolds x
		self.Rex = self.ue*x/nu

        # Boundary layer thichkness with x
		self.delta = np.interp(0.994, self.fd, eta)*(2**0.5)*x/(self.Rex**0.5)

        # Boundary layer displacement thichkness with x
        # eta- f when eta = inf
		self.delta_s = x *((2-beta)/self.Rex)**0.5 * (self.f[-1]-self.f[-2])/(eta[-1]-eta[-2])

        # Boundary layer momentum thichkness with x
		self.delta_ss = x *((2-beta)/self.Rex)**0.5 *((self.fdd[0] - beta*(self.f[-1]-self.f[-2])/(eta[-1]-eta[-2]))/(1+beta)) 

        # Shape factor
		self.H = self.delta_s/self.delta_ss

        # Shear Stress
		self.tw = self.fdd[0]*mu*self.ue*(((self.ue)*(m+1))/(2*nu*x))**0.5
		self.tw[0] = self.tw[1] - ((self.tw[1]-self.tw[2])/(x[1]-x[2]))*x[1]

        # Friction coefficient with x
		self.cf = 2*(self.fdd[0]/(2**0.5))/(((2-beta)*self.Rex**0.5))
		self.cf[0] = self.cf[1] - ((self.cf[1]-self.cf[2])/(x[1]-x[2]))*x[1]

        # Nusselt number
		self.Nu = - (self.Rex/(2-beta))**0.5 * self.thetad[0]


        # x transition: Micjel Criteria (Only valid for Re between 4x10e5 - 7x10e6)
		for i in range(len(x)):
			if (self.delta_ss[i] * self.ue[i])/nu > 1.535 * ((self.ue[i]*x[i])/nu)**0.444:
				self.xt = x[i]
				break

	def grid_values(self,ETA,eta,beta,ue,nu,T,Tp,X):

        # U for each point of the grid
		self.u = self.ue * np.interp(ETA, eta, self.fd)

        # V for each point of the grid
		self.v = ((nu*self.ue)/((2-beta)*X))**0.5 * ((1-beta)*ETA*np.interp(ETA, eta, self.fd) - np.interp(ETA, eta, self.f))

        # Thermal non dimensional temperature
		THETA = np.interp(ETA,eta,self.theta)

        # Real Temperarure (K)
		Temp     = ((Tp - T) * THETA) + T 

	def solve(self,fp,fs,vw,beta):

		# Free Stream values
		ue,rho,mu,nu,Re,Pr,T = fs.ue,fs.rho,fs.mu,fs.nu,fs.Re,fs.Pr,fs.T
		# Flat plate values
		L,w,Tp,x,X,Y,eta,ETA = fp.L,fp.w,fp.Tp,fp.x,fp.X,fp.Y,fp.eta,fp.ETA

		# Falkner-Skan velocity
		self.ue = ue*x**(-beta/(beta-2))

        ## Guessed values after seeing zeros at Phi plot
		sm   = shooting_method()
		sm.s = [beta+10e-2,2]
		sm.q = 1
		#sm.s_guesses  = np.linspace(beta+10e-25,5)
		#sm.plot_error_function(self.fblasius,eta,self.f0,self.finf,(beta,))

        ## Blasius Solution 
		y = sm.execute(self.fblasius,eta,self.f0,self.finf,(beta,))
		self.f, self.fd, self.fdd = y[:,0], y[:,1], y[:,2]

		sm.s = [-0.1,0.1]
		sm.q = 0

		## Thermal Solution 
		y = sm.execute(self.fthermal,eta,self.theta0,self.thetainf,(Pr,self.f,eta,))
		self.theta, self.thetad = y[:,0], y[:,1]

		self.x_values(eta,beta,ue,rho,mu,nu,x,Pr)

		self.Drag(L,w,ue,rho,x)

		self.grid_values(ETA,eta,beta,ue,nu,T,Tp,X)


class transition_boundary_layer(boundary_layer):

	def __init__(self,*args):

		super(transition_boundary_layer,self).__init__(*args)

		# Auto_solve
		self.solve(*args)


	def log_wake (self,ue,rho,mu,nu,x,l_bl):

		id_xt   = np.where(x == l_bl.xt)
		x_lam   = x[0:id_xt[0][0]]
		x_turb  = x[id_xt[0][0]:len(x)] - l_bl.xt
		x_turb[0] = 10e-5

		#id_x1   = np.argwhere(np.diff(np.sign(t_bl.delta_s-l_bl.delta_s[id_xt] * np.ones(len(x)))) !=0).reshape(-1) + 0 

		# Virtual Origin 
		#x0  = l_bl.xt - x[id_x1] 

		#id_x0 = np.argwhere(np.diff(np.sign(x-(x[-1]-x0)*np.ones(len(x)))) !=0).reshape(-1) + 0

		B = 5.3
		A = 0.62
		k = 0.41

		def fturbulent (Red,Rex,k,B,A):
			return  1/((1/k)*np.log(Red)+B-A)**2


		def fu_s (Red,Rex,k,B,A):
			return  1/((1/k)*np.log(Red)+B-A)
		

		self.Rex = ue*x/nu

		# Reynolds x turbulent 
		Re_x_turb = ue*x_turb/nu

		#delta_0 = fsolve(dfdelta0,0.01,args=(l_bl.delta[id_xt],ue,nu,k,B,A))
		delta_0 	= l_bl.delta[id_xt]
		delta_s_0 	= l_bl.delta_s[id_xt]
		
		## Initial cnodition approximation
		Re_delta0 = (x_turb[0]*0.048/Re_x_turb[0]**(1/5))*ue/nu

		## Differential equation for turbulent boundary layers
		Re_delta = odeint(fturbulent,Re_delta0,Re_x_turb,args=(k,B,A,),tcrit=[x[0]])

		# Boundary layer displacement thichkness with x
		l_bl.delta_s[id_xt[0][0]:] = nu*Re_delta.flatten()/ue + l_bl.delta_s[id_xt]
		self.delta_s = l_bl.delta_s
		Re_delta_new = l_bl.delta_s[id_xt[0][0]:]*ue/nu

		## Momentum thickness
		l_bl.delta_ss[id_xt[0][0]:] = nu*Re_delta.flatten()/ue + l_bl.delta_ss[id_xt]
		self.delta_ss = l_bl.delta_ss

		# Non dimensional friction velocity
		u_s   = ue * fu_s (Re_delta_new.flatten(),Re_x_turb,k,B,A)

		# Boundary layer thichkness with x
		l_bl.delta[id_xt[0][0]:]   = nu*Re_delta.flatten()/u_s + l_bl.delta[id_xt]
		self.delta = l_bl.delta

		# Shape factor
		self.H = self.delta_s/self.delta_s

		# # Shear Stress
		l_bl.tw[id_xt[0][0]:] = u_s**2*rho
		self.tw 	= l_bl.tw
		
		# Cf
		self.cf = self.tw/ (0.5*rho*ue**2)

	def solve(self,fp,fs,l_bl):

		ue,rho,mu,nu,x,L,w,Y = fs.ue,fs.rho,fs.mu,fs.nu,fp.x,fp.L,fp.w,fp.Y

		self.log_wake(ue,rho,mu,nu,x,l_bl)

		self.Drag(L,w,ue,rho,x)

	


