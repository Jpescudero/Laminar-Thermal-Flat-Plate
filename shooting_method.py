import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import odeint,RK45


## Problem Specification
## ===========================================================

U     = 10
nu    = 1.7894*10e-5

## Non dimensional Suction Velocity (~-1 or -2) maximum vw = 0.619 provokes separation 
vw    = 0 

# Definition of the flat plate X = 1 grid
L    = 1
n    = 11
x    = np.linspace(0,L,n)
x[0] = 10e-20

## Phisic Model 
## ===========================================================
def fblasius(f, eta):
    return [f[1],f[2], -f[0]*f[2]]

## Shooting Method: Error Fucntion (sec)
## ===========================================================
def error_function(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0

## Initial Vectors
## ===========================================================

n     = 200
eta_f =	50
eta   = np.linspace(0, eta_f, n+1)

## Boundary conditions
## ===========================================================

## Initial values in eta = 0
#{f0  = 0}
#{f1  = 0}
#{f2  = s}  initial guess 

finit =[-(2**0.5)*vw,
		0,
		0]


## Initial values in eta = inf
#{f0  = }
#{f1  = 1}
#{f2  = }

## Boundary value for eta=infty
beta  = 1.0 

## Initial Guesses
## ===========================================================

## Plotting Error Function to know the zeros
s_guesses = np.linspace(0.01,5.0)

phi = []
for s_guess in s_guesses:
	finit[2] = s_guess
	f = odeint(fblasius,finit,eta)
	phi.append(f[-1][1] - beta)

fig1 = plt.figure(1)
plt.plot(s_guesses,phi)
plt.title('Error Phi-function for the Blasius equation shooting')
plt.ylabel('phi')
plt.xlabel('s')
plt.grid(b=True, which='both')

## Guessed values after seeing zeros at Phi plot
s=[0.1,0.8]

## ODE solver : Shooting Method
## ===========================================================

## Compute error of the first guess phi0
finit[2] = s[0]
f        = odeint(fblasius,finit,eta)
phi0     = f[-1][1] - beta

# Number of max
nmax=10
eps = 1.0e-3


for n in range(nmax):

	# Updating next guess
    finit[2] = s[1]

    # Solving the differenctial equation
    f = odeint(fblasius,finit,eta)

    # Computing Phi
    phi1 = f[-1,1] - beta
    
    # Camputing Error
    ds = error_function(phi0,phi1,s[0],s[1])

    # Updating initial guesses
    s[0]  = s[1]
    s[1]  += ds
    phi0 = phi1

    # Convergence criteria
    if (abs(ds)<=eps):
        break

f, fd, fdd = f[:,0], f[:,1], f[:,2]


## Post analysis
## ===========================================================

# Velocity
u = U * fd

# Reynolds x
Rex = U*x/nu

# Boundary layer thichkness with x
delta = 5.2*x/(Rex**0.5)

# Boundary layer displacement thichkness with x
delta_s = 1.72*x/(Rex**0.5)

# Boundary layer momentum thichkness with x
delta_ss = 0.664*x/(Rex**0.5)

# Friction coefficient with x
cf = 0.664/(Rex**0.5)
cf[0] = 0

## Mesh Y Grid boundary layer height
y = []
eta_grid = []
u_grid   = []

Y = {}
U_y = {}

for i in range(len(x)):

	# Definition of a vector o y until boundary layer height (list of a list)
	y.append(np.linspace(0,delta[i],100))

	# Dictionary
	Y[str(x[i])] = np.linspace(0,delta[i],100)

	# Definition of equivalent eta for each point of the grid
	eta_grid.append((U/(2*nu*x[i]))**0.5 * y[-1])

	# Velocity Calculation
	u_grid.append(U * np.interp(eta_grid[-1], eta, fd))

	# Dictionary
	U_y[str(x[i])] = U * np.interp(eta_grid[-1], eta, fd)


## Plots
## ===========================================================
## f-eta plot
fig2 = plt.figure(2)
plt.plot(eta,f ,'r-',linewidth=2,label='f')
plt.plot(eta,fd ,'b-',linewidth=2,label='f´')
plt.plot(eta,fdd ,'g-',linewidth=2,label='f´´' )
plt.title('Blasius Flat plate functions')
plt.xlabel(r' $\eta$')
plt.xlim(0,4)
plt.ylim(0,2)
plt.legend()
plt.grid()


## Velocity-eta
fig3 = plt.figure(3)
plt.plot(eta, u/U,'k-',linewidth=2,label=r' $U/u_e$')
plt.title('Blasius Flat plate velocity')
plt.xlabel(r' $\eta$')
plt.ylabel(r' $U/u_e$')
plt.xlim(0,4)
plt.ylim(0,1)
plt.legend(loc = 'upper left')
plt.grid()
#plt.show()


## f-eta plot
fig4 = plt.figure(4)
plt.plot(x, delta,'b-',linewidth=2,label=r' $\delta$')
plt.plot(x, delta_s,'r-',linewidth=2,label=r' $\delta^*$')
plt.plot(x, delta_ss,'g-',linewidth=2,label=r' $\delta^**$')
plt.title('Blasius Flat plate boundary layer')
plt.xlabel(r' $x$')
#plt.ylabel(r' $U/u_e$')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend()
plt.grid()
#plt.show()

## Cf-x plot
fig5 = plt.figure(5)
plt.plot(x, cf,'k-',linewidth=1,label=r' $C_f$')
plt.title('Blasius Flat plate Skin Friction')
plt.xlabel(r' $x$')
plt.ylabel(r' $C_f$')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend()
plt.grid()
#plt.show()

## velocity profile plot
fig6 = plt.figure(6)
plt.plot((u_grid[3])/U,y[3]/delta[3],'k-',linewidth=2,label=r' $x = 0.33$')
plt.title('Blasius Flat plate velocity profile x=0.33m')
plt.xlabel(r' $u/U$')
plt.ylabel(r' $y/\delta$')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(loc = 'upper left')
plt.grid()
plt.show()