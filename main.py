import numpy as np
from pycse import bvp 
import pandas as pd
import matplotlib.pyplot as plt

## Constants
## ===========================================================

U     = 1
nu    = 1.846*10e-5

## Boundary Conditions
## ===========================================================

def bcfun(Y):
    fa, fb = Y[0, :], Y[-1, :]
    return [fa[0],        # f1(0) =  0
            fa[1],        # f2(0) = 0
            1.0 - fb[1]]  # f2(inf) = 1

## Phisic Model 
## ===========================================================

## Blasius Equation (f´´´ + f*f´´ = 0)
def Blasius(F,eta):
	f1, f2, f3 = F.T
	return np.column_stack([f2,f3,-f1 * f3])

## ODE solver 
## ===========================================================

# Initial Vectors
eta = np.linspace(0, 4, 100)
f1init = eta
f2init = np.exp(-eta)
f3init = np.exp(-eta)

Finit = np.column_stack([f1init, f2init, f3init])

## Solve ODE
sol = bvp(Blasius, bcfun, eta, Finit)
f, fd, fdd = sol.T


## Post analysis
## ===========================================================

# Definition of the flat plate X grid
L    = 1
n    = 11
x    = np.linspace(0,L,n)
x[0] = 10e-20

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

## Export data
## ===========================================================

# X data
X = {'x': x, 'Rex': Rex, 'delta': delta, 'delta_s': delta_s, 'delta_ss': delta_ss,'cf': cf}
X = pd.DataFrame(data=X)
X.set_index('x',inplace=True)
X.to_csv("x_BL_values.csv")

# Y data for each station
Y = pd.DataFrame(data=Y)
Y.to_csv("y_xstations.csv")

# U(Y) data for each station
U_y = pd.DataFrame(data=U_y)
U_y.to_csv("u(y)_xstations.csv")

## Plots
## ===========================================================

## f-eta plot
fig1 = plt.figure(1)
plt.plot(eta,f ,'r-',linewidth=2,label='f')
plt.plot(eta,fd ,'b-',linewidth=2,label='f´')
plt.plot(eta,fdd ,'g-',linewidth=2,label='f´´' )
plt.title('Blasius Flat plate functions')
plt.xlabel(r' $\eta$')
plt.xlim(0,4)
plt.ylim(0,2)
plt.legend()
plt.grid()
#plt.show()


## Velocity-eta
fig2 = plt.figure(2)
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
fig3 = plt.figure(3)
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
fig4 = plt.figure(4)
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
fig5 = plt.figure(5)
plt.plot((u_grid[3]),y[3]/delta[3],'k-',linewidth=2,label=r' $x = 0.33$')
plt.title('Blasius Flat plate velocity profile x=0.33m')
plt.xlabel(r' $u/U$')
plt.ylabel(r' $y/\delta$')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(loc = 'upper left')
plt.grid()
plt.show()