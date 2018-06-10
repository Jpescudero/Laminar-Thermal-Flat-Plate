import numpy as np
from pycse import bvp 
from scipy.interpolate import interp1d
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

n = 100

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

# Definition of the flat plate and the grid
L = 1
x = np.linspace(10e-25,L,100)

u = U * fd

Rex = U*x/nu

delta = 5.2*x/(Rex**0.5)

delta_s = 1.72*x/(Rex**0.5)

theta = 0.664*x/(Rex**0.5)

cf = 0.664/(Rex**0.5)
cf[0] = 0

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
fig1 = plt.figure(2)
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
fig1 = plt.figure(3)
plt.plot(x, delta,'b-',linewidth=2,label=r' $\delta$')
plt.plot(x, delta_s,'r-',linewidth=2,label=r' $\delta^*$')
plt.plot(x, theta,'g-',linewidth=2,label=r' $\theta$')
plt.title('Blasius Flat plate boundary layer')
plt.xlabel(r' $x$')
#plt.ylabel(r' $U/u_e$')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend()
plt.grid()
#plt.show()

## f-eta plot
fig1 = plt.figure(4)
plt.plot(x, cf,'k-',linewidth=2,label=r' $C_f$')
plt.title('Blasius Flat plate Skin Friction')
plt.xlabel(r' $x$')
plt.ylabel(r' $C_f$')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend()
plt.grid()
plt.show()