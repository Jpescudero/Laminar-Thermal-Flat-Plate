import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import odeint,RK45


## Problem Specification
## ===========================================================

ue     = 1
nu    = 1.7894*10e-5

## Non dimensional Suction Velocity (~-1 or -2) 
vw    = 0.5

# Length of the plate
L    = 1

# Definition of the flat plate X = 1 grid
n    = 100
x    = np.linspace(0,L,n+1)
x[0] = 10e-20

## Mesh Y Grid boundary layer height
y = np.linspace(0,0.2,200)

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

n     = 300
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
s_guesses = np.linspace(0.01,0.8)

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
s=[0.01,0.8]

## ODE solver : Shooting Method
## ===========================================================

## Compute error of the first guess phi0
finit[2] = s[0]
f        = odeint(fblasius,finit,eta)
phi0     = f[-1][1] - beta

# Number of max
nmax= 100
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

## X analysis
## ***********************************************************

# Velocity
u = ue * fd

# Reynolds x
Rex = ue*x/nu

# Boundary layer thichkness with x

delta = np.interp(0.994, fd, eta)*(2**0.5)*x/(Rex**0.5)

# Boundary layer displacement thichkness with x
# eta- f when eta = inf
delta_s = (((2)**0.5) * (f[-1]-f[-2])/(eta[-1]-eta[-2])*x)/(Rex**0.5)

# Boundary layer momentum thichkness with x
delta_ss = 2*(fdd[0]/(2**0.5))*x/(Rex**0.5)

# Friction coefficient with x
cf = 2*(fdd[0]/(2**0.5))/(Rex**0.5)
cf[0] = 0

## Grid Real Velocities
## ***********************************************************

# Grid Generation
X, Y = np.meshgrid(x, y)

# Eta for each point of the grid
ETA = (ue/(2*nu*X))**0.5 * Y

# U for each point of the grid
U = ue * np.interp(ETA, eta, fd)

# V for each point of the grid
V = ((nu*ue)/(2*X))**0.5 * (ETA*np.interp(ETA, eta, fd) - np.interp(ETA, eta, f))


## Plots
## ===========================================================
## f-eta plot
fig2 = plt.figure(2)
plt.plot(eta,f ,'r-',linewidth=2,label='f')
plt.plot(eta,fd ,'b-',linewidth=2,label='f´')
plt.plot(eta,fdd ,'g-',linewidth=2,label='f´´' )
plt.title('Blasius Flat plate functions')
plt.xlabel(r' $\eta$')
plt.xlim(0,10)
plt.ylim(0,5)
plt.legend()
plt.grid()


## Velocity-eta
fig3 = plt.figure(3)
plt.plot(eta, u/ue,'k-',linewidth=2,label=r' $U/u_e$')
plt.title('Blasius Flat plate velocity')
plt.xlabel(r' $\eta$')
plt.ylabel(r' $U/u_e$')
plt.xlim(0,4)
plt.ylim(0,1)
plt.legend(loc = 'upper left')
plt.grid()
#plt.show()


## Boundary Layer plot
fig4 = plt.figure(4)
plt.pcolor(X, Y, U)
plt.plot(x, delta,'b-',linewidth=2,label=r' $\delta$')
plt.plot(x, delta_s,'r-',linewidth=2,label=r' $\delta^*$')
plt.plot(x, delta_ss,'k-',linewidth=2,label=r' $\delta^**$')
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

plt.show()


## Write Data
## ===========================================================

## Non Dimensional Data
## ***********************************************************

# vw = str(vw)

# # f
# f_write = {'eta': eta, 'f': f}
# f_write = pd.DataFrame(data=f)
# f_write.to_csv("f_vw=" + str(vw) +".csv")

# # fd
# fd_write = {'eta': eta, 'fd': fd}
# fd_write = pd.DataFrame(data=fd)
# fd_write.to_csv("fd_vw=" + str(vw) + ".csv")

# # fdd
# fdd_write = {'eta': eta, 'fdd': fdd}
# fdd_write = pd.DataFrame(data=fdd)
# fdd_write.to_csv("fdd_vw=" + str(vw) + ".csv")


# ## Data along x
# ## ***********************************************************

# ## Suction = 0
# X = {'x': x, 'Rex': Rex, 'delta': delta, 'delta_s': delta_s, 'delta_ss': delta_ss,'cf': cf}
# X = pd.DataFrame(data=X)
# X.to_csv("x_BL_values"+ "_U="+ str(U) +".csv")
