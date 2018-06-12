import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint,RK45,quad
import csv

## Problem Specification
## ===========================================================

ue   = 10
P    = 101325
T    = 300 

# Flat plate Temperature
Tp   = 270

## Non dimensional Suction Velocity (~-1 or -2) 
vw   = 0

# Length and with of the plate
L    = 1
w    = 1

# Air thermal conductivity 
k    = 0.02

# Air constant specific heat
cp   = 0.3

# Viscosity Sutherland´s Law
mu   = (1.458*10e-6 * T**(3/2)) / (T + 110.4)

# Ideal Gas
rho  = P/(287*T)

# Kinematic Viscosity
nu   = mu/rho

# Reynolds 
Re   = ue*L/nu
print("Re =")
print(round(Re,0))
print()

# Prant number
Pr   = cp*mu/k

print("Pr =")
print(round(Pr,3))
print()

## Grid Definition
## ===========================================================

# Definition of the flat plate X = 1 grid
n    = 100
x    = np.linspace(0,L,n+1)
x[0] = 10e-20

## Mesh Y Grid boundary layer height
y = np.linspace(0,0.4,300)

## Lineal Vector of etas
n     = 200
eta_f =	30
eta   = np.linspace(0, eta_f, n+1)

# Grid Generation
X, Y = np.meshgrid(x, y)

# Eta for each point of the grid
ETA = (ue/(2*nu*X))**0.5 * Y


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
#{f1  = beta}
#{f2  = }

## Boundary value for eta=infty blasius
beta  = 1.0 

## Initial values in eta = 0
#{theta0  = 0}
#{theta1  = s}  initial guess 

thetainit = [1,
			0]

## Initial values in eta = inf
#{theta0  = alfa}
#{theta1  = }

## Boundary value for eta=infty temperature
alfa  = 0 


## Shooting Method: Error Fucntion (sec)
## ===========================================================
def error_function(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0


## Blasius falt plate boundary layer: Phisic Model 
## ===========================================================
def fblasius(f, eta):
    return [f[1],f[2], -f[0]*f[2]]

## Initial Guesses
## ===========================================================

## Plotting Error Function to know the zeros
#s_guesses = np.linspace(0.01,0.8)

# phi = []
# for s_guess in s_guesses:
# 	finit[2] = s_guess
# 	f = odeint(fblasius,finit,eta)
# 	phi.append(f[-1][1] - beta)

# fig1 = plt.figure(1)
# plt.plot(s_guesses,phi)
# plt.title('Error Phi-function for the Blasius equation shooting')
# plt.ylabel('phi')
# plt.xlabel('s')
# plt.grid(b=True, which='both')

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

## Thermal falt plate boundary layer: Phisic Model 
## ===========================================================

def ftemp(theta,x,Pr,f,eta):
	return [theta[1], -Pr*np.interp(x,eta,f)*theta[1]]

# Initial Guesses
# ===========================================================

# # Plotting Error Function to know the zeros
# s_guesses = np.linspace(-5,5)

# phi = []
# for s_guess in s_guesses:
# 	thetainit[1] = s_guess
# 	theta = odeint(ftemp,thetainit,eta,args=(Pr,f,eta,))
# 	phi.append(theta[-1][0] - alfa)

# fig1 = plt.figure(1)
# plt.plot(s_guesses,phi)
# plt.title('Error Phi-function for the Blasius equation shooting')
# plt.ylabel('phi')
# plt.xlabel('s')
# plt.grid(b=True, which='both')

## Guessed values after seeing zeros at Phi plot
s=[-0.1,0.1]

## ODE solver : Shooting Method
## ===========================================================

## Compute error of the first guess phi0
thetainit[1] = s[0]
theta        = odeint(ftemp,thetainit,eta,args=(Pr,f,eta,))
phi0         = theta[-1][0] - alfa

# Number of max
nmax= 100
eps = 1.0e-3


for n in range(nmax):

	# Updating next guess
    thetainit[1] = s[1]

    # Solving the differenctial equation
    theta = odeint(ftemp,thetainit,eta,args=(Pr,f,eta,))

    # Computing Phi
    phi1 = theta[-1,0] - alfa
    
    # Camputing Error
    ds = error_function(phi0,phi1,s[0],s[1])

    # Updating initial guesses
    s[0]  = s[1]
    s[1]  += ds
    phi0 = phi1

    # Convergence criteria
    if (abs(ds)<=eps):
        break

theta, thetad = theta[:,0], theta[:,1]


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

# Shear Stress
tp = fdd[0]*mu*ue*((ue)/(2*nu*x))**0.5
tp[0] = tp[1] - ((tp[1]-tp[2])/(x[1]-x[2]))*x[1]

# Friction coefficient with x
cf = 2*(fdd[0]/(2**0.5))/(Rex**0.5)
cf[0] = cf[1] - ((cf[1]-cf[2])/(x[1]-x[2]))*x[1]

# Nusselt number
Nu = - (Rex/2)**0.5 * thetad[0]

## Drag 
D  = 2 * w * np.trapz(tp, x)
Cd = (2*D)/(rho*ue**2*w*L)

print("Drag Cd =")
print(round(Cd,4))

## Grid Real Velocities
## ***********************************************************

# U for each point of the grid
U = ue * np.interp(ETA, eta, fd)

# V for each point of the grid
V = ((nu*ue)/(2*X))**0.5 * (ETA*np.interp(ETA, eta, fd) - np.interp(ETA, eta, f))

# Thermal non dimensional temperature
THETA = np.interp(ETA,eta,theta)


Temp     = ((Tp - T) * THETA) + T 

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
plt.plot(x, delta_ss,'g-',linewidth=2,label=r' $\delta^**$')
plt.title('Blasius Flat plate boundary layer')
plt.xlabel(r' $x$')
#plt.ylabel(r' $U/u_e$')
plt.xlim(0,1)
plt.ylim(0,0.1)
plt.legend()
plt.grid()
#plt.show()

## Cf-x plot
fig5 = plt.figure(5)
plt.plot(x, cf,'k-',linewidth=2,label=r' $C_f$')
plt.title('Blasius Flat plate Skin Friction')
plt.xlabel(r' $x$')
plt.ylabel(r' $C_f$')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend()
plt.grid()

fig6 = plt.figure(6)
plt.plot(x, cf,'k-',linewidth=2,label=r' $C_f$')
plt.title('Blasius Flat plate Skin Friction')
plt.xlabel(r' $x$')
plt.ylabel(r' $C_f$')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend()
plt.grid()

fig7 = plt.figure(7)
plt.plot(eta,theta,linewidth=2,label="Pr = " + str(round(Pr,3)))
plt.title('Thermal Boundary Layer')
plt.xlabel(r' $\eta$')
plt.ylabel(r' $\theta$')
plt.xlim(0,30)
plt.ylim(0,1)
plt.legend()
plt.grid()

#plt.show()

fig8 = plt.figure(8)
plt.pcolor(X, Y, Temp)
plt.title('Thermal Boundary Layer: Temperature Distribution (K)')
plt.colorbar()
plt.xlabel(r' $x$')
plt.ylabel(r' $y$')


## Write Data
## ===========================================================

## Non Dimensional Data
## ***********************************************************

# f
f_write = {'eta': eta, 'f': f, 'fd': fd, 'fdd': fdd, 'theta': theta, 'thetad': thetad}
f_write = pd.DataFrame(data=f)
f_write.to_csv("eta_f_fd_fdd_theta_thetad.csv")


## Data along x
## ***********************************************************

X = {'x': x, 'Rex': Rex, 'delta': delta, 'delta_s': delta_s, 'delta_ss': delta_ss,'cf': cf, 'tp': tp, 'Nu': Nu}
X = pd.DataFrame(data=X)
X.to_csv("x_Rex_delta_deltas_deltass_cf_tp" +".csv")

plt.show()