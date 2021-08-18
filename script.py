import numpy as np
import matplotlib.pyplot as plt

# Set the parameters
## Number of grid points
npoints = 101

## Minimum and Maximum extent of the domain
xmin = 0.0
xmax = 1.0

## Maximum number of iterations to run
maxiter = 100

## wave speed
a = 1.0

## CFL number
cfl = 0.25

# Define grid
x, dx= np.linspace(xmin, xmax, npoints, endpoint=True, retstep=True)


# Declare the four different solution array
## usoln->Exact solution for comparision
## uftbs->FTBS scheme
## uftcs->FTCS scheme
## uftfs->FTFS scheme
usoln = np.ones(npoints)
uftbs = np.ones(npoints)
uftcs = np.ones(npoints)
uftfs = np.ones(npoints)

def squrewave(x1,x2,x, value):
  '''
  function return a square wave for given inputs
  Inputs:
  x -> grid on which square wave is defined
  x1 -> minimum extent of non-zero value
  x2 -> maximum extent of non-zero value
  value->
  '''
  u = np.zeros(x.size)
  i = np.abs(x-x1).argmin()
  j = np.abs(x-x2).argmin()
  u[i:j]=value
  return u

#inialize to square wave
usoln = squrewave(0.2,0.4,x,2)
uftbs = squrewave(0.2,0.4,x,2)
uftcs = squrewave(0.2,0.4,x,2)
uftfs = squrewave(0.2,0.4,x,2)
ulxwd = squrewave(0.2,0.4,x,2)

#comment out the below section if
# you just want to analyze square wave
#sine wave
k = 8 #wave number
usoln = np.sin(2*np.pi*k*x)
uftbs = np.sin(2*np.pi*k*x)
uftcs = np.sin(2*np.pi*k*x)
uftfs = np.sin(2*np.pi*k*x)
ulxwd = np.sin(2*np.pi*k*x)

# interactive plot
#plt.ion()
fig, axs = plt.subplots(5)
fig.suptitle("Linear Wave equation")
lnsoln, = axs[0].plot(x,usoln, 'b-', label='True Solution')
axs[0].set_ylabel('exact solution')
lnftbs, = axs[1].plot(x,uftbs, 'b-', label='FTBS')
axs[1].set_ylabel('FTBS')
lnftcs, = axs[3].plot(x,uftcs, 'b-', label='FTCS')
axs[3].set_ylabel('FTCS')
lnftfs, = axs[4].plot(x,uftcs, 'b-', label='FTFS')
axs[4].set_ylabel('FTFS')
lnlxwd, = axs[2].plot(x,ulxwd, 'b-', label='Lax-Wendroff')
axs[2].set_ylabel('Lax-Wendroff')
axs[4].set_xlabel('Domain')


#linear wave equation
dt = cfl * dx / np.abs(a)
for iter in range(maxiter):
    usoln[:] =  np.sin(2*np.pi*k*(x-a*(iter*dt)))
    #usoln = squrewave(0.2+a*iter*dt, 0.4+a*iter*dt,x,2)
    uftbs[1:] = uftbs[1:] - a*dt/dx*(uftbs[1:] - uftbs[:-1])
    uftbs[0] = uftbs[-1]
    uftcs[1:-1] = uftcs[1:-1] - 0.5*a*dt/dx*(uftcs[2:] - uftcs[:-2])
    uftcs[0] = uftcs[0] - 0.5*a*dt/dx*(uftcs[1]-uftcs[-2])
    uftcs[-1] = uftcs[0]
    uftfs[:-1] = uftfs[:-1] - a*dt/dx*(uftfs[1:] - uftfs[:-1])
    uftfs[-1] = uftfs[0]
    ulxwd[1:-1] = 0.5*(cfl*(1+cfl)*ulxwd[:-2] + 2*(1-cfl**2)*ulxwd[1:-1] - cfl*(1-cfl)*ulxwd[2:])
    ulxwd[0] = 0.5*(cfl*(1+cfl)*ulxwd[-2] + 2*(1-cfl**2)*ulxwd[0] - cfl*(1-cfl)*ulxwd[1])
    ulxwd[-1] = ulxwd[0]
    lnsoln.set_ydata(usoln)
    lnftbs.set_ydata(uftbs)
    lnftcs.set_ydata(uftcs)
    lnftfs.set_ydata(uftfs)
    lnlxwd.set_ydata(ulxwd)
    # reset the yrange while plotting
    #for i in range(4):
    #   axs[i].relim()
    #   axs[i].autoscale_view()
    plt.pause(0.01)
    #fig.savefig('image'+str(iter).zfill(3)+'.png')


plt.show()
