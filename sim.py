'''
This is intended to be a representation of a single ommatidium with steadily rising levels of Notch signaling and prepatterning of Rough, Lozenge, Sens, Ato
'''

#imports
import numpy as np
import emcee
from data import *
from scipy.integrate import odeint
from scipy import interpolate
from scipy.optimize import minimize as spopt

# what model are we using?
from simplemodel import *

# basically using this as a notch input over time
progfxn = interpolate.interp1d(progs[0],progs[1])

#set up precluster
precluster = np.zeros((nmol,5))
precluster[0,:] = 0.27 # set initial low levels of yan

# set up prepattern
sens = np.zeros(5)
sens[2] = 1.0
rough = np.zeros(5)
rough[1],rough[3] = 1.0,1.0

timerange = range(0,249) #range to solve on

# solve it
# can change the function to use a different one from funcs.py
args = (nmol,optparams,progfxn,sens,rough)
sol = odeint(dXdts,precluster.flatten(),timerange,args)
resols = sol.reshape([len(timerange),nmol,5])
print "done"

# plot simulation and data
plt.clf()
plt.plot(timerange,resols[:,0,2],color = 'r')
plt.plot(timerange,resols[:,0,1],color = 'g')
plt.plot(timerange,resols[:,0,0],color = 'b')
plt.plot(progs[0],progs[1],'ko-')
plt.plot(r8s[0],r8s[1],'ro')
plt.plot(r25s[0],r25s[1],'go')
plt.plot(r34s[0],r34s[1],'bo')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")

# parameter tuning

t_r8s = r8s.transpose()
t_r25s = r25s.transpose()
t_r34s = r34s.transpose()

def err_fxn(params,nmol):
    lparams = [float(val) for val in params] # this might not be necessary
    args = (nmol,lparams,progfxn,sens,rough)
    sol = odeint(dXdts,precluster.flatten(),timerange,args)
    resols = sol.reshape([len(timerange),nmol,5])
    error = 0
    # calculate error in yan data
    for datapoint in t_r8s:
        delta = resols[datapoint[0],0,2] - datapoint[1]
        error += delta**2
    for datapoint in t_r25s:
        delta = resols[datapoint[0],0,1] - datapoint[1]
        error += delta**2
    for datapoint in t_r34s:
        delta = resols[datapoint[0],0,0] - datapoint[1]
        error += delta**2
    return error

# trying scipy.optimize.minimize

res = spopt(err_fxn,params,
            args=[nmol],
            method = "L-BFGS-B",
            bounds = [(0,None) for i in xrange(0,15)],
            tol = 0.075,
            options = {"maxiter":3,"disp":True})
