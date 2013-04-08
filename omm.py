'''
This is intended to be a representation of a single ommatidium with steadily rising levels of Notch signaling and prepatterning of Rough, Lozenge, Sens, Ato
'''

#imports
from support import *
import numpy as np
import emcee
from data import *
from scipy import interpolate
from scipy.optimize import minimize as spopt
#1D diffusions function

def oneDdiff(concs):
    '''takes a 1D array of concentrations and returns the diffusion rate for each'''
    diff_rates = np.zeros(concs.shape)
    for i in range(concs.size):
        # try-except statements catch the boundaries of the array
        try:
            lside = (concs[i-1] - concs[i])
        except IndexError:
            lside = 0.0
        try:
            rside = (concs[i+1] - concs[i])
        except IndexError:
            rside = 0.0
        diff_rates[i] = lside + rside
    return diff_rates

# basically using this as a notch input over time
progfxn = interpolate.interp1d(progs[0],progs[1])

#parameter definitions

y_k = 0.5
y_basal = 0.2
Y = 1.0
ny = 2.0

e_k = 0.01
E = 0.33
ne = 4.0
e_prod = 1.3
D_e = 0.01

r25_k = 1.0
N_25 = 0.1
n25 = 4.0

r34_k = 1.0
N_34 = 0.1
n34 = 4.0

params = [y_k,y_basal,Y,ny,e_k,E,ne,e_prod,D_e,r25_k,N_25,n25,r34_k,N_34,n34]

nmol = 4

# The function that determines the rates of change, which we solve to integrate
def f(y, t, nmol, params):
    # get the parameters back out to a usable form
    y_k = params[0]
    y_basal = params[1]
    Y = params[2]
    ny = params[3]
    e_k = params[4]
    E = params[5]
    ne = params[6]
    e_prod = params[7]
    D_e = params[8]
    r25_k = params[9]
    N_25 = params[10]
    n25 = params[11]
    r34_k = params[12]
    N_34 = params[13]
    n34 = params[14] 
    
    #reshape flat array into shape (nmol,5)
    c = np.reshape(y,[nmol,5])
    #make empty array to hold rates of change
    xprime = np.empty(c.shape)

    yan = c[0,:]
    egfr = c[1,:]
    R_25 = c[2,:]
    R_34 = c[3,:]

    # rates of change for yan
    # add repression by r25 and r34 here
    xprime[0,:] = y_k * (y_basal + (1-hill(egfr,E,ne))*(1-hill(yan,Y,ny))*(progfxn(t)-y_basal) - yan)
    # rates of change for egfr
    xprime[1,:] = e_k * (e_prod*(sens+R_25+R_34) - egfr) + D_e*oneDdiff(egfr)
    # rate of change for R_25
    xprime[2,:] = r25_k * (rough*hill(egfr,N_25,n25)*(1-hill(yan,Y,ny)) - 1.0*R_25)
    # rate of change for R_34
    xprime[3,:] = r34_k * (hill(egfr,N_34,n34)*(1-hill(yan,Y,ny))*(1-hill(rough+sens,0.1,2.0)) - 1.0*R_34)
    return xprime.flatten()

#set up precluster
precluster = np.zeros((nmol,5))
precluster[0,:] = 0.27 # set initial low levels of yan
# set up prepattern
sens = np.zeros(5)
sens[2] = 1.0
rough = np.zeros(5)
rough[1],rough[3] = 1.0,1.0

timerange = range(0,249) #range to solve on

#solve it
args = tuple([nmol,params])
sol = odeint(f,precluster.flatten(),timerange,args)
resols = sol.reshape([len(timerange),nmol,5])
print "done"

# plot simulation and data
plt.clf()
plt.plot(timerange,resols[:,0,0])
plt.plot(timerange,resols[:,0,1])
plt.plot(timerange,resols[:,0,2])
plt.plot(progs[0],progs[1], marker = 'o')
plt.plot(r8s[0],r8s[1], marker = 'o')
plt.plot(r25s[0],r25s[1], marker = 'o')
plt.plot(r34s[0],r34s[1], marker = 'o')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")


# parameter tuning

t_r8s = r8s.transpose()
t_r25s = r25s.transpose()
t_r34s = r34s.transpose()

res = None

def err_fxn(params,nmol):
    lparams = [float(val) for val in params]
    args = tuple([nmol,lparams])
    sol = odeint(f,precluster.flatten(),timerange,args)
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

def callback(ps):
    res = ps

# trying scipy.optimize.minimize

res = spopt(err_fxn,params,
            args=[nmol],
            method = "L-BFGS-B",
            bounds = [(0,None) for i in xrange(0,15)],
            tol = 0.075,
            options = {"maxiter":3,"disp":True},
            callback = callback)
