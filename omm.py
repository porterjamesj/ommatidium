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

params = [y_k,y_basal,e_k,E,ne,D_e,r25_k,Y,ny]

nmol = 3

# The function that determines the rates of change, which we solve to integrate
def f(y, t, nmol, params):
    # get the parameters back out to a usable form
    y_k = params[0]
    y_basal = params[1]
    e_k = params[2]
    E = params[3]
    ne = params[4]
    D_e = params[5]
    r25_k = params[6]
    Y = params[7]
    ny = params[8] 
    
    #reshape flat array into shape (nmol,5)
    c = np.reshape(y,[nmol,5])
    #make empty array to hold rates of change
    xprime = np.empty(c.shape)

    yan = c[0,:]
    egfr = c[1,:]
    R_25 = c[2,:]
    #R_34 = c[3,:]

    # rates of change for yan
    xprime[0,:] = y_k * (y_basal + (1-hill(egfr,E,ne))*(1-hill(yan,Y,ny))*(progfxn(t)-y_basal) - yan)
    # rates of change for egfr
    xprime[1,:] = e_k * (e_prod*(sens+R_25 - egfr)) + D_e*oneDdiff(egfr)
    # rate of change for R_25
    xprime[2,:] = r25_k * (rough*hill(egfr,0.2,4.0)*(1-hill(yan,1.0,2.0)) - 1.0*R_25)
    # rate of change for R_34
    #xprime[6,:] = 0.0 #r34_prod + hill(p1+pp2/E_r,1.0,npr)*(1-hill((sens+rough)/RS_r,1.0,nrr))*(1-hill(yan/Y_r,1.0,nyr)) - r34_deg * R_34
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
plt.plot(timerange,resols[:,0,1])
#plt.plot(progs[0],progs[1], marker = 'o')
plt.plot(r8s[0],r8s[1], marker = 'o')
plt.plot(r25s[0],r25s[1], marker = 'o')
#plt.plot(r34s[0],r34s[1], marker = 'o')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")




'''
# emcee parameter tuning

# The function we are trying to minimize
# Takes yan data into account as well as biology
# (e.g. pnt should be high in differentiated cells, correct markers in place)
def invrmse(params,nmol,progs):
    lparams = [float(val) for val in params]
    args = tuple([nmol,lparams])
    sol = odeint(f,precluster.flatten(),timerange,args)
    resols = sol.reshape([len(timerange),nmol,5])
    error = 0
    # calculate error in yan data
    for datapoint in progs:
        delta = resols[datapoint[0],0,2] - datapoint[1]
        error += delta**2
    return error

# set up sampler parameters
ndim,nwalkers = 6, 100

#inital vector of the system
p0 = [np.random.normal(val,0.1*val) for val in params]

sampler = emcee.EnsembleSampler(nwalkers, ndim, invrmse, args=[nmol,progs], threads = 2, live_dangerously = False)

# trying scipy.optimize.minimize

progprod = interplolate
'''
