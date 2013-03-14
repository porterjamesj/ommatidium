'''
This is intended to be a representation of a single ommatidium with steadily rising levels of Notch signaling and prepatterning of Rough, Lozenge, Sens, Ato
'''

#imports
from support import *
import numpy as np
import scipy as sp
import emcee

#1D diffusions function

def oneDdiff(concs):
    #takes a 1D array of concentrations and returns the diffusion rate for each
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

#parameter definitions

y_prod = 0.27
E_y = 0.5
npy = 2.0
N = 62.0
y_deg = 1.0
k_y = 1.0

E = 2.0

params = [y_prod,E_y,npy,N,y_deg,k_y]

nmol = 1

# The function that determines the rates of change, which we solve to integrate
def f(y, t, nmol, params):
    # get the parameters back out to a usable form
    y_prod = params[0]
    E_y = params[1]
    npy = params[2]
    N = params[3]
    y_deg = params[4]
    k_y = params[5]
    
    #reshape flat array into shape (nmol,5)
    c = np.reshape(y,[nmol,5])
    #make empty array to hold rates of change
    xprime = np.empty(c.shape)

    yan = c[0,:]
    #egfr = c[1,:]
    #R_25 = c[2,:]
    #R_34 = c[3,:]

    # time dependance of Yan induction
    if t < 120:
        ind = float(t)
    else:
        ind = 120.0 * (120.0 / float(t))

    #(1-hill((egfr)/E_y,k_y,npy))*
    #rates of change for yan
    xprime[0,:] = y_prod + hill(ind/N,k_y,nny) - y_deg *yan 
    # rates of change for egfr
    #xprime[4,:] = 0.0 #E * (sens+R_25+R_34) + D_e * oneDdiff(egfr) - egfr * e_deg
    # rate of change for R_25
    #xprime[5,:] = 0.0 #r25_prod + hill(egfr/E_r,1.0,npr)*rough*(1-hill(yan/Y_r,1.0,nyr)) - r25_deg * R_25
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

timerange = range(0,300) #range to solve on

#solve it
args = tuple([nmol,params])
sol = odeint(f,precluster.flatten(),timerange,args)
resols = sol.reshape([len(timerange),nmol,5])
print "done"

# plot simulation and data
plt.clf()
plt.plot(timerange,resols[:,0,2])
plt.plot(progs.transpose()[0],progs.transpose()[1], marker = 'o')
#plt.plot(r8s.transpose()[0],r8s.transpose()[1], marker = 'o')
#plt.plot(r25s.transpose()[0],r25s.transpose()[1], marker = 'o')
#plt.plot(r34s.transpose()[0],r34s.transpose()[1], marker = 'o')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")


# nic's data

progs = np.array([[0,58],
                 [31,61],
                 [62,68],
                 [93,74],
                 [111,91],
                 [140,139],
                 [169,200],
                 [195,242],
                 [227,257],
                 [254,256],
                 [283,251],
                 [312,242],
                 [340,221],
                 [369,190],
                 [398,171],
                 [425,151],
                 [454,136],
                 [483,128],
                 [512,122],
                 [541,113]],dtype = float)

r8s = np.array([[0,58],
                [31,61],
                [62,68],
                [93,74],
                [192,74],
                [254,71],
                [296,60],
                [324,65],
                [358,60],
                [414,49],
                [475,43]],dtype = float)

r25s = np.array([[0,58],
                 [31,61],
                 [62,68],
                 [93,74],
                 [111,91],
                 [221,91],
                 [271,74],
                 [323,65],
                 [390,63],
                 [434,60],
                 [476,57]],dtype = float)

r34s = np.array([[0,58],
                 [31,61],
                 [62,68],
                 [93,74],
                 [111,91],
                 [138,133],
                 [229,147],
                 [257,133],
                 [283,107],
                 [312,90],
                 [340,80],
                 [368,75],
                 [397,71],
                 [426,66],
                 [491,58]], dtype = float)


# scale down to something closer to the simulation

progs[:,1] = progs[:,1] * 0.0046875
r8s[:,1] = r8s[:,1] * 0.0046875
r25s[:,1] = r25s[:,1] * 0.0046875
r34s[:,1] = r34s[:,1] * 0.0046875

progs[:,0] = progs[:,0] * 0.5
r8s[:,0] = r8s[:,0] * 0.5
r25s[:,0] = r25s[:,0] * 0.5
r34s[:,0] = r34s[:,0] * 0.5

# look at data

plt.plot(progs.transpose()[0],progs.transpose()[1], marker = 'o')
plt.plot(r8s.transpose()[0],r8s.transpose()[1], marker = 'o')
plt.plot(r25s.transpose()[0],r25s.transpose()[1], marker = 'o')
plt.plot(r34s.transpose()[0],r34s.transpose()[1], marker = 'o')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")

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
    return -((10*error)**2)

# set up sampler parameters
ndim,nwalkers = 6, 100

#inital vector of the system
p0 = [ [np.random.normal(val,0.1*val) for val in params] for i in xrange(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, invrmse, args=[nmol,progs], threads = 2, live_dangerously = False)
