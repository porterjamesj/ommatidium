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
N = 80.0
nny = 4.0

y_prod = 0.15
y_deg = 0.5
P_y = 0.5
npy = 4.0

p1_deg = 1.0
p1_prod = 0.01
p2_deg = 1.0
p2_prod = 2.0
pp2_deg = 1.0

Y_p = 0.5
nyp = 4.0

E_p = 1.0
npp = 2.0

phos = 0.2

E = 0.2
nse = 2.0
e_deg = 0.01
D_e = 0.2
e_y = 0.2

r25_prod = 0.1
r34_prod = 0.1
r25_deg = 0.5
r34_deg = 0.5

P_r = 1.0
Y_r = 1.0
npr = 2.0
nyr = 2.0

RS_r = 0.5
nrr = 2.0

params = [N, nny, y_prod, y_deg, P_y, npy, p1_deg, p1_prod, p2_deg, p2_prod, pp2_deg, Y_p, nyp, E_p, npp, phos, E, nse, e_deg, e_y, r25_prod, r34_prod, r25_deg, r34_deg, P_r, Y_r, npr, nyr, RS_r, nrr, e_y]

nmol = 7

# The function that determines the rates of change, which we solve to integrate
def f(y, t, nmol, params):
    # get the parameters back out to a usable form
    N = params[0]
    nny = params[1]
    y_prod = params[2]
    y_deg = params[3]
    P_y = params[4]
    npy = params[5]
    p1_deg = params[6]
    p1_prod = params[7]
    p2_deg = params[8]
    p2_prod = params[9]
    pp2_deg = params[10]
    Y_p = params[11]
    nyp = params[12]
    E_p = params[13]
    npp = params[14]
    phos = params[15]
    E = params[16]
    nse = params[17]
    e_deg = params[18]
    D_e = params[19]
    r25_prod = params[20]
    r34_prod = params[21]
    r25_deg = params[22]
    r34_deg = params[23]
    P_r = params[24]
    Y_r = params[25]
    npr = params[26]
    nyr = params[27]
    RS_r = params[28]
    nrr = params[29]
    e_y = params[30]
    
    #reshape flat array into shape (nmol,5)
    c = np.reshape(y,[nmol,5])
    #make empty array to hold rates of change
    xprime = np.empty(c.shape)

    yan = c[0,:]
    p1 = c[1,:]
    p2= c[2,:]
    pp2 = c[3,:]
    egfr = c[4,:]
    R_25 = c[5,:]
    R_34 = c[6,:]

    #rates of change for yan
    xprime[0,:] = y_prod + (1-hill((p1+pp2)/P_y,1.0,npy))*hill(t/N,1.0,nny) - y_deg *yan - e_y * egfr * yan
    #rates of change for p1
    xprime[1,:] = p1_prod + (1-hill(yan/Y_p,1.0,nyp))*hill(p2/E_p,1.0,npp) - p1_deg * p1
    #rates of change for p2
    xprime[2,:] = hill(t/N,1.0,nny) - phos * egfr * p2 - p2_deg * p2
    # rates of change for pp2
    xprime[3,:] = phos * egfr * p2 - pp2_deg * pp2
    # rates of change for egfr
    xprime[4,:] =  E * sens + D_e * oneDdiff(egfr) - egfr * e_deg
    # rate of change for R_25
    xprime[5,:] = r25_prod + hill((p1+pp2)/P_r,1.0,npr)*rough*(1-hill(yan/Y_r,1.0,nyr)) - r25_deg * R_25
    # rate of change for R_34
    xprime[6,:] = r34_prod + hill((p1+pp2)/P_r,1.0,npr)*(1-hill((sens+rough)/RS_r,1.0,nrr))*(1-hill(yan/Y_r,1.0,nyr)) - r34_deg * R_34
    return xprime.flatten()

#set up precluster
precluster = np.zeros((nmol,5))
precluster[0,:] = 0.27 # set initial low levels of yan
# set up prepattern
sens = np.zeros(5)
sens[2] = 1.0
rough = np.zeros(5)
rough[1],rough[3] = 1.0,1.0

timerange = range(0,250) #range to solve on

#solve it
args = tuple([nmol,params])
sol = odeint(f,precluster.flatten(),timerange,args)
resols = sol.reshape([len(timerange),nmol,5])
print "done"

# plot the concentration of a factor in a single cell over time
plt.clf()
plt.plot(timerange,resols[:,0,2])
plt.savefig("fig.png")

# nic's data

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

r8s[:,1] = r8s[:,1] * 0.0046875
r25s[:,1] = r25s[:,1] * 0.0046875
r34s[:,1] = r34s[:,1] * 0.0046875

r8s[:,0] = r8s[:,0] * 0.5
r25s[:,0] = r25s[:,0] * 0.5
r34s[:,0] = r34s[:,0] * 0.5

# look at data

plt.plot(r8s.transpose()[0],r8s.transpose()[1], marker = 'o')
plt.plot(r25s.transpose()[0],r25s.transpose()[1], marker = 'o')
plt.plot(r34s.transpose()[0],r34s.transpose()[1], marker = 'o')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")

# emcee parameter tuning

# The function we are trying to minimize
# Takes yan data into account as well as biology
# (e.g. pnt should be high in differentiated cells, correct markers in place)
def invrmse(params,nmol,r8s,r25s,r34s):
    args = tuple([nmol,params])
    sol = odeint(f,precluster.flatten(),timerange,args)
    resols = sol.reshape([len(timerange),nmol,5])
    error = 0
    # calculate error in yan data
    for datapoint in r8s:
        delta = resols[datapoint[0],0,2] - datapoint[1]
        error += delta**2
    for datapoint in r25s:
        delta = resols[datapoint[0],0,1] - datapoint[1]
        delta += resols[datapoint[0],0,3] - datapoint[1]
        error += delta**2
    for datapoint in r34s:
        delta = resols[datapoint[0],0,0] - datapoint[1]
        delta += resols[datapoint[0],0,4] - datapoint[1]
        error += delta**2
    for param in params:
        if param < 0 or param > 500:
            error += 100
    for num in np.nditer(resols):
        if np.isnan(num):
            error += 100
    return -np.sqrt(error / len(r8s)*5)

# 30 dimensions in the parameter space, 100 walkers
ndim,nwalkers = 31, 10

#inital vector of the system
p0 = [ [np.random.normal(val,0.1*val) for val in params] for i in xrange(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, invrmse, args=[nmol,r8s,r25s,r34s], threads = 2, live_dangerously = True)

