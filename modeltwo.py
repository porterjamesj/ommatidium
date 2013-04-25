"""
Simplest model I could think of.
"""
import numpy as np
from support import *

# how many molecules?
nmol = 4

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
e_y =0.01

P = 0.33
np = 4.0

r25_k = 1.0
N_25 = 0.1
n25 = 4.0

r34_k = 1.0
N_34 = 0.1
n34 = 4.0

params = [y_k,y_basal,Y,ny,e_k,E,ne,e_prod,D_e,r25_k,N_25,n25,r34_k,N_34,n34]

def dXdts(y, t, nmol, params, progfxn, sens, rough):
    # get the parameters back out to a usable form
    y_k,y_basal,Y,ny,e_k,E,ne,e_prod,D_e,r25_k,N_25,n25,r34_k,N_34,n34 = params
    #reshape flat array into shape (nmol,5)
    c = np.reshape(y,[nmol,5])
    #make empty array to hold rates of change
    xprime = np.empty(c.shape)

    yan = c[0,:]
    egfr = c[1,:]
    pnt = c[2,:
    #R_25 = c[3,:]
    #R_34 = c[4,:]

    # rates of change for yan
    xprime[0,:] = y_k * (y_basal + (1-hill(pnt,P,np))*(1-hill(yan,Y,ny))*(progfxn(t)-y_basal) - e_y *egfr*yan  - yan)
    # rates of change for egfr
    xprime[1,:] = e_k * (e_prod*(sens+R_25+R_34) - egfr) + D_e*oneDdiff(egfr)
    # rate of change for R_25
    #xprime[2,:] = r25_k * (rough*hill(egfr,N_25,n25)*(1-hill(yan,Y,ny)) - 1.0*R_25)
    # rate of change for R_34
    #xprime[3,:] = r34_k * (hill(egfr,N_34,n34)*(1-hill(yan,Y,ny))*(1-hill(rough+sens,0.1,2.0)) - 1.0*R_34)
    return xprime.flatten()
