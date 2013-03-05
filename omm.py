'''
This is intended to be a representation of a single ommatidium with steadily rising levels of Notch signaling and prepatterning of Rough, Lozenge, Sens, Ato
'''

#//

#imports
from support import *
import numpy as np
import scipy as sp

#//


#the signaling function
#signals can be recieved from anywhere else in the tissue
def sig(conc_array, triang):
    conc_list = conc_array.flatten(order ='A')
    val_array = np.zeros(conc_list.size)
    for pindex in range(len(triang.points)):
        for oindex in range(len(triang.points)):
            if oindex != pindex:
                val_array[pindex] += (1 / np.linalg.norm(triang.points[pindex] - triang.points[oindex])**2) * (conc_list[oindex])
    return np.reshape(val_array, np.shape(conc_array))

#//
#parameter definitions

N_y = 60.0
nny = 4.0
y_deg = 1.0

ney = 2.0

nmol = 1
xdim = 5
ydim = 5

#//
# note here that nmol does not include the molecules
# in the lubensky model (ato, sense, u, and h)
def f(y, t, nmol, triang):
    #reshape flat array into shape (nmol,x,y)
    c = np.reshape( y, [ nmol, xdim , ydim ])
    #make empty array to hold rates of change
    xprime = np.empty(c.shape)

    yan = c[0,:,:]
    #R2_5 = c[1,:,:]

    # rates of change for yan
    xprime[0,:,:] = (1-hill(sig(sens,triang),1.0,ney)) * hill(t/N_y,1.0,nny) - y_deg * yan
    # rates of change for R2_5
    #xprime[1,:,:] = (1-hill(yan,1.0,ny25)) * rough * hill(sig(sens,triang),1.0,ns25) - deg25 * R2_5 
    
    
    return xprime.flatten()

#//

distlattice = mk_rand_lattice(xdim,ydim) #set up the cells
triang = package(distlattice) #triangulate

#//

#set up initial conditions
initial = np.zeros((nmol,xdim,ydim))
timerange = range(0,150) #range to solve on

#//

#view initial conditions
plt.clf()
plt.scatter(distlattice[1].flatten(), distlattice[0].flatten(), c = lz, vmin = 0, vmax = 2, s = 1500)
plt.axes().set_aspect('equal')
plt.colorbar()
plt.savefig("fig")

#//

# set up prepattern
sens = np.copy(initial[0])
sens[2,2] = 1
rough = np.copy(initial[0])
rough[3,2],rough[1,2],rough[2,2] = 1,1,1
lz = np.ones((xdim,ydim))
lz[3,2],lz[1,2],lz[2,2],lz[1,3],lz[2,3] = 0,0,0,0,0

#//

#solve it
args = (nmol, triang)
sol = odeint(f,initial.flatten(),timerange,args)
print "done"

resols = [np.reshape(i, [nmol,xdim,ydim]) for i in sol]

#//

#make a pretty plot of results, control the time point and molecule using 'c'
plt.clf()
plt.scatter(distlattice[1].flatten(), distlattice[0].flatten(), c = resols[60][0], vmin = 0, vmax = 1, s = 500)
plt.axes().set_aspect('equal')
plt.colorbar()
plt.savefig("fig")

#//

# plot the concentration of a factor in a single cell over time
plt.clf()
plt.plot(timerange,np.vstack(resols)[:,3,2])
plt.savefig("fig.png")

#//
