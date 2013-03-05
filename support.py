#//

#imports
import numpy as np
np.set_printoptions(suppress=True)
from scipy.integrate import odeint
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

#//

#functions

#hill function
def hill(x,k,n):
    return (x ** n) / (x ** n + k ** n)

#function for computing a lattice site from a pair of integers as in lubensky
def mksite(in_x, in_y):
        x = float(in_x) + .5 * float(in_y)
        y = (np.sqrt(3)/2) * float(in_y)
        return [x,y]
    
#function for computing a lattice site from a pair of integers, in a way that makes more sense for the rest of the program
def mksite2(in_x, in_y):
        x = float(in_x) + (.5 * (in_y%2))
        y = (np.sqrt(3)/2) * float(in_y)
        return [x,y]
    
#use that last function to compute the positions on a hexagonal lattice for x and y indicies
def mklattice(x,y):
    pos_array = np.empty([2,x,y])
    it = np.nditer(pos_array[0], flags=['multi_index'])
    while not it.finished:
        pos_array[0,it.multi_index[0],it.multi_index[1]] = mksite2(it.multi_index[0],it.multi_index[1])[0]
        pos_array[1,it.multi_index[0],it.multi_index[1]] = mksite2(it.multi_index[0],it.multi_index[1])[1]
        it.iternext()
    return pos_array

#same as above but with random jitter added to each position
def mk_rand_lattice(x,y):
    pos_array = np.empty([2,x,y])
    it = np.nditer(pos_array[0], flags=['multi_index'])
    while not it.finished:
        pos_array[0,it.multi_index[0],it.multi_index[1]] = mksite2(it.multi_index[0],it.multi_index[1])[0]
        pos_array[0,it.multi_index[0],it.multi_index[1]] +=np.random.normal(0,.15,1)
        pos_array[1,it.multi_index[0],it.multi_index[1]] = mksite2(it.multi_index[0],it.multi_index[1])[1]
        pos_array[1,it.multi_index[0],it.multi_index[1]] +=np.random.normal(0,.15,1)
        it.iternext()
    return pos_array


#find the neighbors of points in a Dealunay triangulation
def find_neighbors(pindex, triang):
    neighbors = list()
    for simplex in triang.vertices:
        if pindex in simplex:
            neighbors.extend([simplex[i] for i in range(len(simplex)) if simplex[i] != pindex])
            '''
            this is a one liner for 'if this simplex contains the point we're interested in,
            extend the neighbors list by appending all the *other* point indices in the simplex
            '''
    #Now we just have to strip out all the duplicate indices and return the list of negihbors:
    return list(set(neighbors))

#packages up a lattice as a triang with neighbors indexed in triang.neighbor_indicies, which is a dictionary for fast lookup
def package(lattice):
    triang = Delaunay(np.reshape(lattice.transpose(), [lattice[0,:,:].size, 2], order = "A"))
    triang.neighbor_indices = dict()
    for index,point in enumerate(triang.points):
        triang.neighbor_indices[index] = find_neighbors(index,triang)
    return triang

#the diffusion opeartor, as defined in lubensky
def diff(conc_array, triang):
    conc_list = conc_array.flatten(order ='A')
    val_array = np.zeros(conc_list.size)
    #an alternate implementation of the diffusion operator, this one using a precomputed triang
    for pindex in range(len(triang.points)):
        for nindex in triang.neighbor_indices[pindex]:
            val_array[pindex] += (1 / np.linalg.norm(triang.points[pindex] - triang.points[nindex])**2) * (conc_list[nindex] - conc_list[pindex])
    return np.reshape(val_array, np.shape(conc_array))

#//
