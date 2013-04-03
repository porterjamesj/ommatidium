# data.py
# just a file where i typed in nic's data

import numpy as np

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

progs = np.transpose(progs)
r8s = np.transpose(r8s)
r25s = np.transpose(r25s)
r34s = np.transpose(r34s)

# look at data
'''
plt.plot(progs.transpose()[0],progs.transpose()[1], marker = 'o')
plt.plot(r8s.transpose()[0],r8s.transpose()[1], marker = 'o')
plt.plot(r25s.transpose()[0],r25s.transpose()[1], marker = 'o')
plt.plot(r34s.transpose()[0],r34s.transpose()[1], marker = 'o')
plt.axis([0,300,0,1.5])
plt.savefig("fig.png")
'''
