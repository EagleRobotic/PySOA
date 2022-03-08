import numpy as np
from IPython import embed

k = 5
M = np.eye(6*k) #Mass matrix
m = np.ones(5) #masses
m = m* 0.17

pgraph = np.array([1,2,3,4,5])
print(len(pgraph))

i,j = 0,1

link_len = np.tile([1, 0, 0], (k+1,1))
link_len = 0.2 * link_len
lci = 0.5*link_len

I = np.array([[8.5e-6, 0, 0],[0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2), 0],[0, 0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2)],
            [8.5e-6, 0, 0],[0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2), 0],[0, 0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2)],
            [8.5e-6, 0, 0],[0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2), 0],[0, 0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2)],
            [8.5e-6, 0, 0],[0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2), 0],[0, 0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2)],
            [8.5e-6, 0, 0],[0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2), 0],[0, 0, float(5.72e-4)+np.power(np.linalg.norm(lci[0,0:3]),2)]])


#embed()

for i in range (k):
    M[6*i:6*i+3,6*i:6*i+3]=I[3*i:3*(i+1),0:3]
    M[6*i+3:6*(i+1),6*i+3:6*(i+1)] =m[i]*M[6*i+3:6*(i+1),6*i+3:6*(i+1)]