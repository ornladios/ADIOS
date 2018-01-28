import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import errno
import numpy as np

from threading import Thread

def myfunc(pipename):
    print('Opening pipes...{0}'.format(pipename))
    pipe_id  = os.open(pipename, os.O_RDONLY)
    print('Pipe opened...{0}'.format(pipename) )
    data = os.read(pipe_id, 1024*1024*32)
    print(len(data))
    os.close(pipe_id)

    global field_data
    global R_data
    global Z_data
    global mesh_data
    if (pipename == field):
        field_data = np.frombuffer(data, dtype=np.float64)
    elif (pipename == R):
        R_data = np.frombuffer(data, dtype=np.float64)
    elif (pipename == Z):
        Z_data = np.frombuffer(data, dtype=np.float64)
    elif (pipename == mesh):
        print('converting mesh data')
        mesh_data = np.frombuffer(data, dtype=np.int32)
    else:
        print("pipe name is wrong")

    print("pipe closed...{0}".format(pipename))
   

field = '/tmp/MdtmManPipes/field'
R     = '/tmp/MdtmManPipes/R'
Z     = '/tmp/MdtmManPipes/Z'
mesh  = '/tmp/MdtmManPipes/mesh'

try:
    os.mkfifo(field)
    os.mkfifo(R)
    os.mkfifo(Z)
    os.mkfifo(mesh)
except OSError as oe: 
    if oe.errno != errno.EEXIST:
        raise

t1 = Thread(target=myfunc, args=(field,))
t1.start()
t2 = Thread(target=myfunc, args=(R,))
t2.start()
t3 = Thread(target=myfunc, args=(Z,))
t3.start()
t4 = Thread(target=myfunc, args=(mesh,))
t4.start()

t1.join()
t2.join()
t3.join()
t4.join()

"""
R_data = R_data.astype(np.float64)
Z_data = Z_data.astype(np.float64)
field_data = field_data.astype(np.float64)
mesh_data = mesh_data.astype(np.int32)
"""
mesh_data = mesh_data.reshape(mesh_data.size // 3,3)
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(R_data, Z_data, mesh_data, field_data, cmap=plt.cm.jet, levels=np.linspace(-120,120,num=1000))
#plt.plot(R_data, Z_data)
plt.title('electrostatic potential')
plt.xlabel('R')
plt.ylabel('Z')
plt.savefig('foo.png')

os.remove(field)
os.remove(R)
os.remove(Z)
os.remove(mesh)

