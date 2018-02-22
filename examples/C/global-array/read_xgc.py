import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import errno
import numpy as np

from threading import Thread

import zmq

def read_pipe_data(pipename, bytes_to_read):
    print('Opening pipe...{0}'.format(pipename) )
    pipe_id  = os.open(pipename, os.O_RDONLY)
    print('Pipe opened...{0}'.format(pipename) )

    while True:
        data = os.read(pipe_id, 1024*16)
        print('Read {0} bytes for {1}'.format(len(data), pipename) )

        global field_data
        global R_data
        global Z_data
        global mesh_data
        global field_bytes
        global R_bytes
        global Z_bytes
        global mesh_bytes
        if (pipename == field):
            field_data = np.append(field_data, np.frombuffer(data, dtype=np.float64))
        elif (pipename == R):
            R_data = np.append(R_data, np.frombuffer(data, dtype=np.float64))
        elif (pipename == Z):
            Z_data = np.append(Z_data, np.frombuffer(data, dtype=np.float64))
        elif (pipename == mesh):
            mesh_data = np.append(mesh_data, np.frombuffer(data, dtype=np.int32))
        else:
            print("pipe name is wrong")

        bytes_to_read = bytes_to_read - len(data)
        if bytes_to_read == 0:
            break

    os.close(pipe_id)
    print("pipe closed...{0}.".format(pipename))
   
def read_pipe_header(pipename):
    pipe_id  = os.open(pipename, os.O_RDONLY)
    print('Pipe opened for header...{0}'.format(pipename) )
    data = os.read(pipe_id, 4)
    os.close(pipe_id)

    global field_bytes
    global R_bytes
    global Z_bytes
    global mesh_bytes
    if (pipename == field):
        field_bytes = np.frombuffer(data, dtype=np.int32)
        print('To read {0} bytes from {1}..'.format(field_bytes, pipename))
    elif (pipename == R):
        R_bytes = np.frombuffer(data, dtype=np.int32)
        print('To read {0} bytes from {1}..'.format(R_bytes, pipename))
    elif (pipename == Z):
        Z_bytes = np.frombuffer(data, dtype=np.int32)
        print('To read {0} bytes from {1}..'.format(Z_bytes, pipename))
    elif (pipename == mesh):
        mesh_bytes = np.frombuffer(data, dtype=np.int32)
        print('To read {0} bytes from {1}..'.format(mesh_bytes, pipename))
    else:
        print("pipe name is wrong")

step = 0
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:5555")

while True:
    print('zmq waiting to recv header.')
    message = socket.recv()
    header = np.frombuffer(message, dtype=np.int32)
    bytes_to_read = header[0]
    bytes_to_read_mesh = header[1]

    print('To read {0} bytes for R/Z/dpot, and {1} bytes for mesh'.format(bytes_to_read, bytes_to_read_mesh))
    field = '/tmp/MdtmManPipes/field'
    R     = '/tmp/MdtmManPipes/R'
    Z     = '/tmp/MdtmManPipes/Z'
    mesh  = '/tmp/MdtmManPipes/mesh'

    field_data = []
    R_data = []
    Z_data = []
    mesh_data = []

    try:
        os.mkfifo(field)
        os.mkfifo(R)
        os.mkfifo(Z)
        os.mkfifo(mesh)
    except OSError as oe: 
        if oe.errno != errno.EEXIST:
            raise

    t1 = Thread(target=read_pipe_data, args=(field,bytes_to_read))
    t1.start()
    t2 = Thread(target=read_pipe_data, args=(R,bytes_to_read))
    t2.start()
    t3 = Thread(target=read_pipe_data, args=(Z,bytes_to_read))
    t3.start()
    t4 = Thread(target=read_pipe_data, args=(mesh,bytes_to_read_mesh))
    t4.start()

    t1.join()
    t2.join()
    t3.join()
    t4.join()

    mesh_data = mesh_data.reshape(mesh_data.size // 3,3)
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.tricontourf(R_data, Z_data, mesh_data, field_data, cmap=plt.cm.jet, levels=np.linspace(-120,120,num=1000))
#plt.plot(R_data, Z_data)
    plt.title('electrostatic potential')
    plt.xlabel('R')
    plt.ylabel('Z')
    figure_name = "foo_" +  str(step) + ".png"
    plt.savefig(figure_name)

    print('filename = {0}'.format(figure_name))
#    os.remove(field)
#    os.remove(R)
#    os.remove(Z)
#    os.remove(mesh)

    step = step + 1
    message = socket.send_string("ok")
