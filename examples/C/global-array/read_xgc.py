import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import sys
import errno
import numpy as np
import logging as log
from threading import Thread
import cv2
import zmq
import argparse
import time
import math
import queue as Queue

parser = argparse.ArgumentParser()
parser.add_argument(
    "-v", "--verbose", help="set verbosity",
    action="store_true", dest="verbosity"
    )
parser.add_argument(
    "-m", "--movie", help="make movie",
    action="store_true", dest="ifmovie"
    )
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument(
    "-p", "--port", type=int, 
    help="Set ZMQ port", required=True
    )
args = parser.parse_args()

if args.verbosity:
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
else:
    log.basicConfig(format="%(levelname)s: %(message)s")

log.info("port = {0}".format(args.port))

field_queue = Queue.Queue()
R_queue = Queue.Queue()
Z_queue = Queue.Queue()
mesh_queue = Queue.Queue()

def plot_xgc():
    step = 0
    while True:
        field_data = field_queue.get()
        R_data = R_queue.get()
        Z_data = Z_queue.get()
        mesh_data = mesh_queue.get()

        mesh_data = mesh_data.reshape(mesh_data.size // 3,3)
        plt.figure(step)
        plt.subplot(221)
        plt.gca().set_aspect('equal')
        plt.tricontourf(R_data, Z_data, mesh_data, field_data, cmap=plt.cm.jet, levels=np.linspace(-110,105,num=50))
#plt.plot(R_data, Z_data)
        plt.title('electrostatic potential')
        plt.xlabel('R')
        plt.ylabel('Z')

        plt.subplot(222)
        plt.title('Mesh')
        plt.xlabel('R')
        plt.ylabel('Z', labelpad=-5)
        axes = plt.gca()
        axes.set_ylim([0.2,0.7])
        axes.set_xlim([1.8,2.4])
        plt.triplot(R_data, Z_data, mesh_data, linewidth=0.1)

        grad_data = np.zeros(field_data.size)
        grad_R_data = np.zeros(field_data.size)
        grad_Z_data = np.zeros(field_data.size)
        for m in range(0, mesh_data.size // 3):
            n1 = int(mesh_data[m][0])
            n2 = int(mesh_data[m][1])
            n3 = int(mesh_data[m][2])

            grad_z = field_data[n1] * (Z_data[n2] - Z_data[n3]) + field_data[n2] * (Z_data[n3] - Z_data[n1]) + field_data[n3]* (Z_data[n1] - Z_data[n2])
            grad_r = field_data[n1] * (R_data[n3] - R_data[n2]) + field_data[n2] * (R_data[n1] - R_data[n3]) + field_data[n3]* (R_data[n2] - R_data[n1])
            grad_mag = math.sqrt (grad_z*grad_z +  grad_r*grad_r)

            grad_data[n1] = grad_mag
            grad_data[n2] = grad_mag
            grad_data[n3] = grad_mag

            grad_R_data[n1] = grad_r
            grad_R_data[n2] = grad_r
            grad_R_data[n3] = grad_r

            grad_Z_data[n1] = grad_z
            grad_Z_data[n2] = grad_z
            grad_Z_data[n3] = grad_z

        plt.subplot(223)
        plt.gca().set_aspect('equal')
        plt.tricontourf(R_data, Z_data, mesh_data, grad_R_data, cmap=plt.cm.jet)
#plt.plot(R_data, Z_data)
        plt.title('gradient (R)')
        plt.xlabel('R')
        plt.ylabel('Z')

        plt.subplot(224)
        plt.gca().set_aspect('equal')
        plt.tricontourf(R_data, Z_data, mesh_data, grad_Z_data, cmap=plt.cm.jet)
#plt.plot(R_data, Z_data)
        plt.title('gradient (Z)')
        plt.xlabel('R')
        plt.ylabel('Z')

        plt.tight_layout()

        figure_name = "dpot_" + format(step, '02d') + ".png"
        plt.savefig(figure_name)

        print('Plotting {0}'.format(figure_name))

        step = step + 1


def read_pipe_data(pipename, bytes_to_read):
    log.info('Opening pipe...{0}'.format(pipename))
 
    pipe_id  = os.open(pipename, os.O_RDONLY)

    log.info('Pipe opened...{0}'.format(pipename))

    global field_queue
    global R_queue
    global Z_queue
    global mesh_queue

    while True:
        data = os.read(pipe_id, 1024*16)
        log.info('Read {0} bytes for {1}'.format(len(data), pipename))

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
            log.error("pipe name is wrong")

        bytes_to_read = bytes_to_read - len(data)
        if bytes_to_read == 0:
            break

    if (pipename == field):
        field_queue.put(field_data)
        log.info("Put field into the queue.")
    elif (pipename == R):
        R_queue.put(R_data)
        log.info("Put R into the queue.")
    elif (pipename == Z):
        Z_queue.put(Z_data)
        log.info("Put Z into the queue.")
    elif (pipename == mesh):
        mesh_queue.put(mesh_data)
        log.info("Put mesh into the queue.")
    else:
        log.error("pipe name is wrong")

    os.close(pipe_id)
    log.info("pipe closed...{0}.".format(pipename))
   
def pdist(pt1, pt2):
    x = pt1[0] - pt2[0]
    y = pt1[1] - pt2[1]
    return math.sqrt(math.pow(x, 2) + math.pow(y, 2))

p1 = Thread(target=plot_xgc, args=())
p1.start()

step = 0
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:{0}".format(args.port))

print("XGC analysis is ready to process data.")
while True:
    try:
        message = socket.recv()
    except KeyboardInterrupt:
        if args.ifmovie:
            print("All plots are achived into avi.")
            os.system("ffmpeg -f image2 -framerate 1 -i dpot_%02d.png -s:v 1280x720 -c:v mpeg4 -crf 20 -pix_fmt yuv420p dpot.mp4")
            os.system("mplayer -vo x11 dpot.mp4")
        socket.close()
        context.term()
        sys.exit()

    header = np.frombuffer(message, dtype=np.int32)
    bytes_to_read = header[0]
    bytes_to_read_mesh = header[1]

    log.info("To read {0} bytes for R/Z/dpot, and {1} bytes for mesh".format(bytes_to_read, bytes_to_read_mesh))
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

    step = step + 1
    message = socket.send_string("ok")
