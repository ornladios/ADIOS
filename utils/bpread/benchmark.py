#! /usr/bin/python
import os
import sys
from time import sleep

def main(argv=None):
	#task = [1024, 2048, 4096, 8192, 16384, 32768]
	task = [16, 32]
	line = ""
	line = line + "#!/bin/sh\n"
	line = line + "#PBS -N test\n"
	line = line + "#PBS -j oe\n"
	line0 = line 
	
	for i in task:
		line = line0 
		wfile = open("job"+str(i),"w")
		line = line + "#PBS -l walltime=00:05:00,nodes="+str(i/2)+":ppn=2\n"
		line = line + "mkdir "+str(i)+"\n"	
		line = line + "\ncd $PBS_O_WORKDIR/"+str(i)+"\n"
		line = line + "mpirun -n " + str(i) \
			+ " adios_test_c >& result" + str(i) \
			+ "\n"
		wfile.write(line)
		wfile.close()
	
if __name__=="__main__":
	main()
