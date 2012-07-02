#!/usr/bin/env python
import adios
import numpy as np

adios.init("config.xml")
fd = adios.open("temperature", "adios_test.bp", "w")

NX = 10
size = 1
rank = 0
groupsize =  4 + 4 + 4 + 8 * 1 * NX
t = np.array(range(NX), dtype=np.float64)
adios.set_group_size(fd, groupsize)
adios.write_int(fd, "NX", NX)
adios.write_int(fd, "size", size)
adios.write_int(fd, "rank", rank)
adios.write(fd, "temperature", t)
adios.close(fd)

adios.finalize()

