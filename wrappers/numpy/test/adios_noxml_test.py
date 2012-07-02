#!/usr/bin/env python
from adios import *
import numpy as np

init_noxml()
allocate_buffer(BUFFER_ALLOC_WHEN.NOW, 10)
gid = declare_group ("temperature", "", FLAG.YES)
select_method (gid, "POSIX", "", "")
define_var (gid, "NX", "", DATATYPE.integer, "", "", "")
define_var (gid, "size", "", DATATYPE.integer, "", "", "")
define_var (gid, "rank", "", DATATYPE.integer, "", "", "")
define_var (gid, "temperature", "", DATATYPE.double, "1,NX", "size,NX", "")

fd = open("temperature", "adios_noxml.bp", "w")

NX = 10
size = 1
rank = 0
groupsize =  4 + 4 + 4 + 8 * 1 * NX
t = np.array(range(NX), dtype=np.float64)
set_group_size(fd, groupsize)
write_int(fd, "NX", NX)
write_int(fd, "size", size)
write_int(fd, "rank", rank)
write(fd, "temperature", t)
close(fd)

finalize()

