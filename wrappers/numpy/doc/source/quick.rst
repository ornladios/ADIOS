.. _quick:

Quick Start Guide
=================

An Adios bp file is a container of variables and attributes, which can
be either scalars or arrays. With the Adios Python wrapper, one can
access them as NumPy's array.

Reading data
------------

In this quick start guide, we assume we have an Adios bp file (name:
adios_test.bp) contains the following three variables and one attribute:
::
   $ bpls -lva adios_test.bp
   integer  NX                        scalar = 10
   integer  size                      scalar = 2
   double   temperature               {2, 10} = 0 / 19 / 9.5 / 5.76628
   string   /temperature/description  attr   = "Global array written from 'size' processes"


Let's start with importing Adios Python module as follows:

>>> import adios as ad

.. note:: A parallel version of Adios Python module is also available,
          called adios_mpi. We will discuss details later.


Then, open the Adios bp file (adios_test.bp) and briefly check its
contents:

>>> f = ad.file('adios_test.bp')
>>> f
AdiosFile (path='adios_test.bp', nvars=3, var=['NX', 'temperature', 'size'], nattrs=1, attr=['/temperature/description'], current_step=0, last_step=0, file_size=1549)

Now let's read a scalar variable, 'NX'. We can access each variable by
using Python's dictionary-style interface:

>>> v = f['NX']
>>> v
AdiosVar (varid=0, type=dtype('int32'), ndim=0, dims=(), nsteps=1)
>>> v.read()
array(10, dtype=int32)

Equivalently, we can use Numpy-style interface:

>>> f['NX'][...]
array(10, dtype=int32)


Then, let's read a multi-dimensional array, 'temperature'.

>>> f['temperature'][:]
array([[  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.],
       [ 10.,  11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.]])

We can read by slices, for example:

>>> f['temperature'][0:2,0:5]
array([[  0.,   1.,   2.,   3.,   4.],
       [ 10.,  11.,  12.,  13.,  14.]])

We can use the most of the NumPy's slice syntax.

Attribute reading is similar:

>>> at = f.attr['/temperature/description']
>>> at
AdiosAttr (name='/temperature/description', type=dtype('S43'))
>>> at.value
array(["Global array written from 'size' processes"],
      dtype='|S43')

Unless attribute's name is not conflict with any variable name in the
file, which is totally fine with Adios, we can access through the
dictionary-style:

>>> f['/temperature/description'].value

.. note:: While variables are read by "read" function or slice
          interface, attributes are accessed through "value" property.

Writing data
------------

Now, we will show how we can create the Adios BP file used in the
previous section:
::

   $ bpls -lva adios_test.bp
   integer  NX                        scalar = 10
   integer  size                      scalar = 2
   double   temperature               {2, 10} = 0 / 19 / 9.5 / 5.76628
   string   /temperature/description  attr   = "Global array written from 'size' processes"


First, we load necessary modules and prepare our Numpy data to save:

>>> import adios as ad
>>> import numpy as np

>>> NX = 10
>>> size = 2
>>> t = np.array(range(NX*size), dtype=np.float64)
>>> tt = t.reshape((size, NX))

I.e., we have two scalar variables (NX and size) and one 2-D array (tt).

Then, we initialize Adios and specify a buffer size (10MB) which Adios
can work with:

>>> ad.init_noxml()
>>> ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);

Then, we give a file name to create and specify a group with Adios method:

>>> fw = ad.writer(fname)
>>> fw.declare_group('group', method='POSIX1')

"POSIX1" is one of many Adios's write methods. Others are "MPI",
"MPI_AGGREGATE", "FLEXPATH", "DATASPACES", etc. More detailed
descriptions are in the Adios manual.

Now, we assign our values:

>>> fw['NX'] = NX
>>> fw['size'] = size
>>> fw['temperature'] = tt

To write an attribute, we can do as follows:

>>> fw.attr['/temperature/description'] = "Global array written from 'size' processes"

Finally, we let Adios to write a file by calling "close"

>>> fw.close()
