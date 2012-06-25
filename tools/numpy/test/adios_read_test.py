#!/usr/bin/env python
from adios import *
import numpy as np

f = AdiosFile("adios_test.bp")
f.printself()
g = f.group["temperature"]
g.printself()
v = g.var["/temperature"]
v.printself()

print v.read()

f.close()
