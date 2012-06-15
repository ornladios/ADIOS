#!/usr/bin/env python
from adios import *
import numpy as np

f = AdiosFile("adios_test.bp")
f.printself()
g = f.group("temperature")
g.printself()
v = g.variable("temperature")
v.printself()

print v.read()

v.close()
g.close()
f.close()
