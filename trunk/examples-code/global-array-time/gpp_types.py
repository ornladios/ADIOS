import sys
import os

adios_size={}
adios_size["byte"]=1
adios_size["integer*1"]=1

adios_size["short"]=2
adios_size["integer*2"]=2

adios_size["integer"]=4
adios_size["integer*4"]=4

adios_size["string"]=1
adios_size["long"]=8
adios_size["integer*8"]=8

adios_size["unsigned byte"]=1
adios_size["unsigned integer*1"]=1

adios_size["unsigned short"]=2
adios_size["unsigned integer*2"]=2

adios_size["unsigned integer"]=4
adios_size["unsigned integer*4"]=4

adios_size["unsigned long"]=8
adios_size["unsigned integer*8"]=8

adios_size["real"]=4
adios_size["real*4"]=4
adios_size["float"]=4

adios_size["unsigned real"]=4
adios_size["unsigned real*4"]=4
adios_size["unsigned float"]=4

adios_size["real*8"]=8
adios_size["double"]=8

adios_size["unsigned real*8"]=8
adios_size["unsigned double"]=8

adios_size["complex"]=8
adios_size["double complex"]=16
