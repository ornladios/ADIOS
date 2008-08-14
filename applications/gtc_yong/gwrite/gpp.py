#! /usr/bin/python
import os
import sys
import xml.dom.minidom
import dircache
import xmlparser 
from time import sleep

currentdir=os.curdir
def callback(arg, dirname,fnames):
    for filename in fnames: #os.listdir(dirname):
        if(filename.endswith(".F90") or filename.endswith(".f90") or filename.endswith(".F95")
            or filename.endswith(".f95") or filename.endswith(".c")):
           arg.append(filename)

def main(argv=None):
    global vardict
    vardict = xmlparser.getVarlistFromXML(sys.argv[1])
    for fname in vardict:
        line=xmlparser.sizestr[fname]
        if(xmlparser.language_group_dict[fname]==1):
           wfile= open("gwrite_"+fname+".fh","w") 
           rfile= open("gread_"+fname+".fh","w")
           wfile.write(line); 
           rfile.write(line); 
        elif(xmlparser.language_group_dict[fname]==2):
           wfile= open("gwrite_"+fname+".ch","w") 
           rfile= open("gread_"+fname+".ch","w")
           wfile.write(line); 
           rfile.write(line); 
        wfile.write(vardict[fname].replace("adios_op","adios_write")) 
        rfile.write(vardict[fname].replace("adios_op","adios_read")) 
        wfile.close()
        rfile.close() 

if __name__=="__main__":
   main()
