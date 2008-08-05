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
        sizestr=xmlparser.sizestr[fname]
        if(len(sizestr)):
           line=sizestr[0]
        else:
           return
        for i in range(1,len(sizestr)):
            line=line+' + '+sizestr[i]
        if(xmlparser.language_group_dict[fname]==1):
           wfile= open("gwrite_"+fname+".fh","w") 
           rfile= open("gread_"+fname+".fh","w")
           #line="print(0,*)"+line+'\n'+"call adios_group_size("+str(len(sizestr))+", "+line+')\n'
           line="call adios_group_size(_handle, "+ line +", totalsize,"+ ' comm, err)\n'
        elif(xmlparser.language_group_dict[fname]==2):
           line="adios_group_size(_handle, "+ line +", &totalsize," + ' &comm);\n'
           wfile= open("gwrite_"+fname+".ch","w") 
           rfile= open("gread_"+fname+".ch","w")
        wfile.write(line); 
        wfile.write(vardict[fname].replace("adios_op","adios_write")) 
        rfile.write(vardict[fname].replace("adios_op","adios_read")) 
        wfile.close()
        rfile.close() 

if __name__=="__main__":
   main()
