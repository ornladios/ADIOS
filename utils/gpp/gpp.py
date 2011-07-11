#!/usr/bin/env python
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


def checkXML (config_file, path):
    if path == '':
        adios_lint = 'adios_lint '
    else:
        adios_lint = path + '/adios_lint '
    rv = os.system (adios_lint + config_file)
    if rv == 0:
        return 'success'
    elif rv == 32512:  # System unable to find adios_lint command
        print "Unable to find adios_lint. Proceeding with code generation."
        return 'success'
    else:
        print "gpp.py failed."
        return 'failure'


def main(argv=None):

    if len (sys.argv) != 2:
        print 'Usage: gpp.py <config file>\n'
        return 1

    check_val = checkXML (sys.argv[1], os.path.dirname(sys.argv[0])) 
    if check_val != 'success':
        return 1

    global vardict
    vardict = xmlparser.getVarlistFromXML(sys.argv[1])
    for fname in vardict:
        #line=xmlparser.sizestr[fname]
        line = ""
        if(xmlparser.language_group_dict[fname]==1):
           wfile= open("gwrite_"+fname+".fh","w") 
           rfile= open("gread_"+fname+".fh","w")
	   rfile.write("adios_groupsize = 0\n"); 
	   rfile.write("adios_totalsize = 0\n"); 
        elif(xmlparser.language_group_dict[fname]==2):
           wfile= open("gwrite_"+fname+".ch","w") 
           rfile= open("gread_"+fname+".ch","w")
	   rfile.write("adios_groupsize = 0;\n"); 
	   rfile.write("adios_totalsize = 0;\n"); 
        wfile.write(line); 
        rfile.write(line); 
        linerw = vardict[fname]
        wfile.write(linerw[0])
        rfile.write(linerw[1]) 
        wfile.close()
        rfile.close() 

if __name__=="__main__":
   main()
