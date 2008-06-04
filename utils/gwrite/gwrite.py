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
        if(filename.endswith(".F90") or filename.endswith(".f90") or filename.endswith(".F95") or filename.endswith(".f95")):
           #basename,extension=filename.split('.')
           arg.append(filename)
def replace(outfile,newline,vardict):
    leftbr=newline.split("(")
    buf=leftbr[1].split(",")
    leftbr=newline.split(",")
    leftbr[1]=leftbr[1].rstrip()
    word=leftbr[1].rstrip(")")
    word=leftbr[1].rstrip(");")
    gname=word.strip("\"")
    gname=gname.rstrip("\"")
    value=vardict[gname].replace("gid",buf[0])
    outfile.write(value)
    return
 
def search_gwrite(ifname,ofname,vardict):
    #basename,extension=filename.split('.')
    #newfilename=basename+"_g."+extension
    outfile = open(ofname,"w")
    infile = open(ifname,"r")
    for line in infile:
          if (line.find("//adios_gwrite")>=0 or line.find("!call adios_gwrite")>=0):
              return
          if line.find("adios_gwrite")>=0:
             newline=line.replace("adios_gwrite","adios_write")
             replace(outfile,newline,vardict)
          else: 
             outfile.write(line)
    infile.close() 
    outfile.close() 
    return

def main(argv=None):
    arglist=[]
    # os.path.walk can be used to traverse directories recursively 
    # to apply changes to a whole tree of files. 
    os.path.walk("./",callback,arglist)
    os.system("mkdir -p "+ sys.argv[1])
    vardict = xmlparser.getVarlistFromXML("config.xml")
    for infile in arglist:
        outfile=sys.argv[1]+"/"+infile
#        print outfile
        search_gwrite(infile,outfile,vardict) 
    return

if __name__=="__main__":
   main()

