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
           #basename,extension=filename.split('.')
           arg.append(filename)
def replace(indent,outfile,newline,vardict,flag,type):
    leftbr=newline.split("(")
    buf=leftbr[1].split(",")
    leftbr=newline.split(",")
    leftbr[1]=leftbr[1].rstrip()
    word=leftbr[1].rstrip(")")
    word=leftbr[1].rstrip(");")
    gname=word.strip("\"")
    gname=gname.rstrip("\"")
    value=vardict[gname].replace("aaaabbbb",buf[0])
    if (type==0):
       value1=value.replace("adios_op","adios_read")
    elif(type==1):
       value1=value.replace("adios_op","adios_write")
    if(flag==2):
       value=value1.replace("adios_",indent+"adios_")
    elif(flag==1):
       value=value1.replace("call",indent+"call")
    outfile.write(value)
    return
 
def search_gwrite(ifname,ofname,vardict):
    #basename,extension=filename.split('.')
    #newfilename=basename+"_g."+extension
    outfile = open(ofname,"w")
    infile = open(ifname,"r")
    for line in infile:
          if (line.find("//adios_gwrite")>=0 or line.find("!call adios_gwrite")>=0 
          or line.find("//adios_gread")>=0  or line.find("!call adios_gread")>=0):
             outfile.write(line)
          elif line.find("adios_gwrite")>=0:
             if (xmlparser.language_flag==1): #Fortran
                indent=line.split("call") 
                outfile.write(line.replace("call","!call"))
             elif (xmlparser.language_flag==2):
                indent=line.split("adios_gwrite") 
                outfile.write(line.replace("adios_gwrite","//adios_gwrite"))
             replace(indent[0],outfile,line,vardict,xmlparser.language_flag,1)
          elif line.find("adios_gread")>=0:
             if (xmlparser.language_flag==1):
                indent=line.split("call") 
                outfile.write(line.replace("call","!call"))
             elif (xmlparser.language_flag==2):
                indent=line.split("adios_gwrite") 
                outfile.write(line.replace("adios_gwrite","//adios_gwrite"))
             replace(indent[0],outfile,line,vardict,xmlparser.language_flag,0)
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
    os.system("mkdir -p "+ sys.argv[2])
    vardict = xmlparser.getVarlistFromXML(sys.argv[1])
    for infile in arglist:
        outfile=sys.argv[2]+"/"+infile
        search_gwrite(infile,outfile,vardict)
    return

if __name__=="__main__":
   main()
