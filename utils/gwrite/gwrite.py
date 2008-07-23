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
def gwrite_replace(indent,outfile,newline,flag,type):
    global vardict
    leftbr=newline.split("(")
    buf=leftbr[1].split(",")
    leftbr=newline.split(",")
    leftbr[1]=leftbr[1].rstrip()
    word=leftbr[1].rstrip(")")
    word=leftbr[1].rstrip(");")
    gname=word.strip("\"")
    gname=gname.rstrip("\"")

    if(flag==1):    #fortran
       gname=word.strip("\'")
       gname=gname.rstrip("\'")

    value=vardict[gname].replace("aaaabbbb",buf[0])
    sizestr=xmlparser.sizestr[gname]
    if(len(sizestr)):
       line=sizestr[0]
    else:
       return
    for i in range(1,len(sizestr)):
        line=line+' + '+sizestr[i]
######################
# type:  read=0
#        write=1 
######################
    if (type==0):
       value1=value.replace("adios_op","adios_read")
    elif(type==1):
       value1=value.replace("adios_op","adios_write")
######################
# flag:  c=0
#        fortran=1 
######################
    if(flag==1):
       line=indent+"print(0,*)"+line+'\n'+indent+"call adios_group_size("+buf[0]+", "+str(len(sizestr))+", "+line+')\n'
       value=value1.replace("call",indent+"call")
    elif(flag==2):
       line=indent+'printf(\"%d\\n",'+line+');\n'+indent+"adios_group_size("+buf[0]+", "+str(len(sizestr))+", "+line+');\n'
       value=value1.replace("adios_",indent+"adios_")
    outfile.write(line+value)
    return
 
def gwrite_seach(ifname,ofname):
    #basename,extension=filename.split('.')
    #newfilename=basename+"_g."+extension
    outfile = open(ofname,"w")
    infile = open(ifname,"r")
    indent=[""];
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
             gwrite_replace(indent[0],outfile,line,xmlparser.language_flag,1)
          elif line.find("adios_gread")>=0:
             if (xmlparser.language_flag==1):
                indent=line.split("call") 
                outfile.write(line.replace("call","!call"))
             elif (xmlparser.language_flag==2):
                indent=line.split("adios_gwrite") 
                outfile.write(line.replace("adios_gwrite","//adios_gwrite"))
             gwrite_replace(indent[0],outfile,line,xmlparser.language_flag,0)
          else: 
             outfile.write(line)
    infile.close() 
    outfile.close() 
    return

def main(argv=None):
    global vardict
    arglist=[]
    # os.path.walk can be used to traverse directories recursively 
    # to apply changes to a whole tree of files. 
    os.path.walk("./",callback,arglist)
    os.system("mkdir -p "+ sys.argv[2])
    vardict = xmlparser.getVarlistFromXML(sys.argv[1])
    if (len(vardict)==0):
        print "Fatal: no adios-group in xml file"
        raise SystemExit
    for infile in arglist:
        outfile=sys.argv[2]+"/"+infile
        gwrite_seach(infile,outfile)
    return

if __name__=="__main__":
   main()
