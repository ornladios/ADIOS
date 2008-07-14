#!/usr/bin/env python
import sys
import os
import xml.dom.minidom
from time import sleep

items=""
sizeformular=[]

language_flag=0
getsize={}
getsize["byte"]=1
getsize["integer*1"]=1

getsize["short"]=2
getsize["integer*2"]=2

getsize["integer"]=4
getsize["integer*4"]=4

getsize["long"]=8
getsize["integer*8"]=8

getsize["unsigned byte"]=1
getsize["unsigned integer*1"]=1

getsize["unsigned short"]=2
getsize["unsigned integer*2"]=2

getsize["unsigned integer"]=4
getsize["unsigned integer*4"]=4

getsize["unsigned long"]=8
getsize["unsigned integer*8"]=8

getsize["real"]=4
getsize["real*4"]=4
getsize["float"]=4

getsize["unsigned real"]=4
getsize["unsigned real*4"]=4
getsize["unsigned float"]=4

getsize["real*8"]=8
getsize["double"]=8

getsize["unsigned real*8"]=8
getsize["unsigned double"]=8

getsize["complex"]=16
getsize["double complex"]=16

def processvar(node,language_sw,coord_comm,coord_var):
    global sizeformular
    line=""
    dimsname=node.getAttribute("dimensions")
    typename=node.getAttribute("type")
    strsize=str(getsize[typename])
    if (dimsname!=""):
        return line 
    sizeformular.append(strsize)
    varname=str(node.getAttribute("name"))
    if (language_sw==1):
        if(node.getAttribute("gname")!=""):
           varname_g=node.getAttribute("gname")
        else:
           varname_g=varname
        if(node.getAttribute("goffset")!=""):
           varname_g=varname_g+"("+str(node.getAttribute("goffset")) +")"
        if(coord_comm==varname or coord_var==varname): 
           line="call adios_write(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
        else:
           line="call adios_op(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
            
    elif(language_sw==2):
        if(node.getAttribute("gname")!=""):
           varname_g=node.getAttribute("gname")
        else:
           varname_g=varname
        if(node.getAttribute("goffset")!=""):
           varname_g="&"+varname_g+"["+str(node.getAttribute("goffset")) +"]"
        elif(node.getAttribute("dimensions")!="" or node.getAttribute("copy-on-write")=="yes"):
           varname_g=varname_g
        else:
           varname_g="&"+varname_g
        if(coord_comm==varname or coord_var==varname):
           line="call adios_write(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
        else:
           line="call adios_op(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
    return line+'\n'

def processdset(node,language_sw):
    global sizeformular
    line=""
    dimsname=node.getAttribute("dimensions")
    typename=node.getAttribute("type")
    strsize=str(getsize[typename])
    if (dimsname==""):
        return line
    dimsarr=dimsname.split(',');
    for i in range(0,len(dimsarr)):
        if(i==0):
           line=strsize+'*'+'('+dimsarr[i]+')'
        else:
           line=line+'*'+'('+dimsarr[i]+')'
    sizeformular.append(line)
    line=""
    varname=str(node.getAttribute("name"))
    if (language_sw==1):
        if(node.getAttribute("gname")!=""):
           varname_g=node.getAttribute("gname")
        else:
           varname_g=varname
        if(node.getAttribute("goffset")!=""):
           varname_g=varname_g+"("+str(node.getAttribute("goffset")) +")"
        line="call adios_op(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
            
    elif(language_sw==2):
        if(node.getAttribute("gname")!=""):
           varname_g=node.getAttribute("gname")
        else:
           varname_g=varname
        if(node.getAttribute("goffset")!=""):
           varname_g="&"+varname_g+"["+str(node.getAttribute("goffset")) +"]"
        elif(node.getAttribute("dimensions")!="" or node.getAttribute("copy-on-write")=="yes"):
           varname_g=varname_g
        else:
           varname_g="&"+varname_g
        line="adios_op(aaaabbbb,"+"\""+varname+"\","+varname_g+");"
    return line+'\n'
def processnode(nodelist,language_sw,coord_comm,coord_var):
   global items,sizeformular
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName=="var"):
           items=items+processvar(node,language_sw,coord_comm,coord_var)
       subnodelist=node.childNodes
       processnode(subnodelist,language_sw,coord_comm,coord_var)
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName=="var"):
           items=items+processdset(node,language_sw)
       elif (node.nodeType==node.ELEMENT_NODE and node.nodeName=="gwrite"):
           varname=str(node.getAttribute("src"))
           items=items+varname+"\n"
       subnodelist=node.childNodes
       processnode(subnodelist,language_sw,coord_comm,coord_var)

def getVarlistFromXML(xmlFile):
    global items,sizeformular, language_flag 
#    if ( os.path.isfile (xmlFile)):
    variables={}
    doc=xml.dom.minidom.parse(xmlFile)# parse an XML file by name
    idx=0
    group=doc.getElementsByTagName("adios-config")
    language=(group[0].getAttribute("host-language"))
    
    if (language=="Fortran"):
	 language_flag=1
    elif(language=="C" or language=="c" or language=="cpp" or language=="CPP"):
         language_flag=2
    for group in doc.getElementsByTagName("adios-group"):
         items=""
         indent=""
         gname=group.getAttribute("name")
         coord_comm=group.getAttribute("coordination-communicator")
         coord_var=group.getAttribute("coordination-var")
         nodelist=group.childNodes
         processnode(nodelist,language_flag,coord_comm,coord_var)
         variables[str(gname)]=items
    return variables

def main(argv=None):
    if argv is None:
        argv=sys.argv
    vardict = getVarlistFromXML("config.xml")
    #print vardict.keys()
#    for key,value in vardict.items():
#        if(key=="restart"): 
#	        print key,value
if __name__ == "__main__":
    main()
