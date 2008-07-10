#!/usr/bin/env python
import sys
import os
import xml.dom.minidom
from time import sleep

items=""
nindent=0
language_flag=0
def processvar(node,language_sw,coord_comm,coord_var):
    line="" 
    varname=str(node.getAttribute("name"))
    if (language_sw==1):
        if(node.getAttribute("gname")!=""):
           varname_dim=node.getAttribute("gname")
        else:
           varname_dim=varname
        if(node.getAttribute("goffset")!=""):
           varname_dim=varname_dim+"("+str(node.getAttribute("goffset")) +")"
        if(coord_comm==varname or coord_var==varname):
           line="call adios_write(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_dim+")"
        else:
           line="call adios_op(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_dim+")"
    elif(language_sw==2):
        if(node.getAttribute("gname")!=""):
           varname_dim=node.getAttribute("gname")
        else:
           varname_dim=varname
        if(node.getAttribute("goffset")!=""):
           varname_dim="&"+varname_dim+"["+str(node.getAttribute("goffset")) +"]"
        elif(node.getAttribute("dimensions")!="" or node.getAttribute("copy-on-write")=="yes"):
           varname_dim=varname_dim
        else:
           varname_dim="&"+varname_dim
        if(coord_comm==varname or coord_var==varname):
           line="adios_write(aaaabbbb,"+"\""+varname+"\","+varname_dim+");"
        else:
           line="adios_op(aaaabbbb,"+"\""+varname+"\","+varname_dim+");"
    return line

def processnode(nodelist,language_sw,coord_comm,coord_var):
   global items
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName=="var"):
           items=items+processvar(node,language_sw,coord_comm,coord_var)+"\n"
       elif (node.nodeType==node.ELEMENT_NODE and node.nodeName=="gwrite"):
           varname=str(node.getAttribute("src"))
           items=items+varname+"\n"
       subnodelist=node.childNodes
       processnode(subnodelist,language_sw,coord_comm,coord_var)
 
def getVarlistFromXML(xmlFile):
    global items 
    global language_flag 
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
