#!/usr/bin/env python
import sys
import os
import xml.dom.minidom
from time import sleep

items=""
nindent=0
def processvar(node,language):
    varname=str(node.getAttribute("name"))
    if (language=="Fortran"):
        if(node.getAttribute("gname")!=""):
           varname_dim=node.getAttribute("gname")
        else:
           varname_dim=varname
        if(node.getAttribute("goffset")!=""):
           varname_dim=varname_dim+"("+str(node.getAttribute("goffset")) +")"
        #else:
  	#   varname_dim=varname
    elif(language=="C" or language=="c" or language=="cpp" or language=="CPP"):
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
    line="call adios_write(gid,"+"\""+varname+"\"//char(0),"+varname_dim+")"
    return line


def processnode(nodelist,language):
   global items
#   global nindent
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName=="var"):
#           indent=""
#           for i in range(nindent): 
#               indent=indent + "\t"
#           items=items+"\n"+indent+processvar(node,language)
           items=items+processvar(node,language)+"\n"
       elif (node.nodeType==node.ELEMENT_NODE and node.nodeName=="condition"):
           varname=str(node.getAttribute("expression"))
#           indent=""
#           varname=str(node.getAttribute("begin_expression"))
#           if (varname!=""):
#                for i in range(nindent):
#                    indent=indent+"\t"
#                nindent=nindent+1
#           else:
#                varname=str(node.getAttribute("end_expression"))
#                if (varname!=""):
#                    nindent=nindent-1
#                    for i in range(nindent):
#                        indent=indent+"\t"
#           items=items+"\n"+indent+varname
           items=items+varname+"\n"
       subnodelist=node.childNodes
       processnode(subnodelist,language)
 
def getVarlistFromXML(xmlFile):
#    print xmlFile
    global items 
#    if ( os.path.isfile (xmlFile)):
    variables={}
    #items={}
    doc=xml.dom.minidom.parse(xmlFile)# parse an XML file by name
    idx=0
    group=doc.getElementsByTagName("adios-config")
    language=(group[0].getAttribute("host-language"))
    for group in doc.getElementsByTagName("adios-group"):
         items=""
         indent=""
         gname=group.getAttribute("name")
         nodelist=group.childNodes
         #for node in nodelist:
         processnode(nodelist,language)
         variables[str(gname)]=items
         #print items
    return variables
#	     if (node.nodeType==node.ELEMENT_NODE and node.nodeName=="var"):
#		 items=items+"\n"+processvar(node,language)
#	     elif (node.nodeType==node.ELEMENT_NODE and node.nodeName=="condition"):
#                varname=str(node.getAttribute("expression"))
#                items=items+"\n"+varname
#                subnodelist=node.childNodes
#                for subnode in subnodelist:
#         	    if (subnode.nodeType==subnode.ELEMENT_NODE and subnode.nodeName=="var"):
#		        items=items+"\n\t"+processvar(node,language)

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
