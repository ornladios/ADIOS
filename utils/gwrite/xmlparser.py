#!/usr/bin/env python
import sys
import os
import xml.dom.minidom
import gwrite_types 
from time import sleep

items=""
sizeformular=[]
sizestr={}
getsize=gwrite_types.adios_size
language_flag=0
groupdict={}
def processvar(node,language_sw,coord_comm,coord_var):
    global sizeformular
    global groupdict 
    line=""
########################################################          
# Modified Code: Case Insensitive
########################################################          
    attkeys=node.attributes.keys()
    varname=""
    typename="" 
    pathname="" 
    dimsname=""
    varname_g=""
    for akey in attkeys:
        akeystr=str(akey).lower()
        if(akeystr=="dimensions"):
           dimsname=node.attributes[akey].value
        elif(akeystr=="type"):
           typename=node.attributes[akey].value
        elif(akeystr=="path"):
           pathname=node.attributes[akey].value
        elif(akeystr=="name"):
           varname=node.attributes[akey].value
        elif(akeystr=="gname"):
           varname_g=node.attributes[akey].value
    if(varname=="" or typename==""):
       print "Warning: empty varname or type for adiosgroup: "+node.nodeName
    if (dimsname!=""):
        return line
    #print varname,dimsname,typename,varname_g 
    strsize=str(getsize[typename])
    sizeformular.append(strsize)
    if(varname_g==""):
       varname_g=varname
    if(groupdict.has_key(varname)):
       if(pathname[len(pathname)-1]!='/'):
          varname=pathname+'/'+varname 
       groupdict[varname]=varname_g
    else:
       groupdict[varname]=varname_g

########################################################          
# Original Code: Case Sensitive
########################################################          
    #dimsname=node.getAttribute("dimensions")
    #typename=node.getAttribute("type")
    #varname=str(node.getAttribute("name"))
    #if(node.getAttribute("gname")!=""):
    #   varname_g=node.getAttribute("gname")
    #else:
    #   varname_g=varname
########################################################          
           
    if (language_sw==1):
        if(coord_comm==varname or coord_var==varname): 
           line="call adios_write(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
        else:
           line="call adios_op(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
    elif(language_sw==2):
        if(coord_comm==varname or coord_var==varname):
           line="adios_write(aaaabbbb,"+"\""+varname+"\","+"&"+varname_g+");"
        else:
           for c in varname_g:
               if(c=='+' or c=='-' or c=='*' or c=='/' or c=='^' or c=='%'):
                  print "Fatal: var --"+varname_g+"-- cannot be written"
                  raise SystemExit
           line="adios_op(aaaabbbb,"+"\""+varname+"\","+"&"+varname_g+");"
    return line+'\n'

def processattr(node,language_sw):
    line=""
    return line

def processdset(node,language_sw):
    global sizeformular
    global groupdict
    line=""

########################################################          
# Modified Code: Case Insensitive
########################################################          
    attkeys=node.attributes.keys()
    varname=""
    typename="" 
    dimsname=""
    varname_g=""
    copyflag=""
    for akey in attkeys:
        akeystr=str(akey).lower()
        if(akeystr=="dimensions"):
           dimsname=node.attributes[akey].value
        elif(akeystr=="type"):
           typename=node.attributes[akey].value
        elif(akeystr=="name"):
           varname=node.attributes[akey].value
        elif(akeystr=="gname"):
           varname_g=node.attributes[akey].value
        elif(akeystr=="copy-on-write"):
           copyflag=node.attributes[akey].value
    if (dimsname==""):
        return line 
    if(varname_g==""):
       varname_g=varname
########################################################          
# Original Code: Case Sensitive
########################################################          
#    dimsname=node.getAttribute("dimensions")
#    typename=node.getAttribute("type")
########################################################          
    strsize=str(getsize[typename])
    dimsarr=dimsname.split(',');
    for dimsele in dimsarr:
           
        if(dimsele==dimsarr[0]):
           if(dimsele.isdigit()):
              line=strsize+'*'+'('+dimsele+')'
           else:
              line=strsize+'*'+'('+groupdict[dimsele]+')'
        else:
           if(dimsele.isdigit()):
              line=line+'*'+'('+dimsele+')'
           else:
              line=line+'*'+'('+groupdict[dimsele]+')'
   
    if(str(node.getAttribute("copy-on-write")).lower()=="yes"):
           line ='2*'+'('+line+')'
    
    sizeformular.append(line)
    line=""
    if (language_sw==1):
        line="call adios_op(aaaabbbb,"+"\""+varname+"\"//char(0),"+varname_g+")"
    elif(language_sw==2):
        line="adios_op(aaaabbbb,"+"\""+varname+"\","+varname_g+");"
    return line+'\n'

def processnode(nodelist,language_sw,coord_comm,coord_var):
   global items,sizeformular
# process scalar-type var first
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="var"):
           items=items+processvar(node,language_sw,coord_comm,coord_var)
       subnodelist=node.childNodes
       processnode(subnodelist,language_sw,coord_comm,coord_var)
# process dataset-type var  
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="var"):
           items=items+processdset(node,language_sw)
       elif (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="attribute"):
           items=items+processattr(node,language_sw)
       elif (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="gwrite"):
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
    group=doc.childNodes
    if(group.length!=1):
       print "Fatal: adios-config should be one and only one!"
       raise SystemExit
    root=group[0]
    if(root.tagName.lower()!="adios-config"):
       print "Fatal: wrong root element, switch to adios-config"  
       raise SystemExit
    #group=doc.getElementsByTagName("adios-config")
    attnum=len(root.attributes.keys())
    if (attnum==0) :
        language_flag=1
    elif(attnum==1):
        attkeys=root.attributes.keys()
        if(str(attkeys[0]).lower()=="host-language"):
           language=group[0].getAttribute(attkeys[0])
        if (language.lower()=="fortran"):
    	    language_flag=1
        elif(language.lower()=="c" or language.lower()=="cpp"):
            language_flag=2
        else:
            language_flag==-1
            print "Fatal: Unknown language supported!"
            raise SystemExit
    gname=""
    coord_comm=""
    coord_var=""
    for group in root.childNodes:
      if (group.nodeType==group.ELEMENT_NODE): 
        if (group.tagName.lower()=="adios-group"): 
            items=""
            indent=""
            attkeys=group.attributes.keys()
            for akey in attkeys:
                akeystr=str(akey).lower()
                if (akeystr=="name"):
                    gname=group.attributes[akey].value
                elif (akeystr=="coordination-communicator"):
                    coord_comm=group.attributes[akey].value
                elif (akeystr=="coordination-var"):
                    coord_var=group.attributes[akey].value
                else:
                    print "Warning: Unknown attribute --"+str(akey)+"-- for adios-group: "+gname
            if (gname ==""):
                print "Fatal: name for adios-group is required!"
                raise SystemExit
            nodelist=group.childNodes
            processnode(nodelist,language_flag,coord_comm,coord_var)
            groupdict={}
            variables[str(gname)]=items
            sizestr[str(gname)] = sizeformular    
            sizeformular=[]
    return variables

def processgroup(group,gname,coord_comm,coord_var):
    attkeys=group.attributes.keys()
    for akey in attkeys:
        akeystr=str(akey).lower()
        if (akeystr=="name"):
            gname=group.attributes[akey].value
        elif (akeystr=="coordination-communicator"):
            coord_comm=group.attributes[akey].value
        elif (akeystr=="coordination-var"):
            coord_var=group.attributes[akey].value
        else:
            print "Warning: Unknown attribute --"+str(akey)+"-- for adios-group: "+gname
    
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
