#!/usr/bin/env python
import sys
import os
import xml.dom.minidom
import gpp_types 
from time import sleep

sizeformular=[]
sizestr={}
getsize=gpp_types.adios_size
var_size_dict={}

def processvar(node,language_sw,coord_comm,coord_var,time_var):
    global sizeformular
    global var_gname_dict 
    global var_size_dict 
    line=""
    liner=""
    linew=""
    linerw={}
########################################################          
# Modified Code: Case Insensitive
########################################################          
    attkeys=node.attributes.keys()
    varname=""
    typename="" 
    pathname="" 
    dimsname=""
    gwname=""
    grname=""
    readyn = True 
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
        elif(akeystr=="gwrite"):
           gwname=node.attributes[akey].value
        elif(akeystr=="gread"):
           grname=node.attributes[akey].value
        elif(akeystr=="read"):
           if(node.attributes[akey].value=='no'):
              readyn = False
    if(varname=="" or typename==""):
       print "Warning: empty varname or type for adiosgroup: "+node.nodeName
    if(gwname==""):
       gwname=varname
    if(grname==""):
       grname=varname
    for c in grname:
       if (c == '+' or c == '-' or c == '*' or c == '/' or c == '^' or c == '**'):
          readyn = False
    if(var_gname_dict.has_key(varname)):
       if(pathname!="" and pathname[len(pathname)-1]!='/'):
          varname=pathname+'/'+varname
       var_gname_dict[varname] = gwname 
    else:
       var_gname_dict[varname] = gwname 

# Add var-size-mapping   
    line="" 
    if (dimsname==""):
       if (typename=="string"):
           if (language_sw==1):
               line="len_trim("+gwname+")" 
               liner="len("+grname+")" 
           else:        
               line="strlen("+gwname+")" 
               # FIXME: we do not know how much space do we have in the 
               #  variable when reading a string into it. 
               #  we just give a stupidly big number (2GB) here and if the
               #  space allocated is less then what is read in,
               #  that is considered to be the responsibility of the programmer
               liner="2147483648" 
       else:
           line=str(getsize[typename])
           if (readyn==True):
                liner=line 
    else:
       str_varsize=str(getsize[typename])
       dimsarr=dimsname.split(',');
       for dimsele in dimsarr:
         if(dimsele==dimsarr[0]):
            if(dimsele.isdigit()):
               line=str_varsize+' * '+'('+dimsele+')'
            elif(var_gname_dict.has_key(dimsele)):
               line=str_varsize+' * '+'('+var_gname_dict[dimsele]+')'
            else:
               line=str_varsize
         else:
            if(dimsele.isdigit()):
               line = line+' * '+'('+dimsele+')'
            elif(dimsele == time_var):
               line = line
            elif(var_gname_dict.has_key(dimsele)):
               line = line + ' * (' + var_gname_dict[dimsele]+')'
            else: 
               line = line
       if (readyn==True):
           liner = line
    sizeformular.append(line)    
    var_size_dict[varname]=line

    if (language_sw==1):
        linew = "call adios_write (adios_handle, "   \
              + "\""+varname+"\", "         \
              + gwname +", adios_err)\n"         
        if (readyn): 
           liner = "adios_buf_size = "+liner                 \
                 + "\ncall adios_read (adios_handle, " \
                 + "\"" + varname                      \
                 + "\", " + grname            \
                 + ", adios_buf_size, adios_err)\n"
    elif(language_sw==2):
        if (dimsname==""):
           if(typename=="string"): 
              linew = "adios_write (adios_handle, "          \
                 + "\"" + varname + "\", "                   \
                 + gwname + ");\n"                         
           else: 
              linew = "adios_write (adios_handle, "          \
                 + "\"" + varname + "\", &"                  \
                 + gwname + ");\n"                         
           if (readyn):
             if(typename=="string"): 
               liner = "adios_buf_size = "+liner                  \
                     + ";\nadios_read (adios_handle, "      \
                     + "\"" + varname + "\", "             \
                     + grname + ", adios_buf_size);\n"
             else: 
               liner = "adios_buf_size = "+liner                  \
                     + ";\nadios_read (adios_handle, "      \
                     + "\"" + varname + "\", &"             \
                     + grname + ", adios_buf_size);\n"
        else: 
           linew = "adios_write (adios_handle, "         \
                 + "\"" + varname + "\", "               \
                 + gwname + ");\n"                         
           if (readyn): 
              liner = "adios_buf_size = " +liner                 \
                     + ";\nadios_read (adios_handle, " \
                     + "\"" + varname + "\", "               \
                     + grname + ", adios_buf_size);\n"
    linerw=[linew,liner]
    return linerw

def processattr(node, language_sw):
    global sizeformular
    global var_gname_dict
    attkeys=node.attributes.keys()
    attrname=""
    varname=""
    pathname=""
    valuename=""
    typename="" 
    for akey in attkeys:
        akeystr=str(akey).lower()
        if ( akeystr == 'name'):
           attrname = node.attributes[akey].value 
        elif ( akeystr == 'path'):
           pathname = node.attributes[akey].value 
        elif ( akeystr == 'value'):
           valuename = node.attributes[akey].value 
        elif (akeystr == 'type'):
           typename = node.attributes[akey].value 
        elif ( akeystr == 'var'):
           varname=node.attributes[akey].value
        else:
           print "\""+pathname+attrname+"\" is not declared before the attribute \"" \
                     +node.attributes["name"].value+"\""
    if (valuename != ''):
       if (typename != ''):
           var_gname_dict[pathname+attrname] = valuename
           var_gname_dict[attrname] = valuename
    elif (varname != ""):
       var_gname_dict[attrname] = var_gname_dict[varname]
         
def processnode(nodelist,language_sw,coord_comm,coord_var,time_var):
   global items,sizeformular
   for node in nodelist:
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="var"):
           itemsrw=(processvar(node,language_sw,coord_comm,coord_var,time_var))
           items[0] = items[0] + itemsrw[0]
           items[1] = items[1] + itemsrw[1]
       if (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="attribute"):
           processattr(node,language_sw)
       elif (node.nodeType==node.ELEMENT_NODE and node.nodeName.lower()=="gwrite"):
           srccode =str(node.getAttribute("src"))
           items [0] = items[0] + srccode + "\n"
           items [1] = items[1] + srccode + "\n"
       subnodelist=node.childNodes
       processnode(subnodelist,language_sw,coord_comm,coord_var,time_var)

def getVarlistFromXML(xmlFile):
    global items,sizeformular, language_group_dict,var_gname_dict,var_size_dict 
    var_gname_dict={}
    items={}
    items=["",""]
    language_group_dict={}
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
    time_var = ""
    coord_var = ""
    coord_comm = "MPI_COMM_WORLD"
    for group in root.childNodes:
      if (group.nodeType==group.ELEMENT_NODE): 
        if (group.tagName.lower()=="adios-group"): 
            indent=""
            attkeys=group.attributes.keys()
            glanguage=language_flag
            for akey in attkeys:
                akeystr=str(akey).lower()
                if (akeystr=="name"):
                    gname=group.attributes[akey].value
                elif (akeystr == "coordination-communicator"):
                    coord_comm=group.attributes[akey].value
                elif (akeystr == "coordination-var"):
                    coord_var=group.attributes[akey].value
                elif (akeystr=="host-language"):
                    language=group.attributes[akey].value
                    if (language.lower()=="fortran"):
                               glanguage=1
                    elif(language.lower()=="c" or language.lower()=="cpp"):
                        glanguage=2
                    else:
                        glanguage==-1
                        print "Fatal: Unknown language supported!"
                        raise SystemExit
                elif (akeystr == "time-index"):
                    time_var = group.attributes[akey].value
                else:
                    print "Warning: Unknown attribute --"+str(akey)+"-- for adios-group: "+gname
            if (gname ==""):
                print "Fatal: name for adios-group is required!"
                raise SystemExit

            nodelist=group.childNodes
             
            #line = ""
            #if(glanguage==1):
            #        if (coord_var!= ""):
            #            line = line                              \
            #                + "call adios_write (adios_handle, "\
            #                + "\""+coord_var+"\"//char(0), "    \
            #                + coord_var +")\n"         
            #        if (time_var != ""):
            #             line = line                              \
            #                  + "call adios_write (adios_handle, "\
            #                  + "\""+time_var+"\"//char(0), "     \
            #                  + time_var +")\n"        
            #elif(glanguage):
            #        if (coord_var!= ""):
            #            line = line                              \
            #                  + "adios_write (adios_handle, "\
            #                  + "\""+coord_var+"\", &"    \
            #                  + coord_var +");\n"         
            #    if (time_var != ""):
            #            line = line                              \
            #                  + "adios_write (adios_handle, "\
            #                  + "\""+time_var+"\", &"     \
            #                  + time_var +");\n"
            #items[0] = items[0]+line
            items=["",""]
            processnode(nodelist,glanguage,coord_comm,coord_var,time_var)

            if(glanguage==1):
               if (coord_comm == ''):
                  line="call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)\n"
               else:
                  line="call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)\n"
            else:
               if (coord_comm == ''):
                   line="adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);\n"
               else:
                   line="adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);\n"
            items[0]=line+items[0]
            items[1]=line+items[1]
            var_gname_dict={}
            line = sizeformular[0]
            if (glanguage == 1): 
                for i in range(1,len(sizeformular)):
                    line =line +' &\n                + '+sizeformular[i]
                items[0]="adios_groupsize = " + line+"\n"+items[0]
            elif (glanguage == 2): 
                for i in range(1,len(sizeformular)):
                    line =line +' \\\n                + '+sizeformular[i]
                items[0]="adios_groupsize = " + line+";\n"+items[0]

            language_group_dict[gname]=glanguage
             
            variables[str(gname)]=items
            sizestr[str(gname)] = line
            sizeformular=[]
    return variables

def main(argv=None):
    if argv is None:
        argv=sys.argv
    vardict = getVarlistFromXML("config.xml")
    #print vardict.keys()
#    for key,value in vardict.items():
#        if(key=="restart"): 
#            print key,value
if __name__ == "__main__":
     main()
