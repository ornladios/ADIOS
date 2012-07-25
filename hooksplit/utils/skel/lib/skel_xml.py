#!/usr/bin/env python

import argparse
import xml.dom.minidom

import skel_settings


def parse_command_line():
    parser = argparse.ArgumentParser (description='Rewrite elements of the adios config file that would produce invalid syntax in the generated code, e.g. gwrite attributes containing + or %')
    parser.add_argument ('project', metavar='project', help='Name of the skel project')

    return parser.parse_args()


def cleanse (str):

    # Trim the name at the first open paren to deal with gts issue
    str = str.split('(')[0]
    
    # replace unfriendly characters with safer versions
    str = str.replace ("+", "__plus__")
    str = str.replace ("%", "__pct__")
    str = str.replace ("/", "__div__")
    str = str.replace ("-", "__minus__")
    str = str.replace ("**", "__pow__")
    str = str.replace ("*", "__mul__")


    return str




def create_skel_xml (project):

    outfilename = project + '_skel.xml'

    skel_file = open (outfilename, 'w')

    doc = xml.dom.minidom.parse (project + '.xml')
  
    #remove diagnostic groups (xgc1)

    groups = doc.getElementsByTagName ('adios-group')
    for g in groups:
        if g.getAttribute('name').startswith ('diag'):
            # remove this group
            doc.getElementsByTagName ('adios-config')[0].removeChild(g)

    methods = doc.getElementsByTagName ('method')
    for m in methods:
        if m.getAttribute ('group').startswith ('diag'):
            # remove this method
            doc.getElementsByTagName ('adios-config')[0].removeChild(m)

    for m in methods:
        token = doc.createTextNode ('***skel-parameters***')
        m.appendChild (token)



	

    #deal with the language
    lang = doc.getElementsByTagName ('adios-config')[0].getAttribute ('host-language')

#    #for now, assume the skeletal will be C
#    if (lang == 'Fortran'):
#        flip_dims = 'True'
#        doc.getElementsByTagName ('adios-config')[0].setAttribute ('host-language', 'C')
#    else:
#        flip_dims = 'False'

    #Do not change the language
    flip_dims = 'False'

    
    #clean up the variables
    vars = doc.getElementsByTagName ('var')
    for v in vars:
        gwrite = v.getAttributeNode ('gwrite')
        dims = v.getAttributeNode ('dimensions')
        name = v.getAttributeNode ('name')
        if gwrite != None:
            gwrite.value = cleanse (gwrite.value)
        if dims != None:
            dims.value = cleanse (dims.value)
            if flip_dims == 'True':
                dims.value = ','.join(dims.value.split(',')[::-1] )
        if name != None:
            name.value = cleanse (name.value)


    global_bounds = doc.getElementsByTagName ('global-bounds')
    for gb in global_bounds:
        dims = gb.getAttribute ('dimensions')
        offsets = gb.getAttribute ('offsets')

        if dims != None:
            dims = cleanse (dims)
            if flip_dims == 'True':
                dims = ','.join (dims.split(',')[::-1] )
            gb.setAttribute ('dimensions', dims)

        if offsets != None:
            offsets = cleanse (offsets)
            if flip_dims == 'True':
                offsets = ','.join (offsets.split(',')[::-1] )
            gb.setAttribute ('offsets', offsets)

    xmlstr = doc.toxml() 

    # ADIOS stumbles when there is no newline after the xml declaration, so we'll add one here
    xmlstr = xmlstr.replace ("?>", "?>\n", 1)

    skel_file.write (xmlstr)



def main ():

    skel_settings.create_settings_dir_if_needed()

    args = parse_command_line()

    create_skel_xml (args.project)


if __name__ == "__main__":
    main()


