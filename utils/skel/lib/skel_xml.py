#!/usr/bin/env python

import argparse
import os
import xml.dom.minidom

import skel_bpy
import skel_settings


def pparse_command_line (parent_parser):
    parser = argparse.ArgumentParser (
                parents = [parent_parser],
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prog='skel',
                #add_help=False,
                description='''\
        skel xml
            create an xml file to define the I/O pattern for the target skeletal application''')

    parser.add_argument ('project', metavar='project', help='Name of the skel project')
    parser.add_argument ('-y', '--yaml-file', dest='yamlfile', help='yaml file to use for I/O pattern')
    parser.add_argument ('-f', '--force', dest='force', action='store_true', help='overwrite existing XML file')
    parser.set_defaults(force=False)

    return parser.parse_args()


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




def create_skel_xml (parent_parser):

    args = pparse_command_line (parent_parser)
    if args.yamlfile is not None:
        create_from_yaml (args.project, args)
    else:
        create_from_xml (args.project, args)


def create_from_yaml (project,args):
    #print "using yaml file"

    outfilename = project + '_skel.xml'

    # Only proceed if outfilename does not already exist, or if -f was used
    if os.path.exists (outfilename) and not args.force:
        print "%s exists, aborting. Delete the file or use -f to overwrite." % outfilename
        return 999

    skel_file = open (outfilename, 'w')

    bpy = skel_bpy.skel_bpy (args.yamlfile)


    # Okay, it's time to try out Cheetah.
    from Cheetah.Template import Template
    template_file = open (os.path.expanduser("~/.skel/templates/xml.tmpl"), 'r')
    t = Template(file=template_file)
    t.bpy = bpy
    skel_file.write (str(t) )
    # All done. That was easy.


def create_from_xml (project, args):

    outfilename = project + '_skel.xml'

    # Only proceed if outfilename does not already exist, or if -f was used
    if os.path.exists (outfilename) and not args.force:
        print "%s exists, aborting. Delete the file or use -f to overwrite." % outfilename
        return 999

    skel_file = open (outfilename, 'w')

    



    doc = xml.dom.minidom.parse (project + '.xml')
 
    # TODO: remove this application specific kludge.
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
        # TODO: Before adding this node, we should remove any text that is already here...
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


