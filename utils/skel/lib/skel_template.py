#!/usr/bin/env python

import argparse
import os
import stat
import subprocess
import sys
import yaml

import adios
import skel_settings


def pparse_command_line (parent_parser):
    parser = argparse.ArgumentParser (
                parents = [parent_parser],
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prog='skel',
                #add_help=False,
                description='''\
        skel template 
            Apply the list of items in the given yaml file to a template''')

    parser.add_argument ('-y', '--yaml-file', dest='yamlfile', required=True, help='yaml file to provide template parameters')
    parser.add_argument ('-o', '--output-file', dest='outfile', required=True, help='output file')
    parser.add_argument ('-t', '--template-file', dest='templatefile', required=True, help='template file')
    parser.add_argument ('-f', '--force', dest='force', action='store_true', help='overwrite existing source files')
    parser.set_defaults(force=False)

    return parser.parse_args()



def fill_template (parent_parser):
    args = pparse_command_line (parent_parser)

    #print "filling template"


    output_file = open (args.outfile, "w")
    
    #Get the yaml items visible for cheetah...
    stream = file (args.yamlfile, 'r')
    items = yaml.load(stream)

    #Fill the template 
    from Cheetah.Template import Template
    template_file = open (args.templatefile, 'r')
    t = Template(file=template_file)

    # No, I don't like these either.
    t.i = items 
    #print items


    output_file.write (str(t) )

    output_file.close()





