#!/usr/bin/env python
import sys
import os
import argparse

import adios
import skel_settings


def generate_param_file (app, outfile, config, groupname):

    param_file = open (outfile, 'w')

    #Write the file header
    param_file.write ('<?xml version="1.0"?>')
    param_file.write ('\n<skel-config application="' + app + '">')

    param_file.write ('\n\n<!--')
    param_file.write ('\n  Within each group, use the scalar elements to control things like array sizes and offsets.')
    param_file.write ('\n  Simply adjust the value attribute as needed. The type is provided for convenience.')
    param_file.write ('\n  Note that there are 2 special values that you can use:')
    param_file.write ('\n    skel_mpi_size refers to the number of processes participating in this run, and')
    param_file.write ('\n    skel_mpi_rank is used to indicate the rank of the local process')
    param_file.write ('\n  -->\n')


    #Write a section for each group of interest
    for group in config.get_groups():

        # if we've specified a particular group, ignore all of the other groups
        if (groupname != None and groupname != group.get_name() ): 
            continue

        param_file.write ('\n\n  <adios-group name="' + group.get_name() + '">')

        all_scalars = set()
        all_arrays = set()

        for var in group.get_vars():
            if var.is_scalar():
                all_scalars.add ('\n    <scalar name="' + var.get_name() + '" type="' + var.get_type() + '" value="128" />')
            else:
                dims = var.get_dimensions()
                dim_str ='dims="'
                for dim in dims:
                    dim_str = dim_str + dim + ','
                dim_str = dim_str.rstrip(',')
                dim_str = dim_str + '"'

                all_arrays.add ('\n    <array name="' + var.get_gwrite() + '" type="' + var.get_type() + '" ' + dim_str + ' fill-method="rank"></array>')

        for s in all_scalars:
            param_file.write (s)
        for a in all_arrays:
            param_file.write (a)

        param_file.write ('\n  </adios-group>')

    # Make a test run for all of the writes
    param_file.write ('\n\n  <batch name="writes" cores="128" walltime="0:30:00">')
    for group in config.get_groups():

        param_file.write ('\n    <test type="write" group="' + group.get_name() + '" method="POSIX" iterations="10" rm="pre" tags="name1:val1,name2:val2" />')

    param_file.write ('\n  </batch>')


    #Write the footer

    param_file.write ('\n\n</skel-config>')
    param_file.close()



def parse_command_line():

    parser = argparse.ArgumentParser (description='Create a parameter file for the given skel project')
    parser.add_argument ('project', metavar='project', help='Name of the skel project')
    parser.add_argument ('-g', '--group', help='If specified, produce output only for this group')    
    return parser.parse_args()



def main(argv=None):

    skel_settings.create_settings_dir_if_needed()

    args = parse_command_line()

    config = adios.adiosConfig (args.project + '_skel.xml')

    # Determine outfile name
    outfilename = args.project + '_params.xml.default'

    generate_param_file (args.project, outfilename, config, args.group)
        

if __name__ == "__main__":
    main()

 
