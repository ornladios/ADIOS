#!/usr/bin/env python
import sys
import os
import argparse

import adios
import skel_settings
import skel_bpls

# Command line parsing is chained together. This is stage two. The first stage happens in ../bin/skel
def pparse_command_line (parent_parser):

    parser = argparse.ArgumentParser (
                parents=[parent_parser],
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prog='skel',
                #add_help=False,
                description='''\
        skel params
            create a parameter file to define skeletal application behavior''')

    parser.add_argument ('project', metavar='project', help='Name of the skel project')
    parser.add_argument ('-g', '--group', help='adios group')
    parser.add_argument ('-b', '--bpls', help='file containing bpls output')

    parser.add_argument ('-f', '--force', dest='force', action='store_true', help='overwrite existing params file')
    parser.set_defaults(force=False)

    return parser.parse_args()


def generate_param_file_with_args (parent_parser):
    args = pparse_command_line (parent_parser)
   
    try:
        config = adios.adiosConfig (args.project + '_skel.xml')
    except (IOError):
        print "XXError reading " + args.project + "_skel.xml. Try running skel xml " + args.project + " first."
        return 1

 
    outfilename = args.project + '_params.xml'

    # Only proceed if outfilename does not already exist, or if -f was used
    if os.path.exists (outfilename) and not args.force:
        print "%s exists, aborting. Delete the file or use -f to overwrite." % outfilename
        return 999

    try:
        config = adios.adiosConfig (args.project + '_skel.xml')
    except (IOError):
        print "Error reading " + args.project + "_skel.xml. Try running skel xml " + args.project + " first."
        return 1

    generate_param_file (args.project, outfilename, config, args.group, args.bpls)


def generate_param_file (app, outfile, config, groupname, bplsfile=None):

    param_file = open (outfile, 'w')

    if bplsfile is not None:
        print "Using bpls data in %s" % bplsfile
        bpdata = skel_bpls.bpls (open (bplsfile, 'r') )

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
                if bplsfile is None:
                    all_scalars.add ('\n    <scalar name="' + var.get_name() + '" type="' + var.get_type() + '" value="128" />')
                else:
                    scalar_value = None

                    first_use_name, first_use_dim_num = var.find_first_use () # Get the name and dimension number of the first array that uses this scalar, or None if it is not used
                    if first_use_name is not None:
                        dims = bpdata.get_dims (first_use_name)
                        if dims is None:
                            # Try adding a leading slash to deal with the way that bpls reports variable names without one
                            dims = bpdata.get_dims ("/%s" % first_use_name)
                        if dims is not None:
                            scalar_value = dims[first_use_dim_num]
                    

                    if scalar_value is None:
                        scalar_value = 0 # Should be used only for variables that do not appear in any array dimensions

                    all_scalars.add ('\n    <scalar name="' + var.get_name() + '" type="' + var.get_type() + '" value="%s" />' % scalar_value)
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


# TODO: Get rid of this in favor of chained version, above.
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

    # Only proceed if outfilename does not already exist.
    if os.path.exists (outfilename):
        print "%s exists, aborting. Delete the file or use '-f' to overwrite."
        return 999

    generate_param_file (args.project, outfilename, config, args.group)
        

if __name__ == "__main__":
    main()

 
