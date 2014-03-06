#!/usr/bin/env python

import argparse
import os
import stat
import subprocess
import sys

import adios
import skelconf
import skel_bpy
import skel_settings


def pparse_command_line (parent_parser):
    parser = argparse.ArgumentParser (
                parents = [parent_parser],
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prog='skel',
                #add_help=False,
                description='''\
        skel replay 
            Generate an entire skeletal application that duplicates a specified I/O''')

    parser.add_argument ('project', metavar='project', help='Name of the skel project')
    parser.add_argument ('-y', '--yaml-file', dest='yamlfile', help='yaml file to load I/O pattern')
    parser.add_argument ('-b', '--bp-file', dest='bpfile', help='bp file to extract I/O pattern')
    parser.add_argument ('-f', '--force', dest='force', action='store_true', help='overwrite existing source files')
    parser.add_argument ('-n', '--noxml', dest='noxml', action='store_true', help='generate noxml code')
    parser.set_defaults(force=False)
    parser.set_defaults(noxml=False)

    return parser.parse_args()



def do_replay_with_args (parent_parser):
    args = pparse_command_line (parent_parser)

    if args.bpfile:
        do_replay_from_bpfile (args)
        return

    if args.yamlfile:
        do_replay_from_yaml (args)
        return

    print "No bp file or yaml file specified, exiting"
    return

#    else:
#        try:
#            params = skelconf.skelConfig (args.project + '_params.xml')
#        except (IOError):
#            print "Error reading " + args.project + "_params.xml. Try running skel params " + args.project + " first,"
#            print "then check that " + args.project + "_params.xml exists."
#            return

#        generate_makefiles (params, config)




def do_replay_from_bpfile (args):
    print "Replaying using %s" % args.bpfile
    
    # First, call skeldump to get the yamlfile
    sdcmd = "skeldump %s > %s.yaml" % (args.bpfile, args.project)
    print (sdcmd)
    os.system (sdcmd)

    # Now just do the replay from the yamlfile
    args.yamlfile = "%s.yaml" % args.project
    do_replay_from_yaml (args)



def do_replay_from_yaml (args):
    print "Replaying using %s" % args.yamlfile

    replay_file_name = "%s_replay.sh" % args.project
    replay_file = open (replay_file_name, "w")
    bpy = skel_bpy.skel_bpy (args.yamlfile)

    # Generate replay_yaml.sh shell script
    from Cheetah.Template import Template
    template_file = open (os.path.expanduser("~/.skel/templates/replay_yaml.tmpl"), 'r')
    t = Template(file=template_file)
    t.bpy = bpy
    t.noxml = args.noxml

    # No, I don't like these either.
    t.yamlfile = args.yamlfile
    t.project = args.project
    t.force = args.force

    replay_file.write (str(t) )

    replay_file.close()

    # Adjust the permissions of the replay script to make it runnable by user
    os.chmod (replay_file_name, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)

    # Run it
    subprocess.check_call ("./%s" % replay_file_name)




