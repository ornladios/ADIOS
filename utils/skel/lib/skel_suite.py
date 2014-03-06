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
import skel_test_plan


def pparse_command_line (parent_parser):
    parser = argparse.ArgumentParser (
                parents = [parent_parser],
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prog='skel',
                #add_help=False,
                description='''\
        skel suite 
            instantiate a suite of tests based on a given test plan''')

    parser.add_argument ('project', metavar='project', help='Name of the skel project')
    parser.add_argument ('-y', '--yaml-file', dest='yamlfile', help='yaml file to load I/O pattern')
    parser.add_argument ('-f', '--force', dest='force', action='store_true', help='overwrite existing source files')
    parser.set_defaults(force=False)

    return parser.parse_args()



def gen_suite_with_args (parent_parser):
    args = pparse_command_line (parent_parser)

    print "Generating test suite using %s" % args.yamlfile

    suite_gen_file_name = "%s_gen_suite.sh" % args.project
    suite_gen_file = open (suite_gen_file_name, "w")
    plan = skel_test_plan.skel_test_plan (args.yamlfile)

    # Generate gen_suite.sh shell script
    from Cheetah.Template import Template
    template_file = open (os.path.expanduser("~/.skel/templates/create_suite.tmpl"), 'r')
    t = Template(file=template_file)
    t.test_plan = plan

    # No, I don't like these either.
    t.yamlfile = args.yamlfile
    t.project = args.project
    t.force = args.force

    suite_gen_file.write (str(t) )

    suite_gen_file.close()

    # Adjust the permissions of the replay script to make it runnable by user
    os.chmod (suite_gen_file_name, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)

    # Run it
    #print "Run ./%s [disabled]" % suite_gen_file_name
    subprocess.check_call ("./%s" % suite_gen_file_name)




