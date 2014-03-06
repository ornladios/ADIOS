#!/usr/bin/env python

import argparse
import os
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
                add_help=True,
                description='''\
        skel makefile 
            create a makefile for building a skeletal application''')

    parser.add_argument ('project', metavar='project', help='Name of the skel project')
    parser.add_argument ('-y', '--yaml-file', dest='yamlfile', help='yaml file to store I/O pattern')
    parser.add_argument ('-b', '--bp-file', dest='bpfile', help='bp file to extract I/O pattern')
    parser.add_argument ('-f', '--force', dest='force', action='store_true', help='overwrite existing source file')
    parser.add_argument ('-n', '--noxml', dest='noxml', action='store_true', help='generate noxml code')
    parser.set_defaults(force=False)
    parser.set_defaults(noxml=False)

    return parser.parse_args()



def generate_makefiles_with_args (parent_parser):
    args = pparse_command_line (parent_parser)

    try:
        config = adios.adiosConfig (args.project + '_skel.xml')
    except (IOError):
        print "XXError reading " + args.project + "_skel.xml. Try running skel xml " + args.project + " first."
        return 1


    if args.yamlfile:
        generate_makefile_from_yaml (args)
    else:

        try:
            params = skelconf.skelConfig (args.project + '_params.xml')
        except (IOError):
            print "Error reading " + args.project + "_params.xml. Try running skel params " + args.project + " first,"
            print "then check that " + args.project + "_params.xml exists."
            return

        generate_makefiles (params, config)


def generate_makefile_from_yaml (args):

    bpy = skel_bpy.skel_bpy (args.yamlfile)

    template_file_name = "~/.skel/templates/Makefile.tmpl"
    outfilename = "Makefile"

    # Only proceed if outfilename does not already exist, or if -f was used
    if os.path.exists (outfilename) and not args.force:
        print "%s exists, aborting. Delete the file or use -f to overwrite." % outfilename
        return 999

    skel_file = open (outfilename, 'w')


    # Now for the Cheetah magic:
    from Cheetah.Template import Template
    template_file = open (os.path.expanduser(template_file_name), 'r')
    t = Template(file=template_file)

    t.bpy = bpy
    t.project = args.project
    t.bpfile = args.bpfile
    t.noxml = args.noxml
    skel_file.write (str(t) )
 


def generate_makefiles (params, config):
    lang = config.get_host_language()
    if "C" == lang or "c" == lang:
        generate_makefiles_c (params)
    else:
        generate_makefiles_fortran (params)



def generate_makefiles_fortran (params):

    platform = params.get_target()

    settings = skel_settings.skel_settings()

    makefile = open ('Makefile', 'w')

    # Makefile generation no longer depends on the target, just using default here.
    #makefile_template_name = '~/.skel/templates/Makefile.' + platform + '.tpl'
    makefile_template_name = '~/.skel/templates/Makefile.default.tpl'
    makefile_template = open(os.path.expanduser(makefile_template_name), 'r')

    include_statement = "" + os.path.dirname (sys.argv[0]) + '/../etc/skel/compiler_fragment.mk'

    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))        

    for template_line in makefile_template:

        # Fill in any replacement vars in this line...
        template_line = template_line.replace ('$$ADIOS_BIN_DIR$$', bindir)
        template_line = template_line.replace ('$$APP$$', params.get_application () )
        template_line = template_line.replace ('$$INCLUDE$$', include_statement)
        template_line = template_line.replace ('$$TARGET$$', platform)
        template_line = template_line.replace ('$$DEPLOY_DIR$$', settings.get_deploy_dir())
        template_line = template_line.replace ('$$CORES_USED$$', '%d'%params.get_batches()[0].get_cores())

        if '$$CTESTS$$' in template_line:
            template_line = template_line.replace ('$$CTESTS$$', '')

        if '$$FTESTS$$' in template_line:
            test_string = ''
            test_set = set()
            for batch in params.get_batches():
                for test in batch.get_tests():
                    test_set.add (params.get_application() + '_skel_' + test.get_group_name() + '_' + test.get_type() + ' ')

            for t in test_set:
                test_string = test_string + t
            template_line = template_line.replace ('$$FTESTS$$', test_string)

        makefile.write (template_line)


def generate_makefiles_c (params):

    platform = params.get_target()

    settings = skel_settings.skel_settings()

    makefile = open ('Makefile', 'w')

    # Makefile generation no longer depends on the target, just using default here.
    #makefile_template_name = os.path.dirname (sys.argv[0]) + '/../etc/skel/templates/Makefile.' + platform + '.tpl'
    makefile_template_name = os.path.dirname (sys.argv[0]) + '/../etc/skel/templates/Makefile.default.tpl'
    makefile_template = open(os.path.expanduser(makefile_template_name), 'r')

    include_statement = "" + os.path.dirname (sys.argv[0]) + '/../etc/skel/compiler_fragment.mk'
    
    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))        

    for template_line in makefile_template:

        # Fill in any replacement vars in this line...
        template_line = template_line.replace ('$$ADIOS_BIN_DIR$$', bindir)
        template_line = template_line.replace ('$$APP$$', params.get_application () )
        template_line = template_line.replace ('$$INCLUDE$$', include_statement)
        template_line = template_line.replace ('$$TARGET$$', platform)
        template_line = template_line.replace ('$$DEPLOY_DIR$$', settings.get_deploy_dir())

        template_line = template_line.replace ('$$CORES_USED$$', '%d'%params.get_batches()[0].get_cores())

        if '$$FTESTS$$' in template_line:
            template_line = template_line.replace ('$$FTESTS$$', '')

        if '$$CTESTS$$' in template_line:
            test_string = ''
            test_set = set()
            for batch in params.get_batches():
                for test in batch.get_tests():
                    test_set.add (params.get_application() + '_skel_' + test.get_group_name() + '_' + test.get_type() + ' ')

            for t in test_set:
                test_string = test_string + t
            template_line = template_line.replace ('$$CTESTS$$', test_string)

        makefile.write (template_line)


def parse_command_line():

    parser = argparse.ArgumentParser (description='Create a Makefile for the given skel project')
    parser.add_argument ('project', metavar='project', help='Name of the skel project')

    return parser.parse_args()


def main(argv=None):

    skel_settings.create_settings_dir_if_needed()

    args = parse_command_line()

    config = adios.adiosConfig (args.project + '_skel.xml')
    params = skelconf.skelConfig (args.project + '_params.xml')

    lang = config.get_host_language ()
    if 'c' == lang or 'C' == lang:
        print 'generating C flavored Makefile'
        generate_makefiles_c (params)
    else:
        print 'generating fortran flavored Makefile'
        generate_makefiles_fortran (params)



if __name__ == "__main__":
    main()



