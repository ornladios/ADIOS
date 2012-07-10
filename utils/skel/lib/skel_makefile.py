#!/usr/bin/env python

import argparse
import os
import sys

import adios
import skelconf
import skel_settings



def generate_makefiles (params, config):
    lang = config.get_host_language()
    if "C" == lang or "c" == lang:
        generate_makefiles_c (params)
    else:
        generate_makefiles_fortran (params)



def generate_makefiles_fortran (params):

    platform = params.get_target()
    print 'generating Makefile for target ' + platform

    settings = skel_settings.skel_settings()

    makefile = open ('Makefile', 'w')

    makefile_template_name = '~/.skel/templates/Makefile.' + platform + '.tpl'
    makefile_template = open(os.path.expanduser(makefile_template_name), 'r')

    include_statement = "" + os.path.dirname (sys.argv[0]) + '/../etc/skel/compiler_fragment.mk'

    for template_line in makefile_template:

        # Fill in any replacement vars in this line...
        template_line = template_line.replace ('$$APP$$', params.get_application () )
        template_line = template_line.replace ('$$INCLUDE$$', include_statement)
        template_line = template_line.replace ('$$TARGET$$', platform)
        template_line = template_line.replace ('$$DEPLOY_DIR$$', settings.get_deploy_dir())
        template_line = template_line.replace ('$$SKEL_HOME$$', settings.get_skel_home())
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
    print 'generating Makefile for target ' + platform

    settings = skel_settings.skel_settings()

    makefile = open ('Makefile', 'w')

    makefile_template_name = os.path.dirname (sys.argv[0]) + '/../etc/skel/templates/Makefile.' + platform + '.tpl'
    makefile_template = open(os.path.expanduser(makefile_template_name), 'r')

    include_statement = "" + os.path.dirname (sys.argv[0]) + '/../etc/skel/compiler_fragment.mk'

    for template_line in makefile_template:

        # Fill in any replacement vars in this line...
        template_line = template_line.replace ('$$APP$$', params.get_application () )
        template_line = template_line.replace ('$$INCLUDE$$', include_statement)
        template_line = template_line.replace ('$$TARGET$$', platform)
        template_line = template_line.replace ('$$DEPLOY_DIR$$', settings.get_deploy_dir())
        template_line = template_line.replace ('$$SKEL_HOME$$', settings.get_skel_home())
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



