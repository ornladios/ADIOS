#!/usr/bin/env python

import argparse
import os

import skelconf
import adios
import skel_settings


# To produce submit scripts, we'll work from a template. There will
# be two types of replacement, simple variables, and macros (for the
# tests)
def generate_submit_scripts (params):

    settings = skel_settings.skel_settings()

    for batch in params.get_batches():
        #platform = params.get_target()
        settings = skel_settings.skel_settings()
        platform = settings.get_submit_target()

        sfile = open ('submit_' + platform + '_' + batch.get_name(), 'w')
        sfile_template = open (os.path.expanduser('~/.skel/templates/submit_' + platform + '.tpl'), 'r')

        i = 0
        template_lines = sfile_template.readlines()
        while i < len (template_lines):

            template_line = template_lines[i]

            if '$$START_TEST$$' in template_line:
            # This is the test macro, run through it for each test
                template_start_index = i + 1
                for test in batch.get_tests():
                    j = template_start_index
                    template_line = template_lines[j]
                    while not '$$END_TEST$$' in template_line:
                        sfile.write (submit_line_template_replace (template_line, params, batch, test, settings))
                        j = j + 1
                        template_line = template_lines[j]
                    # Point at the first line after the macro
                    i = j + 1
            else:
                # Fill in any replacement vars in this line...
                template_line = submit_line_template_replace (template_line, params, batch, None, settings)
                sfile.write (template_line)
                i = i + 1

        sfile_template.close()
        sfile.close()


import re
import math

def submit_line_template_replace (template_line, params, batch, test, settings):

    template_line = template_line.replace ('$$JOB_NAME$$', batch.get_name() + '_%d'%batch.get_cores() + '_skel_' + params.get_application() )
    template_line = template_line.replace ('$$WALLTIME$$', batch.get_walltime() )
    template_line = template_line.replace ('$$APP$$', params.get_application() )
    template_line = template_line.replace ('$$CORES_USED$$', '%d'%batch.get_cores() )
    template_line = template_line.replace ('$$TARGET$$', params.get_target() )
    template_line = template_line.replace ('$$ACCOUNT$$', settings.get_account() )

    if test != None:

        #Test specific replacements
        template_line = template_line.replace ('$$TAGS$$', test.get_tags() )
        template_line = template_line.replace ('$$METHOD$$', test.get_method() )
        template_line = template_line.replace ('$$EXEC$$', params.get_application() + '_skel_' + test.get_group_name() + '_' + test.get_type() )

        template_line = template_line.replace ('$$ITERATIONS$$', test.get_iterations() )

        template_line = template_line.replace ('$$METHOD_PARAMS$$', test.get_method_params() )
        template_line = template_line.replace ('$$EXT$$', test.get_ext() )

        if test.get_rm() == 'pre' or test.get_rm() == 'both':
            prerm = 'rm -rf out*'
        else:
            prerm = ''
        template_line = template_line.replace ('$$PRE_RM$$', prerm)

        if test.get_rm() == 'post' or test.get_rm() == 'both':
            postrm = 'rm -rf out*'
        else:
            postrm = ''
        template_line = template_line.replace ('$$POST_RM$$', postrm)


    if '$$CORES_TOTAL$$' in template_line:
        pattern = re.compile (r"\$\$CORES_TOTAL\$\$[\d]*\$\$")
        match = pattern.search (template_line)
        match_term = match.group()

        # If we split the matched string at the dollar signs, the cores/node will be
        # at index 4
        count = float(match_term.split('$')[4])
        total_cores = int (math.ceil( (batch.get_cores() / count) ) * count)
        template_line = template_line.replace (match_term, '%d'%total_cores)


    if '$$NODES_TOTAL$$' in template_line:
        pattern = re.compile (r"\$\$NODES_TOTAL\$\$[\d]*\$\$")
        match = pattern.search (template_line)
        match_term = match.group()

        count = float(match_term.split('$')[4])
        total_nodes = int (math.ceil( (batch.get_cores() / count) ) )
        template_line = template_line.replace (match_term, '%d'%total_nodes)

    return template_line


def parse_command_line():

    parser = argparse.ArgumentParser (description='Create submission scripts for the given skel project')
    parser.add_argument ('project', metavar='project', help='Name of the skel project')

    return parser.parse_args()


def main(argv=None):

    skel_settings.create_settings_dir_if_needed()
    
    args = parse_command_line()

    config = adios.adiosConfig (args.project + '_skel.xml')
    params = skelconf.skelConfig (args.project + '_params.xml')

    #generate_makefiles_c (params)
    generate_submit_scripts (params)



if __name__ == "__main__":
    main()



