#!/usr/bin/env python

import sys
import os

import ad_config
import type_mapper


# checkXML attempts to call adios_lint to verify that adios will accept the xml file as valid
def checkXML (config_file, path):
    if path == '':
        adios_lint = 'adios_lint '
    elif os.path.exists (path+'/adios_lint'):
        adios_lint = path + '/adios_lint '
    elif os.path.exists (path+'/../adios_lint/adios_lint'):
        adios_lint = path+'/../adios_lint/adios_lint '
    else:
        adios_lint = 'adios_lint '

    rv = os.system (adios_lint + config_file)
    if rv == 0:
        return 'success'
    elif rv == 32512:  # System unable to find adios_lint command
        print "Unable to find adios_lint. Proceeding with code generation."
        return 'success'
    else:
        print "gpp.py failed."
        return 'failure'


def get_c_groupsize_code (group):
    groupsize_code_string = ''
    groupsize_code_string += 'adios_groupsize = '
    for v in group.get_vars():
        if v.is_scalar():
            if v.get_type() != 'string':
                groupsize_code_string += ('%d' % type_mapper.get_size (v.get_type() ) + ' \\\n                + ')
            else:
                groupsize_code_string += ('strlen(' + v.get_gwrite() + ')' + ' \\\n                + ')
        else:
            groupsize_code_string += ('%d * ' % type_mapper.get_size (v.get_type() ) )

            for d in v.get_dimensions():
                # need to check whether this is the timestep
                groupsize_code_string += '(' + d + ') * '

            groupsize_code_string = groupsize_code_string.rstrip ('* ')

            groupsize_code_string += (' \\\n                + ')

    # remove the final +\, and add the ;
    groupsize_code_string = groupsize_code_string.rstrip('+\\\n ') + ';'

    groupsize_code_string += '\nadios_group_size (adios_handle, adios_groupsize, &adios_totalsize);'

    return groupsize_code_string;


def get_fortran_groupsize_code (group):
    #print 'Get Fortran Groupsize code for group "'+group.get_name()+'"'
    groupsize_code_string = ''
    groupsize_code_string += 'adios_groupsize = '
    for v in group.get_vars():
        #print '  count variable "'+v.get_fullpath()+'"'
        if (v.is_scalar() ):
            if v.get_type() != 'string':
                groupsize_code_string += ('%d' % type_mapper.get_size (v.get_type() ) + ' &\n                + ')
            else:
                groupsize_code_string += ('len_trim(' + v.get_gwrite() + ')' + ' &\n                + ')
        else:
            groupsize_code_string += ('%d * ' % type_mapper.get_size (v.get_type() ) )

            for d in v.get_dimensions():
                #print '  count dim "'+d+'"'
                # need to check whether this is the timestep
                groupsize_code_string += '(' + d + ') * '

            groupsize_code_string = groupsize_code_string.rstrip ('* ')

            groupsize_code_string += (' &\n                + ')

    # remove the final &+
    groupsize_code_string = groupsize_code_string.rstrip('+\n &') 

    groupsize_code_string += '\ncall adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)'

    #print 'Done Fortran Groupsize'
    return groupsize_code_string;


def get_fortran_write_statements (group):

    statements = ""

    for item in group.get_ordered_contents():

        if item.__class__.__name__ == 'gwrite':
            statements += '\n' + item.get_src() 
            continue

        if item.__class__.__name__ == 'attr':
            # Ignore the attributes for now
            continue

        # Otherwise, the item must be a variable
        var = item

        statements += '\ncall adios_write (adios_handle, "' + var.get_fullpath() + '", ' + var.get_gwrite() + ', adios_err)'

    statements += '\n'
    return statements



def get_c_write_statements (group):

    statements = ""

    for item in group.get_ordered_contents():

        if item.__class__.__name__ == 'gwrite':
            statements += '\n' + item.get_src() 
            continue

        if item.__class__.__name__ == 'attr':
            # Ignore the attributes for now
            continue

        # Otherwise, the item must be a variable
        var = item

        # The tricky bit here is deciding whether we need the & before the variable name.
        # We omit it in two cases: 1) the variable type is string, or 2) the variable is not a scalar
        if (var.get_c_type() == 'string' or var.get_dimensions() != None):
            var_prefix = ''
        else:
            var_prefix = '&'

        statements += '\nadios_write (adios_handle, "' + var.get_fullpath() + '", ' + var_prefix + var.get_gwrite() + ');'

    statements += '\n'
    return statements


def get_fortran_read_statements (group):

    statements = ""
    # Make a selection to capture writes done by the corresponding process in the writing application
    statements += 'call adios_selection_writeblock (s, rank)\n'

    for var in group.get_vars():
        if var.get_dimensions() == None:
            continue

        statements += 'call adios_schedule_read (fp, s, "' + var.get_name() + '", 0, 1, ' + var.get_gwrite() + ', adios_err)\n'

    statements += '\ncall adios_perform_reads (fp, adios_err)\n'
    statements += 'call adios_selection_delete (s)\n'
    return statements


def get_c_read_statements (group):

    statements = ""
    # Make a selection to capture writes done by the corresponding process in the writing application
    statements += 's = adios_selection_writeblock (rank);\n'
    #statements += 'uint64_t counts[64];\n'

    for var in group.get_vars():
        if var.get_dimensions() == None:
            continue

        # The tricky bit here is deciding whether we need the & before the variable name.
        # We omit it in two cases: 1) the variable type is string, or 2) the variable is not a scalar
        if (var.get_c_type() == 'string' or var.get_dimensions() != None):
            var_prefix = ''
        else:
            var_prefix = '&'

        statements += 'adios_schedule_read (fp, s, "' + var.get_name() + '", 0, 1, ' + var_prefix + var.get_gwrite() + ');\n'

    statements += 'adios_perform_reads (fp, 1);\n'
    statements += 'adios_selection_delete (s);\n'
    return statements



def generate_fortran (config):

    # Output a gwrite_*.fh file for each group in the config
    for group in config.get_groups():
        # Open file - the filename is gwrite_<group>.fh
        gwfile = open ('gwrite_' + group.get_name() + '.fh', 'w')

        # Set groupsize
        gwfile.write (get_fortran_groupsize_code (group) )

        # Write vars
        gwfile.write (get_fortran_write_statements (group) )

        # Close file
        gwfile.close()

    # Output a gread_*.fh file for each group in the config
    for group in config.get_groups():
        # Open file - the filename is gread_<group>.fh
        grfile = open ('gread_' + group.get_name() + '.fh', 'w')

        # Read vars
        grfile.write (get_fortran_read_statements (group) )

        # Close file
        grfile.close()


def generate_c (config):
    #print "c"

    # Output a gwrite_*.ch file for each group in the config
    for group in config.get_groups():
        # Open file - the filename is gwrite_<group>.ch
        gwfile = open ('gwrite_' + group.get_name() + '.ch', 'w')

        # Set groupsize
        gwfile.write (get_c_groupsize_code (group) )

        # Write vars
        gwfile.write (get_c_write_statements (group) )

        # Close ch file
        gwfile.close()

    # Output a gread_*.ch file for each group in the config
    for group in config.get_groups():
        # Open file - the filename is gread_<group>.ch
        grfile = open ('gread_' + group.get_name() + '.ch', 'w')

        # Read vars
        grfile.write (get_c_read_statements (group) )

        # Close ch file
        grfile.close()


def main (argv=None):
  
    # Must be called with one argument, the name of the xml file
    if len (sys.argv) != 2:
        print 'Usage: gpp.py <config file>\n'
        return 1


    # If possible, run adios_lint to catch malformed config files
    check_val = checkXML (sys.argv[1], os.path.dirname(sys.argv[0]))
    if check_val != 'success':
        return 1


    # The ad_config object provides all relevant info from the XML file
    config = ad_config.adiosConfig (sys.argv[1])


    # Now just call the generator for the specified language
    lang = config.get_host_language().lower()
    if lang == 'fortran':
        generate_fortran (config)
    elif lang == 'c' or lang == 'cpp':
        generate_c (config)
    else:
        print "Fatal: Language unknown or unspecified"
        raise SystemExit



if __name__ == "__main__":
    main()

