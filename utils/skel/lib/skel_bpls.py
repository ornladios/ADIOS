#!/usr/bin/env python

import re

class bpls:
    def __init__ (self, file):

        self.vars = {} # This will be a dictionary of dictionaries. The first level is 
                       # Indexed by the variable names, and the second level has all of
                       # the specific information for each variable (initially just type
                       # and dimensions)

        # For now, assume that the input is from a bpls call with no cl arguments given
        # Can add more flexibility later by checking the output first and then parsing
        # accordingly
        for line in file:
            var_dict = {}

            tokens = line.split()

            # The first item is the type
            var_dict ['type'] = tokens[0]
            
            # Now parse the last item, which is either 'scalar' or a comma separated list of 
            # integer dimensions wrapped in curly braces, i.e. {7, 8, 9}
            if tokens[-1] == 'scalar':
                var_dict ['dims'] = None
            else:
                start = line.rindex('{') + 1
                end = line.rindex('}')
                var_dict ['dims'] = line[start:end].split (', ')

            # Now everything that is left, minus external whitespace, is the variable name.
            # There is a small hole here which is that if the variable name ends with
            # a space, there is no way to tell, since bpls fills with extra spaces to 
            # align the columns. It is probably a bad idea to have variable names that
            # end with spaces anyway, so I'm not going to worry about this right now.

            line = line.strip()
            start = line.index (' ') + 1 
            if tokens[-1] == 'scalar':
                end = -6
            else:
                end = line.rindex ('{')

            var_dict ['name'] = line[start:end].strip()


            # Put this var in the top level map according to its name
            self.vars [var_dict['name'] ] = var_dict


    def get_vars (self):
        return self.vars.keys()


    def get_dims (self, var):
        print "getting dims for %s" % var
        if var not in self.vars.keys():
            return None
        return self.vars[var]['dims']




def main(argv=None):


#    args = parse_command_line()

    test = open ("gts.bpls")

    b = bpls (test)

    for var in b.get_vars():
        print '%s    %s' % (var, b.get_dims (var) ) 
        

if __name__ == "__main__":
    main()



