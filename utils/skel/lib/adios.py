import xml.dom.minidom

import typeMapper

class adiosConfig:

    def __init__ (self, config_file_name):

        self.config_file_name = config_file_name

        #This would be a good time to parse the file...
        doc = xml.dom.minidom.parse (config_file_name)
        nodes = doc.childNodes
        if (nodes.length != 1):
            print 'malformed adios config file, should contain only a single adios-config element'
            raise SystemExit
        self.config_node = nodes[0]

        # for each of the groups, instantiate an adiosGroup object, and store in self.adios_groups
        self.adios_groups = []
        self.methods = []
        self.buffer = None

        self.host_language = self.config_node.getAttribute ('host-language')

        for node in self.config_node.getElementsByTagName('adios-group'):
            self.adios_groups.append (adiosGroup (node) )

        for node in self.config_node.getElementsByTagName('method'):
            self.methods.append (method (node) )

        for node in self.config_node.getElementsByTagName ('buffer'):
            # there should be only one of these... this code ignores all but the last one.
            self.buffer = buffer (node)

        #We are currently ignoring any analysis declarations


    def get_filename (self):
        return self.config_file_name


    def get_groups (self):
        #return the group with the specified name
        return self.adios_groups


    def get_buffer (self):
        #return the buffer info
        print 'implement get_buffer'

    def get_host_language (self):
        return self.host_language

class adiosGroup:

    def __init__ (self, group_node):
        self.group_node = group_node

        self.time_index = self.group_node.getAttribute ('time-index')

        self.vars = []
        self.vardict = {}

        for node in self.group_node.childNodes:
            if (node.localName == 'var'):
                newvar = var (node, self, self.time_index)
                self.vars.append (newvar)
                self.vardict [newvar.get_name()] = newvar

        for node in self.group_node.getElementsByTagName('global-bounds'):
            for varnode in node.getElementsByTagName('var'):
                newvar = var (varnode, self, self.time_index, node) 
                self.vars.append (newvar)
                self.vardict [newvar.get_name()] = newvar

    def get_name (self):
        return self.group_node.getAttribute ('name')

    def get_vars (self):
        return self.vars

    def get_var (self, varname):
        return self.vardict [varname]



class method:
    
    def __init__ (self, method_node):
        self.method_node = method_node


class buffer:
    
    def __init__ (self, buffer_node):
        self.buffer_node = buffer_node


class var:

    def __init__ (self, var_node, group, time_index=None, global_bounds_node=None):
        self.var_node = var_node
        self.group = group
        self.time_index = time_index
        self.global_bounds_node = global_bounds_node

    def get_name (self):
        name = self.var_node.getAttribute ('name')
        name = name.replace ("+", "_plus_")
        name = name.replace ("%", "_pct_")
        name = name.split ('(')[0]
        return name

    def get_path (self):
        path = self.var_node.getAttribute ('path')
        return path

    def get_fullpath (self):
        path = self.get_path()
        name = self.get_name()
        if (path == ''):
            fullpath = name
        elif (path == '/'):
            fullpath = '/'+name
        else:
            fullpath = path + '/' + name
        return fullpath

    def get_gwrite (self):
        gw = self.var_node.getAttribute ('gwrite')

        if (gw == ''):
            gw = self.get_name()
    
        # Trim the name at the first open paren to deal with gts issue
        gw = gw.split('(')[0]
        gw = gw.replace ("+", "_plus_")
        gw = gw.replace ("%", "_pct_")

        return gw

    def get_group (self):
        return self.group

    def get_c_type (self):
        return typeMapper.get_c_type (self.var_node.getAttribute ('type') )

    def get_fortran_type (self):
        return typeMapper.get_fortran_type (self.var_node.getAttribute ('type') )

    def get_type (self):
        return self.var_node.getAttribute ('type')

    def get_dimensions (self):
        if (self.var_node.getAttribute ('dimensions') == ''):
            return None
        else:
            # place the dimensions in a list and remove the time-index if it is there.
            dims = filter (lambda x : x != self.time_index, self.var_node.getAttribute ('dimensions').split(',') )
            cleandims = []
            for d in dims:

                if d.isdigit():
                    cleandims.append (d)
                    continue

                # Here we need to translate the variable name for this dimension (if it's a var) into the gwrite 
                # for that variable
                dim_var = self.get_group().get_var (d)				
                if dim_var != None:
                    d = dim_var.get_gwrite()

				# Now clean up any unsightly parts
                cleand = d.replace ("+", "_plus_") 
                cleand = cleand.split('(')[0]
                cleandims.append (cleand)
            return cleandims

    def is_scalar (self):
        return self.get_dimensions() == None

    # TODO: Implement this
    def find_first_use (self):
        # Loop through all of the vars in the group
        for var in self.group.get_vars():
            dim_num = 0;
            if var.get_dimensions() is not None:
                for dim in var.get_dimensions():
                    # if this one uses this variable as a dimension, return the name and dim number
                    if dim == self.get_name():
                        return var.get_name(), dim_num
                    dim_num = dim_num + 1

        # None found, return None,None
        return None,None


class fortranFormatter:
    @staticmethod
    def get_write_line (var):
        retval = '\n  call adios_write (adios_handle, "' + var.get_fullpath() + '", ' + var.get_gwrite() + ', adios_error)'  
        #print retval
        return retval

    @staticmethod
    def get_declaration (var, group_params):
        dims = var.get_dimensions()
        if (dims == None):
            # I think this should be get_name instead of get_gwrite. 
            #init_val = group_params.get_scalar (var.get_gwrite() ).get_value()
            init_val = group_params.get_scalar (var.get_name() ).get_value()
            return '\n  ' + var.get_fortran_type() + ' :: ' + var.get_gwrite() 
        else:
            #fill_method = group_params.get_array (var.get_gwrite() ).get_fill_method()
            dimspec = ''
            for d in dims:
                dimspec = dimspec + ':,'
            dimspec = dimspec.rstrip (',')
            return '\n  ' + var.get_fortran_type() + ', ALLOCATABLE, DIMENSION(' + dimspec + ') :: ' + var.get_gwrite()  

    @staticmethod
    def get_initialization (var, group_params):
        dims = var.get_dimensions()
        if (dims == None):
            # I think this should be get_name instead of get_gwrite. 
            #init_val = group_params.get_scalar (var.get_gwrite() ).get_value()
            init_val = group_params.get_scalar (var.get_name() ).get_value()
            return '\n  ' + var.get_gwrite() + ' = ' + init_val 
        else:
            fill_method = group_params.get_array (var.get_gwrite() ).get_fill_method()

            return '\n  allocate (' + var.get_gwrite() + '(' + fortranFormatter.get_dim_str (var, ',') + ') )'

            #return '\n  ' + var.get_gwrite() + ' = (' + var.get_c_type() + '*) malloc (' + cFormatter.get_dim_str (var, '*') + ' * sizeof (' + var.get_c_type() + ') );\n' + cFormatter.get_fill (var, fill_method)

    @staticmethod
    def get_dim_str (var, sep):
        rv = ''
        for d in var.get_dimensions():
            rv += d
            rv += sep
        return rv.rstrip (sep)

    @staticmethod
    def get_groupsize_code (group):
        groupsize_code_string = ''
        groupsize_code_string += '\n\n! Set the adios group size'
        groupsize_code_string += '\n  adios_groupsize = &'
        for v in group.get_vars():
            if (v.is_scalar() ):
                groupsize_code_string += ('\n                     %d +' % typeMapper.get_size (v.get_type() ) + ' &')
            else:
                groupsize_code_string += ('\n                     %d * ' % typeMapper.get_size (v.get_type() ) )
                
                for d in v.get_dimensions():
                    # need to check whether this is the timestep
                    groupsize_code_string += '(' + d + ') * '

                groupsize_code_string = groupsize_code_string.rstrip ('* ')

                groupsize_code_string += (' + &')

        # remove the final +, and add the ;
        groupsize_code_string = groupsize_code_string.rstrip('+ &')

        groupsize_code_string += '\n  call adios_group_size (adios_handle, adios_groupsize, skel_total_size, adios_error)'

        return groupsize_code_string;



class cFormatter:
    @staticmethod
    def get_write_line (var):
        # The tricky bit here is deciding whether we need the & before the variable name.
        # We omit it in two cases: 1) the variable type is string, or 2) the variable is not a scalar
        if (var.get_c_type() == 'string' or var.get_dimensions() != None):
            var_prefix = ''
        else:
            var_prefix = '&'

        retval = '\nadios_write (adios_handle, "' + var.get_fullpath() + '", ' + var_prefix + var.get_gwrite() + ');'  
        #print retval
        return retval

    @staticmethod
    def get_read_all_line (var):
        if (var.get_c_type() == 'string' or var.get_dimensions() != None):
            var_prefix = ''
        else:
            var_prefix = '&'

        return '\nadios_write (adios_handle, "' + var.get_name() + '", ' + var_prefix + var.get_gwrite() + ');'

    @staticmethod
    def get_dim_str (var, sep):
        rv = ''
        for d in var.get_dimensions():
            rv += d
            rv += sep
        return rv.rstrip (sep)


    @staticmethod
    def get_fill (var, method):
        fill_content = ''

        if (method == 'rank'):
            dims = var.get_dimensions()
            fill_content = ''

            fill_content += 'for (skel_i = 0; skel_i < ' + cFormatter.get_dim_str (var, '*') + '; skel_i++) \n'

            fill_content += '    ' + var.get_gwrite() + '[skel_i] = (' + var.get_c_type() + ') skel_mpi_rank;'

        return fill_content


    @staticmethod
    def get_declaration (var, group_params):
        dims = var.get_dimensions()
        if (dims == None):
            # I think this should be get_name instead of get_gwrite. 
            #init_val = group_params.get_scalar (var.get_gwrite() ).get_value()
            init_val = group_params.get_scalar (var.get_name() ).get_value()
            return '\n' + var.get_c_type() + ' ' + var.get_gwrite() + ';'
        else:
            fill_method = group_params.get_array (var.get_gwrite() ).get_fill_method()

            return '\n' + var.get_c_type() + ' * ' + var.get_gwrite() + ';' 
    
    
    @staticmethod
    def get_initialization (var, group_params):
        dims = var.get_dimensions()
        if (dims == None):
            # I think this should be get_name instead of get_gwrite. 
            #init_val = group_params.get_scalar (var.get_gwrite() ).get_value()
            init_val = group_params.get_scalar (var.get_name() ).get_value()
            return '\n' + var.get_gwrite() + ' = ' + init_val + ';'
        else:
            fill_method = group_params.get_array (var.get_gwrite() ).get_fill_method()

            return '\n' + var.get_gwrite() + ' = (' + var.get_c_type() + '*) malloc (' + cFormatter.get_dim_str (var, '*') + ' * sizeof (' + var.get_c_type() + ') );\n' + cFormatter.get_fill (var, fill_method)



    @staticmethod
    def get_groupsize_code (group):
        groupsize_code_string = ''
        groupsize_code_string += '\n\n// Set the adios group size'
        groupsize_code_string += '\nadios_groupsize ='
        for v in group.get_vars():
            if (v.is_scalar() ):
                groupsize_code_string += ('\n                     %d +' % typeMapper.get_size (v.get_type() ) )
            else:
                groupsize_code_string += ('\n                     %d * ' % typeMapper.get_size (v.get_type() ) )
                
                for d in v.get_dimensions():
                    # need to check whether this is the timestep
                    groupsize_code_string += '(' + d + ') * '

                groupsize_code_string = groupsize_code_string.rstrip ('* ')

                groupsize_code_string += (' +')

        # remove the final +, and add the ;
        groupsize_code_string = groupsize_code_string.rstrip('+') + ';'

        groupsize_code_string += '\nadios_group_size (adios_handle, adios_groupsize, &skel_total_size);'

        return groupsize_code_string;




