import xml.dom.minidom
import type_mapper

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
    
        for node in self.config_node.getElementsByTagName ('adios-group'):
            self.adios_groups.append (adiosGroup (node) )

        for node in self.config_node.getElementsByTagName ('method'):
            self.methods.append (method (node) )

        for node in self.config_node.getElementsByTagName ('buffer'):
            # there should be only one of these... this code ignores all but the last one.
            self.buffer = buffer (node)

        # We are currently ignoring any analysis declarations


    def get_filename (self):
        return self.config_file_name

    def get_groups (self):
        return self.adios_groups

    def get_buffer (self):
        #return the buffer info
        print 'get_buffer is not yet implemented'

    def get_host_language (self):
        return self.config_node.getAttribute ('host-language')


class adiosGroup:

    def __init__ (self, group_node):
        self.group_node = group_node

        self.time_index = self.group_node.getAttribute ('time-index')

        self.vars = []
        self.vardict = {}
        self.vars_and_gwrites_and_attrs = []

        self.attrs = []
        self.attrdict = {}

        for node in self.group_node.childNodes:
            if node.localName == 'var':
                newvar = var (node, self, self.time_index)
                self.vars.append (newvar)
                #print 'Add to dict local var ['+newvar.get_fullpath()+']'
                self.vardict [newvar.get_fullpath()] = newvar
                self.vars_and_gwrites_and_attrs.append (newvar)
            #elif node.localName == 'attribute':
                #handle attribute

            elif node.localName == 'gwrite':
                self.vars_and_gwrites_and_attrs.append (gwrite (node) )

            elif node.localName == 'global-bounds':
                for gb_node in node.childNodes:
                    if gb_node.localName == 'var':
                        newvar = var (gb_node, self, self.time_index)
                        self.vars.append (newvar)
                        #print 'Add to dict global var ['+newvar.get_fullpath()+']'
                        self.vardict [newvar.get_fullpath()] = newvar
                        self.vars_and_gwrites_and_attrs.append (newvar)
                    elif gb_node.localName == 'gwrite':
                        self.vars_and_gwrites_and_attrs.append (gwrite (node) )
            
            elif node.localName == 'attribute':
                newattr = attr (node)
                self.attrs.append (newattr)
                self.attrdict [newattr.get_name()] = newattr
                self.vars_and_gwrites_and_attrs.append (newattr)


    # Returns the name of the group
    def get_name (self):
        return self.group_node.getAttribute ('name')

    # Returns a list of var objects for all of the variables in the group
    def get_vars (self):
        return self.vars

    # Returns the variable from this group with the specified name, or None
    def get_var (self, varfullpath):
        #print '          get_var('+varfullpath+')'
        if self.vardict.has_key (varfullpath):
            return self.vardict [varfullpath]

        return None


    # Returns a list containing all of the vars and gwrites and attributes in the same order
    # as was specified in the xml
    def get_ordered_contents (self):
        return self.vars_and_gwrites_and_attrs


class gwrite:
    def __init__(self, gwrite_node):
        self.gwrite_text = gwrite_node.getAttribute ('src')

    def get_src (self):
        return self.gwrite_text


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

    def get_path (self):
        path = self.var_node.getAttribute ('path')
        return path

    def get_name (self):
        name = self.var_node.getAttribute ('name')
        return name

    def get_fullpath (self):
        path = self.get_path()
        name = self.get_name()
        if (path == ''):
            fullpath = name
        elif (path[-1:] == '/'):
            fullpath = path + name
        else:
            fullpath = path + '/' + name
        return fullpath

    def get_gwrite (self):
        gw = self.var_node.getAttribute ('gwrite')

        if (gw == ''):
            gw = self.get_name()
    
        return gw

    def get_group (self):
        return self.group

    def get_c_type (self):
        return type_mapper.get_c_type (self.var_node.getAttribute ('type') )

    def get_type (self):
        return self.var_node.getAttribute ('type')

    def get_dimensions (self):
        if (self.var_node.getAttribute ('dimensions') == ''):
            return None
        else:
            # place the dimensions in a list and remove the time-index if it is there.
            dims = filter (lambda x : x != self.time_index, self.var_node.getAttribute ('dimensions').split(',') )
            cleandims = []
            #print '       get_dimensions of var '+self.get_fullpath()
            for d in dims:

                #print '          dim "'+str(d)+'"'
                if d.isdigit():
                    cleandims.append (d)
                    continue

                # Here we need to translate the variable name for this dimension (if it's a var) into the gwrite 
                # for that variable
                dim_var = self.get_group().get_var (d)				
                if dim_var != None:
                    #print '            dim var found, get name...'
                    d = dim_var.get_gwrite()
                #else:
                    #print '            dim var NOT found'
                    

                cleandims.append (d)
            return cleandims

    def is_scalar (self):
        return self.get_dimensions() == None


class attr:

    def __init__ (self, attr_node):
        self.attr_node = attr_node

    def get_name (self):
        return self.attr_node.getAttribute ('name')

