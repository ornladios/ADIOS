import xml.dom.minidom

import skel_settings
import skel_have_adios_timing


class skelConfig:
    def __init__ (self, infile):
        doc = xml.dom.minidom.parse (infile)

        nodes = doc.childNodes
        if (nodes.length != 1):
            print 'malformed param file, should contain only a single skel-config element'
            raise SystemExit
        self.config_node = nodes[0]

        self.groups = []
        self.group_dict = {}

        for node in self.config_node.getElementsByTagName('adios-group'):
            this_group = group (node)
            self.groups.append (this_group)
            self.group_dict[this_group.get_name()] = this_group

        self.batches = []

        for node in self.config_node.getElementsByTagName('batch'):
            self.batches.append (batch (node) )

    def get_group (self, name):
        return self.group_dict[name]

    def get_groups (self):
        return self.groups

    def get_batches (self):
        return self.batches

    def get_application (self):
        return self.config_node.getAttribute ('application')

    def get_target (self):
        target = self.config_node.getAttribute ('target')

        if target == None or target == '':
            return "default"

        return target


class group:
    def __init__ (self, group_node):
        self.group_node = group_node

        self.scalars = []
        self.arrays = []

        self.scalar_dict = {}
        self.array_dict = {}

        for node in self.group_node.childNodes:
            if (node.localName == 'scalar'):
                this_scalar = scalar (node)
                self.scalars.append (this_scalar)
                self.scalar_dict[this_scalar.get_name()] = this_scalar
            if (node.localName == 'array'):
                this_array = array (node)
                self.arrays.append (this_array)
                self.array_dict[this_array.get_name()] = this_array

    def get_name (self):
        return self.group_node.getAttribute ('name')

    def get_scalars (self):
        return self.scalars

    def get_scalar (self, name):
        return self.scalar_dict[name]

    def get_arrays (self):
        return self.arrays

    def get_array (self, name):
        return self.array_dict[name]


class batch:
    def __init__ (self, batch_node):
        self.batch_node = batch_node

        self.tests = []
        for node in self.batch_node.getElementsByTagName ('test'):
            self.tests.append (test (node) )

    def get_tests (self):
        return self.tests

    def get_name (self):
        return self.batch_node.getAttribute ('name')

    def get_cores (self):
        return int (self.batch_node.getAttribute ('cores') )

    def get_walltime (self):
        return self.batch_node.getAttribute ('walltime')


class test:
    def __init__ (self, test_node):
        self.test_node = test_node
        self.measure = measure (self.test_node.getAttribute ('measure') )


    def get_type (self):
        return self.test_node.getAttribute ('type')

    def get_method_params (self):
        return self.test_node.getAttribute ('method_params')

    def get_group_name (self):
        return self.test_node.getAttribute ('group')

    def get_method (self):
        return self.test_node.getAttribute ('method')

    def get_rm (self):
        return self.test_node.getAttribute ('rm')

    def get_tags (self):
        return self.test_node.getAttribute ('tags')

    def get_steps (self):
        steps = self.test_node.getAttribute ('steps')
        if steps == None or steps == '':
            return '1'
        else:
            return steps

    def get_iterations (self):
        iter = self.test_node.getAttribute ('iterations')
        if iter == None or iter == '':
            return '1'
        else:
            return iter

    def get_ext (self):
        return self.test_node.getAttribute ('ext')

    def get_measure (self):
        return self.measure

# the measure object hides the logic of which code to emit for the various measurement options
class measure:
    def __init__ (self, measure_string):
        self.settings_dict = skel_settings.skel_settings().get_settings_dict()
        self.measure_string = measure_string

    def use_barrier_before_open (self):
        val = self.settings_dict["barrier_before_open"]
        if val =="no":
            return False

        return True

    def use_barrier_before_access (self):
        val = self.settings_dict["barrier_before_access"]
        if val =="yes":
            return True 

        return False

    def use_barrier_before_close (self):
        val = self.settings_dict["barrier_before_close"]
        if val =="yes":
            return True 

        return False 

    def use_barrier_after_close (self):
        val = self.settings_dict["barrier_after_close"]
        if val =="yes":
            return True 

        return False 

    def use_barrier_before_final_time (self):
        val = self.settings_dict["barrier_after_steps"]
        if val =="no":
            return False 

        return True

    def use_reduce (self):
        #val = self.settings_dict["reduce_times"]
        #if val =="no":
        #    return False 

        # Only reduced output is available from skel. Use adios timing mechanism for more detailed results

        return True

    def report_all (self):
        #val = self.settings_dict["report_all"]
        #if val =="yes":
        #    return True 

        # Only reduced output is available from skel. Use adios timing mechanism for more detailed results

        return False

    def use_adios_timing (self):
        val = self.settings_dict["use_adios_timing"]
        if val =="yes" and skel_have_adios_timing.skel_have_adios_timing == True:
            return True 

        return False

    def use_sleep_before_open (self):
        val = self.settings_dict["sleep_before_open"]
        if val =="yes":
            return True 

        return False 

class scalar:
    def __init__ (self, scalar_node):
        self.scalar_node = scalar_node

    def get_name (self):
        return self.scalar_node.getAttribute ('name')

    def get_value (self):
        return self.scalar_node.getAttribute ('value')


class array:
    def __init__ (self, array_node):
        self.array_node = array_node

    def get_name (self):
        return self.array_node.getAttribute ('name')

    def get_gwrite (self):
        return self.array_node.getAttribute ('gwrite')

    def get_fill_method (self):
        return self.array_node.getAttribute ('fill-method')
