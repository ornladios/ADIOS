#!/usr/bin/env python

import yaml


# A representation of a skel test plan that is stored in a yaml document
class skel_test_plan:
    def __init__ (self, filename):
        stream = file (filename, 'r')
        self.doc = yaml.load(stream)



    # Returns a list of test objects that implement the test plan
    def get_tests (self):

        tests = []
        for i in range (len (self.doc) ): # Assume each element of doc is a test
            tests.append (test (self.doc[i]['procs'], self.doc[i]['name'], self.doc[i]['method'], self.doc[i]['parameters'], self.doc[i]['decomp-procs']) )

        # This is the old way, convenient, but lacking in flexibility
        #tests = []
        #for i in range (len(self.doc['procs']) ):
        #    for j in range (len(self.doc['methods']) ):
        #        name = "T%s_%i" % (chr (ord('a') + j), self.doc['procs'][i])
        #        tests.append (test (self.doc['procs'][i], name, 
        #                      self.doc['methods'][j]['m'],
        #                      self.doc['methods'][j]['p']) )

        return tests        


class test:
    def __init__ (self, num_procs, subdir, method, parameters, decomp_procs):
        self.num_procs = num_procs
        self.subdir = subdir
        self.method = method
        self.parameters = parameters
        self.decomp_procs = decomp_procs

    def get_subdir (self):
        return self.subdir

    def get_method (self):
        return self.method

    def get_parameters (self):
        return self.parameters

    def get_decomp_procs (self):
        return self.decomp_procs

    def to_yaml (self):
        return "procs: %s\nmethod: %s\nparameters: %s\ndecomp_procs: %s" % (self.num_procs, self.method, self.parameters, self.decomp_procs)




def main(argv=None):
    b = skel_test_plan ("test_plan.yaml")

    print b.get_tests()


if __name__ == "__main__":
    main()
