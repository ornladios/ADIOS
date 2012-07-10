#!/usr/bin/env python

import sys
import xml.dom.minidom

#import skel_settings


#def parse_command_line():
#    parser = argparse.ArgumentParser (description='Merge two skel output files into one, keeping a single document root.')
#    parser.add_argument ('f1', metavar='file1', help='destination file')
#    parser.add_argument ('f2', metavar='file2', help='file to merge into destination file')

#    return parser.parse_args()





def concatenate (f1, f2):
    # Read both files, merge, write to f1 

    doc1 = xml.dom.minidom.parse (f1)
    doc2 = xml.dom.minidom.parse (f2)

    timing_to_merge = doc2.getElementsByTagName ('adios_timing')[0]
    doc1.documentElement.appendChild (timing_to_merge)

    file1 = open (f1, 'w')
    file1.write (doc1.toxml() )



def main ():

    #skel_settings.create_settings_dir_if_needed()

    #args = parse_command_line()

    concatenate (sys.argv[1], sys.argv[2])



if __name__ == "__main__":
    main()


