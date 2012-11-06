#!/usr/bin/env python

import argparse
import os
import os.path
import shutil
import sys
import xml.dom.minidom

#import skel_settings


def parse_command_line():
    parser = argparse.ArgumentParser (description='Merge two skel output files into one, keeping a single document root.')
    parser.add_argument ('f1', metavar='file1', help='destination file')
    parser.add_argument ('f2', metavar='file2', help='file to merge into destination file')

    parser.add_argument ('--app', help='The name of the application corresponding to the result being concatenated, e.g. gts.')
    parser.add_argument ('--tags', help='A comma delimited list of name:value pairs to add to the result being concatenated to facilitate analysis')

    return parser.parse_args()





def concatenate (f1, f2, app, tags):
    # Read both files, merge, write to f1 

    doc1 = xml.dom.minidom.parse (f1)
    doc2 = xml.dom.minidom.parse (f2)

    timing_to_merge = doc2.getElementsByTagName ('adios_timing')[0]
    timing_to_merge.setAttribute ('app', app)
    timing_to_merge.setAttribute ('tags', tags)
    doc1.documentElement.appendChild (timing_to_merge)

    file1 = open (f1, 'w')
    file1.write (doc1.toxml() )






def move (f1, f2, app, tags):
    doc = xml.dom.minidom.parse (f2)

    timing_element = doc.getElementsByTagName ('adios_timing')[0]
    timing_element.setAttribute ('app', app)
    timing_element.setAttribute ('tags', tags)

    file1 = open (f1, 'w')
    file1.write (doc.toxml() )

    os.remove (f2)


def main ():

    #skel_settings.create_settings_dir_if_needed()

    args = parse_command_line()

    # If f1 does not exist, just move f2 to f1
    if not os.path.isfile (args.f1):
        #shutil.move (args.f2, args.f1) # Note args are (src, dest)
        move (args.f1, args.f2, args.app, args.tags)
    else:
        concatenate (args.f1, args.f2, args.app, args.tags)
        os.remove (args.f2)


if __name__ == "__main__":
    main()


