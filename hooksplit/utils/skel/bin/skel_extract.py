#!/usr/bin/env python

import argparse
import xml.dom.minidom
import sys

#import skel_settings


def parse_command_line():
    parser = argparse.ArgumentParser (description='Extract data from specified skel output file')
    parser.add_argument ('skel_output', metavar='Skel result', help='file to extract data from')
    parser.add_argument ('-d', dest='dest', help='file to place the results in')
    parser.add_argument ('-s', dest='select', help='measurement(s) to extract (default = *)')
    parser.add_argument ('-r', dest='ranks', default='0', help='Use -r all to include measurements from all ranks for just one iteration.')
    return parser.parse_args()




# By default, print all data for rank 0 as one line, repeat for each iteration
def extract (skel_output, dest, select, ranks):

    doc = xml.dom.minidom.parse (skel_output)
    if dest == None or dest == '':
        outfile = sys.stdout
    else:
        outfile = open (dest, 'w')

    if select == None or select == '':
        select = '*'

    keys = doc.getElementsByTagName('adios_timing')[0].getAttribute('keys').split(',')
    tempkeys = []
    for key in keys:
        tempkeys.append (key.strip())

    keys = tempkeys

    if select == '*':
        selected_fields = keys
    else:
        selected_fields = select.split(',')

    # check the selected fields
    for field in selected_fields:
        if not field in keys:
            print 'Invalid selection, field ' + field
            return

    #Print the header
    header = ''
    for field in selected_fields:
        header = header + field
        header = header +','

    header = header.rstrip (',')
    outfile.write (header + '\n')

    for st in doc.getElementsByTagName ('adios_timing'):

        #Print the data
        for proc in st.getElementsByTagName('proc'):
            if proc.getAttribute ('id') == '0' or ranks == 'all':  
                data = ''
                vals = proc.getAttribute ('vals').split(',')
                for field in selected_fields:
                    data = data + vals[keys.index(field)].strip() + ','
                data = data.rstrip (',') + '\n'
                outfile.write (data)
        if ranks == 'all':
            break
    
    outfile.close () 
    



def main ():


    args = parse_command_line()
    extract (args.skel_output, args.dest, args.select, args.ranks)



if __name__ == "__main__":
    main()


