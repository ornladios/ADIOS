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
    parser.add_argument ('-R', dest='generate_R', action="store_true", default=False, help='Extract into a form to be read by R.')
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
    


# Hilde's function for extracting to a format understood by the R script
def extract_R (skel_output, select, ranks, iteration):
 
    doc = xml.dom.minidom.parse (skel_output)
    #if dest == None or dest == '':
    #    outfile = sys.stdout
    #else:
    outfile = open (skel_output+".txt", 'w')
    #runval=dest.split('.')[0].split('_')[1]
    #outfile.write ("runid = " +runval +'\n')

    if select == None or select == '':
        select = '*'
 


    
    method = doc.getElementsByTagName('adios_timing')[0].getAttribute('method')

    keys = doc.getElementsByTagName('adios_timing')[0].getAttribute('keys').split(',')
    tempkeys = []
    for key in keys:
        tempkeys.append (key.strip())

    keys = tempkeys

    numcores= doc.getElementsByTagName('adios_timing')[0].getAttribute('cores')

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
        header = header +' '

    header = "Iteration "+"Core "+header.rstrip (' ')
    outfile.write (method + '\n')
    outfile.write (numcores +'\n')
    outfile.write (header + '\n')
    
    # Replacing iter with the iteration argument passed in
    iter = 0

    for st in doc.getElementsByTagName ('adios_timing'):

        #Print the data
        for proc in st.getElementsByTagName('proc'):
            if proc.getAttribute ('id') == '0':
                iter = iter +1
            data = ''
            core = proc.getAttribute ('id')
            vals = proc.getAttribute ('vals').split(',')
            for field in selected_fields:
                data = data + vals[keys.index(field)].strip(',') + ' '
            data = str(iteration)+" "+ core +" "+ data.rstrip (' ') + '\n'
            outfile.write (data)
        if ranks == 'all':


            
            break
    
    outfile.close ()

def parse_iteration (filename):
    #assume filename ends with .xml
    if not filename.endswith (".xml"):
        print "Warning: filename does not meet expectations, should end with .xml"

    filename = filename [:-4]

    iteration = filename.rsplit ("_", 1)[1]

    print iteration

    return iteration

def main ():


    args = parse_command_line()
    if (args.generate_R):
        iteration = parse_iteration (args.skel_output)
        extract_R (args.skel_output, args.select, args.ranks, iteration)
    else:    
        extract (args.skel_output, args.dest, args.select, args.ranks)



if __name__ == "__main__":
    main()


