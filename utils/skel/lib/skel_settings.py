#!/usr/bin/env python

import sys
import os.path
import shutil


class skel_settings:
    def __init__(self):
        self.settings_dict = {}

        # Parse the settings file. Ignore comments. Process name=value pairs.
        settings_file_name = '~/.skel/settings'
        settings_file = open (os.path.expanduser (settings_file_name) )
        for line in settings_file:
            line = line.strip()
            if line.startswith ('#'):
                continue
            if line == '':
                continue
            split_line = line.split('=')
            if not len (split_line) == 2:
                print 'Malformed configuration line: ' + line
                print 'Ignoring'
                continue
            self.settings_dict[split_line[0]] = split_line[1]



#    def get_skel_home (self):
#        return self.settings_dict['skel_home']

    def get_deploy_dir (self):
        return self.settings_dict['deploy_dir']


    def get_submit_target (self):
        t = self.settings_dict['submit_target']
        if t == None or t == '':
            return "sith"
        else:
            return t

    def get_account (self):
        return self.settings_dict['account']

    def get_settings_dict (self):
        return self.settings_dict




def create_settings_dir_if_needed():
    skel_settings_dir_name = os.path.expanduser ('~/.skel')
    if not os.path.exists (skel_settings_dir_name):
        bindir = os.path.dirname (sys.argv[0])
        shutil.copytree (bindir + '/../etc/skel', skel_settings_dir_name)
        print 'Created ' + skel_settings_dir_name


def main(argv=None):
    create_settings_dir_if_needed()


if __name__ == "__main__":
    main()



