import os
import sys

def main(argv=None):
    cmd = "make\n"
    cmd += "./genbp\n"
    cmd += "./bp_read_c testbp_c.bp\n"
    cmd += "./bp_read_f\n"
    os.system(cmd)

    return
if __name__=="__main__":
        main()
