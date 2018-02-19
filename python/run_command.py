# run the command in string cmd and put the output into a file name stored in
# out_file  and the error output into a file stored in err_file
#
# if the program run reads from stdin (e.g. and interactive program)
# pyt the answers intot a file, and open it and pass the file from 
# my_file = open('argument_file)
# with the                           stdin=my_file 
# keyword argument to run_command
#

import subprocess as sub
import pdb
import sys
#------------------------------------------------------------
def run_command( cmd, out_file, err_file, **kwargs):
   output_f = open(out_file,"w")
   error_f =  open(err_file,"w")
   p = sub.Popen(cmd, \
                 shell=True, \
                 stdout=output_f, \
                 stderr=error_f, \
                 close_fds=True, \
                    **kwargs)
   p.communicate(None)
   error_f.close()
   output_f.close()
   return p.returncode
#------------------------------------------------------------
