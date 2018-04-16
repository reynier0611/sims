# convert deep.data file to a root file
#
# tun this script
#
import subprocess as sub

import shutil as sh
import os.path as p


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



# SIMC root directory
# SIMDIR = '/data/boeglin.1/HallC/simc_gfortran.1/'
SIMDIR = './'

# simc input file directory
INDIR = 'infiles/'
# simc output parameters
OUTDIR = 'outfiles/'
# simc data file
DATDIR = 'worksim/'

# final results in
RESDIR = './output/'

#----------------------------------------------------------------------
# get normfac
def get_normfac(file):
    data = open(file).readlines()
    found = False
    value = -1.
    for d in data:
        if d.find('normfac') < 0: continue
        fields = d.split('=')
        try:
            value = float( fields[1])
            found = True
        except:
            print 'cannot extract normfac'
        break
    return value
#----------------------------------------------------------------------
# default simc output file name
in_name = 'current.data'
name,ext = p.splitext(in_name)
data_file = SIMDIR + DATDIR + name + ext
hist_file = SIMDIR + OUTDIR + name + '.hist'
root_file = SIMDIR + DATDIR + name + '.root'
root_res_file = SIMDIR + RESDIR + name + '.root'
normfac = get_normfac(hist_file)
print 'convert : ' + data_file + ' to ' + root_file + ' with normfac = ', normfac
# open temp data file
dd = open('./make_root_tmp.dat','w')
dd.write( "%r\n"%(normfac))
dd.close()
print 'temporary file make_root_tmp.dat created converto to root file'
command = 'cat ' + data_file + ' >> ./make_root_tmp.dat'
run_command(command, 'out1', 'err1')
command = SIMDIR+'root/make_tree < make_root_tmp.dat'
run_command(command,'out2','err2')
print 'rename root file to : ', root_res_file
command = 'mv simc.root ' + root_res_file
run_command(command,'out3','err3')
#command = 'rm make_root_tmp.dat'
#rc.run_command(command,'out4','err4')
