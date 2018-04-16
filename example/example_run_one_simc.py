#
# run simc for one example kinematic
#
# run this script from the directory where you ran setup.sh before

import run_command as rc
import shutil as sh
import os
import sys

#----------------------------------------------------------------------
# get normfac
def get_normfac(file):
    data = open(file).readlines()
    value = -1.
    for d in data:
        if d.find('normfac') < 0: continue
        fields = d.split('=')
        try:
            value = float( fields[1])
        except:
            print 'cannot extract normfac'
        break
    return value
#----------------------------------------------------------------------
# change the value for number of successes in the input file (fin) and
# write the new version to fout
def set_nevents(n, fin, fout):
    d = open(fin).readlines()
    for i,l in enumerate(d):
        # find the ngen line
        if l.find('ngen')>=0:
            f = l.split('=')
            new_line = f[0] + ' =  {}'.format(n) + ' ! ' + f[1].split()[-1] + '\n'
            d[i] = new_line
    open(fout, 'w').writelines(d)

#----------------------------------------------------------------------
        

# input files used

# main input file and its extension
# in_file_root = 'werner_h2_example'
in_file_root = 'werner_h2_example_laget'
in_file_ext = '.data'


extra_info_file = 'extra_hydrogen_norad.data'

change_nevents = 50000
#change_nevents = 50000

# local files, where you are running the code

LOCALDIR = './'

# location of common input files
INFILE_DIR = './infiles/' 
# directory with the extra file
EXTRA_DIR = './infiles/'

# check directory structure

INDIR =  LOCALDIR+'./infiles/'
if not os.path.exists(INDIR):
    print INDIR + ' does not exist, will create it !'
    os.mkdir(INDIR)

RESDIR = LOCALDIR+'./worksim/'
if not os.path.exists(RESDIR):
    print RESDIR + ' does not exist, will create it !'
    os.mkdir(RESDIR)

OUTDIR = LOCALDIR+'./outfiles/'
if not os.path.exists(OUTDIR):
    print OUTDIR + ' does not exist, will create it !'
    os.mkdir(OUTDIR)

ERRDIR = LOCALDIR+'./err/'
if not os.path.exists(ERRDIR):
    print ERRDIR + ' does not exist, will create it !'
    os.mkdir(ERRDIR)

LOGDIR = LOCALDIR+'./log/'
if not os.path.exists(LOGDIR):
    print LOGDIR + ' does not exist, will create it !'
    os.mkdir(LOGDIR)

# SIMC locations, this should not change
SIMCDIR ='/data/boeglin.1/HallC/GitHub/deut_simc/'
# use new simc
simc_command = SIMCDIR + './simc'

SIMCINDIR = './infiles/'
SIMCOUTDIR = './outfiles/'
SIMCRESDIR = './worksim/'


# create the correct scipt for the interactive input to simc
cc = open('./run_simc','w')
cc.write('current.data\n')
cc.close()

# common extra info file
extra_file = 'extra_hydrogen.data'


in_file = in_file_root + in_file_ext
root_file = in_file_root + '.root'

print 'calculate : ' + in_file
print 'using extra info from : ' + extra_info_file
print 'copied to : ' + extra_file

# copy current input file to current.data
if change_nevents > 0:
    print "change number of events to ", change_nevents
    set_nevents(change_nevents, INFILE_DIR+in_file, SIMCINDIR+'current.data')
else:
    sh.copy(INFILE_DIR+in_file, SIMCINDIR+'current.data')
# setup extra file
sh.copy(EXTRA_DIR+extra_info_file, SIMCINDIR+extra_file)
# setup error and log files
out_file = LOGDIR + in_file_root+'.out'
err_file = ERRDIR + in_file_root+'.err'
print "Logging in : ", out_file
print "Errors in : ", err_file
# run simc using the script using the contents in run_simc as stdin input
f_stdin = open('./run_simc')
# ret = rc.run_command(SIMCDIR + './simc',out_file, err_file, stdin = f_stdin)
ret = rc.run_command(simc_command,out_file, err_file, stdin = f_stdin)
f_stdin.close()
if ret != 0:
    print "Problem with : ", in_file
    sys.exit(-1)
else:
    print 'finished SIMC for ', in_file
# copy the output files
sh.copy(SIMCOUTDIR+'current.hist', OUTDIR+in_file_root+'.hist')
sh.copy(SIMCOUTDIR+'current.geni', OUTDIR+in_file_root+'.geni')
sh.copy(SIMCOUTDIR+'current.gen',  OUTDIR+in_file_root+'.gen')
# create the root file
normfac = get_normfac(SIMCOUTDIR+'current.hist')
if (normfac < 0):
    print " cannot get normfac, skip ", in_file_root
    sys.exit(-2)
else:
    open('normfact.data', 'w').write("%r"%(normfac))
command = SIMCDIR+'root/fmake_tree'
ir = rc.run_command(command, LOGDIR+'out1', ERRDIR+'err1')
if ir != 0:
    print "problem running " + command + " check out1 and err1 !"
    sys.exit(-3)
print 'rename root file to : ', root_file
command = 'mv '+ SIMCRESDIR + 'simc.root ' + RESDIR + root_file
# you might need to use the commen below
# ir = rc.run_command(command, LOGDIR+'out2', ERRDIR+'err2', shell = True)
ir = rc.run_command(command, LOGDIR+'out2', ERRDIR+'err2')
if ir != 0:
    print "problem running " + command + " check out2 and err2 !"
    sys.exit(-4)
# all done
print "------------- all finished ---------------------"
