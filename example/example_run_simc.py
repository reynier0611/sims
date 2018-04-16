#
# run simc for Hall C proposal
#
# do simc calculation for a range of input files
#
# run this script from the directory where simc resides

import run_command as rc
import shutil as sh
import os
from LT.datafile import dfile

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
        

# select which type of calculation is done via the extension of
# input files used

# kinematics list
kin_list = 'kin_list.data'

# it  ratiation should be included set radiate to True
radiate = False
#radiate = True

change_nevents = 2000
#change_nevents = 50000

# local files, where you are running the code
#LOCALDIR = '/data/boeglin.1/user/deut/hi_pm/simc/'
LOCALDIR = '/data/boeglin.1/HallA/analysis/exp_results_q1/SIMC/WB_1_simc/'

KINDIR = LOCALDIR+'./'

# location of common input files
INFILE_DIR = '/data/boeglin.1/HallA/analysis/exp_results_q1/SIMC/q1_infiles_WB_new/' 
# directory with the extra file
EXTRA_DIR = '/data/boeglin.1/HallA/analysis/exp_results_q1/SIMC/q1_kin_offset_WB_new/'

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
SIMCDIR ='/data/boeglin.1/HallC/simc_gfortran.1/'
# use new simc
simc_command = SIMCDIR + './simc'

SIMCINDIR = './infiles/'
SIMCOUTDIR = './outfiles/'
SIMCRESDIR = './worksim/'


# create the correct scipt for the interactive input to simc
cc = open('./run_simc','w')
cc.write('current.data\n')
cc.close()


# open the kinematics file
kin_d = dfile(KINDIR+kin_list)

ext1 = 'data'

# read the header
if radiate:
    ext2 = 'rad'
    extra_file = 'extra_hydrogen_rad.data'
    extra_ext = '.dat'
else:
    ext2 = 'norad'
    extra_file = 'extra_hydrogen_norad.data'
    extra_ext = '_norad.dat'

for k in kin_d.data:
    # setup file names from kin_list.data file
    kin = k['kin'].strip()   # remove leading and trailing blanks
    name = 'q1_'+ kin
    calc_type = k['calc']
    # create the file name
    file_root = name
    file_name = file_root + '.' + calc_type + '.'+ ext2
    file = file_name + '.' + ext1
    offset_file = 'offset_' + kin + extra_ext 
    root_file = file_name + '.root'
    print 'calculate : ' + file
    print 'using offsets from : ' + offset_file
    print 'copied to : ' + extra_file
    print 'ext2:', ext2 , 'ext1:', ext1
    # copy current input file to current.data
    if change_nevents > 0:
        print "change number of events to ", change_nevents
        set_nevents(change_nevents, INFILE_DIR+file, SIMCINDIR+'current.data')
    else:
        sh.copy(INFILE_DIR+file, SIMCINDIR+'current.data')
    sh.copy(EXTRA_DIR+offset_file, SIMCINDIR+extra_file)
    #sh.copy(INDIR+'/defaults/'+extra_file, SIMCINDIR+'/.')
    out_file = LOGDIR + file+'.out'
    err_file = ERRDIR + file+'.err'
    print "Logging in : ", out_file
    print "Errors in : ", err_file
    # run simc using the script using the contents in run_simc as stdin input
    f_stdin = open('./run_simc')
    # ret = rc.run_command(SIMCDIR + './simc',out_file, err_file, stdin = f_stdin)
    ret = rc.run_command(simc_command,out_file, err_file, stdin = f_stdin)
    f_stdin.close()
    if ret != 0:
        print "Problem with : ", file
        continue
    else:
        print 'finished SIMC for ', file
    # copy the output files
    sh.copy(SIMCOUTDIR+'current.hist', OUTDIR+file_root+'.hist')
    sh.copy(SIMCOUTDIR+'current.geni', OUTDIR+file_root+'.geni')
    sh.copy(SIMCOUTDIR+'current.gen',  OUTDIR+file_root+'.gen')
    # create the root file
    normfac = get_normfac(SIMCOUTDIR+'current.hist')
    if (normfac < 0):
        print " cannot get normfac, skip ", name
        continue
    else:
        open('normfact.data', 'w').write("%r"%(normfac))
    command = SIMCDIR+'root/fmake_tree'
    ir = rc.run_command(command, LOGDIR+'out1', ERRDIR+'err1')
    if ir != 0:
        print "problem running " + command + " check out1 and err1 !"
    print 'rename root file to : ', root_file
    command = 'mv '+ SIMCRESDIR + 'simc.root ' + RESDIR + root_file
    ir = rc.run_command(command, LOGDIR+'out2', ERRDIR+'err2', shell = True)
    if ir != 0:
        print "problem running " + command + " check out2 and err2 !"
# do the next file
