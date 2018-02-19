#
# run simc for E01020
#
# do simc calculation for a range of input files
#
# run this script from the directory where simc resides

import run_command as rc
import shutil as sh


# select which type of calculation is done via the extension of
# input files used

ext1 = '.pwia'
ext2 = '.rad.20'
ext3 = '.data'

SIMCDIR = '/data/boeglin.1/HallC/simc_gfortran.1/'
INDIR = './infiles/'
OUTDIR = './outfiles/'
RESDIR = './worksim/'

# create the correct scipt for the interactive input to simc
cc = open('./run_simc','w')
cc.write('current.data\n')
cc.close()

for name in names:
    # create the file name
    file_root = name + ext1 + ext2
    file = file_root + ext3
    print 'calculate : '+file
    # copy curren input file to current.data
    sh.copy(INDIR+file, INDIR+'current.data')
    out_file = file+'.out'
    err_file = file+'.err'
    # run simc using the script using the contents in ruN_simc as stdin input
    f_stdin = open('./run_simc')
    ret = rc.run_command('./simc',out_file, err_file, stdin = f_stdin)
    f_stdin.close()
    if ret != 0:
        print "Problem with : ", file
        continue
    # copy the output files
    sh.copy(OUTDIR+'current.hist', OUTDIR+file_root+'.hist')
    sh.copy(OUTDIR+'current.geni', OUTDIR+file_root+'.geni')
    sh.copy(OUTDIR+'current.gen', OUTDIR+file_root+'.gen')
    sh.copy(SIMDIR+'current.data', SIMDIR+file_root+'.data')
# do the next file
