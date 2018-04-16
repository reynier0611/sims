#
# run simc for a list of input files
#
# do simc calculation for a range of input files
#
# run this script from the directory where simc resides

import run_command as rc
import shutil as sh

# the last Q2 set will be run !

# Q2 = 3.5
names=\
['pm_0.5_3.5_1.35_11.00',\
 'pm_0.6_3.5_1.35_11.00',\
 'pm_0.7_3.5_1.35_11.00',\
 'pm_0.8_3.5_1.35_11.00',\
 'pm_0.9_3.5_1.35_11.00',\
 'pm_1.0_3.5_1.35_11.00',\
 'pm_1.1_3.5_1.35_11.00',\
 'pm_1.2_3.5_1.35_11.00',\
 'pm_1.3_3.5_1.35_11.00',\
 'pm_1.4_3.5_1.35_11.00']

# Q2 = 5.0
names=\
['pm_0.5_5.0_1.35_11.00',\
 'pm_0.6_5.0_1.35_11.00',\
 'pm_0.7_5.0_1.35_11.00',\
 'pm_0.8_5.0_1.35_11.00',\
 'pm_0.9_5.0_1.35_11.00',\
 'pm_1.0_5.0_1.35_11.00',\
 'pm_1.1_5.0_1.35_11.00',\
 'pm_1.2_5.0_1.35_11.00']


# Q2 = 4.25 now q2
names=\
['pm_0.1_4.2_1.00_11.00',\
 'pm_0.2_4.2_1.00_11.00']

# Q2 = 3.5 low q2
names=\
['pm_0.1_3.5_1.00_11.00',\
 'pm_0.2_3.5_1.00_11.00']

# Q2 = 4.25
names=\
['pm_0.5_4.2_1.35_11.00',\
 'pm_0.6_4.2_1.35_11.00',\
 'pm_0.7_4.2_1.35_11.00',\
 'pm_0.8_4.2_1.35_11.00',\
 'pm_0.9_4.2_1.35_11.00',\
 'pm_1.0_4.2_1.35_11.00',\
 'pm_1.1_4.2_1.35_11.00',\
 'pm_1.2_4.2_1.35_11.00']

# special calculations
names=\
['pm_0.5_4.2_1.35_11.00',\
 'pm_1.0_4.2_1.35_11.00']
# special calculations
names=\
['werner_h2_example']

# select which type of calculation is done via the extension of
# input files used

ext1 = '.pwia'
# norad: no radiation
# rad including radiation

# ext2 = '.norad'
# ext2 = '.rad'

ext2 = '.rad.20'
ext3 = '.data'

INDIR = './infiles/'
OUTDIR = './outfiles/'
SIMDIR = './worksim/'

# create the correct scipt for the interactive input to simc
cc = open('./run_simc','w')
cc.write('current.data\n')
cc.close()

for name in names:
    # create the file name
    file_root = name + ext3
    file = name + ext3
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
