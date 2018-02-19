# convert deep.data file to a root file
#
# tun this script
#
import run_command as rc
import shutil as sh
import os.path as p

# SIMC root directory
SIMDIR = '/data/boeglin.1/HallC/simc_gfortran.1/'

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
print 'got normfact : ', normfac
# open normfact data file
dd = open('./normfact.data','w')
dd.write( "%r\n"%(normfac))
dd.close()

