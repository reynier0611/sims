#!/usr/bin/env python
#
# create namelist datafile from original CTP parameter file
#
# note 'begin' and 'end' statements must start at the beginning 
# of the line
#
import sys
import getopt
import pdb

logicals = [
'MC_SMEAR',
'DOING_PHSP',
'DOING_KAON',
'DOING_PION',
'DOING_DELTA',
'DOING_SEMI',
'DOING_HPLUS',
'DOING_RHO',
'DOING_DECAY',
'MC_SMEAR',
'HARD_CUTS',
'USING_RAD',
'DOING_PHSP',
'USING_ELOSS',
'CORRECT_ELOSS',
'CORRECT_RASTER',
'USING_CIT_GENERATION',
'USING_COULOMB',
'USE_OFFSHELL_RAD',
'DO_FERMI',
'USING_TGT_FIELD']

output = sys.stdout

def is_logical(x):
   item = x.strip()
   # check if the item is in the list
   try:
      logicals.index( item.upper() ) # use upper case
   except:
      is_logical = False
      return is_logical
   is_logical = True
   return is_logical


#------------------------------------------------------------
def usage():   # describe the usage of this program
   print "\nusage: ctp2nml.py [-options] ctp_file \n"
   print "options:      -h,? this message "
   print "              -o   output  file (def: stdout) "

#------------------------------------------------------------

#------------------------------------------------------------
# commandline arguments

args = sys.argv[1:]  # argument list w/o program name
#------------------------------------------------------------
# handle options
options = "o:h?"  # the ones following the colon require an argument

try:
   options, arguments = getopt.getopt(args, options)

except getopt.GetoptError, (msg, opt) :   
   
   # print help information and exit:
   print "ERROR in option : ",  opt, " reason: ", msg
   usage()
   sys.exit(1)

# handle the option values
for o, v in options:

    if o == "-h" or o == "-?":
        # help
        usage()
        sys.exit(0)
    if o == "-o":
        output = open(v,'w')
        
#------------------------------------------------------------

if (len(arguments) != 1):
     print 'you need to enter the name of the CTP file ! \n'
     usage()
     sys.exit(-1)

filename = arguments[0]
d = open(filename).readlines()

found_set = False
first = True
output.write('! This is a namelist file created from :' + filename+' \n') 
for l in d:
    blank = (len(l) == 0)
    if blank: continue
    # replace all ; with the comment character !
    ll = l.replace(';','!')
    # handle the begin param: extract the namelist name and write the data
    if ll.find('begin') == 0:
        line = ' &'
        found_set = True
        f = ll.split()
        name = f[2].upper()
        line += name + '\n'
        output.write(line)
        continue
    # handle the end param part
    if ll.find('end') == 0:
        found_set = False
        output.write('/ \n')
        continue
    # extract the variable assignments including comments
    if found_set and (not blank) :
        line = '   '
        comment = ''
        f = ll.split('!')
        if len(f) > 1:
            comment = '!' + f[1][:-1]
        fv = f[0].split('=')
        var_name = fv[0].replace('.','%')
        if len(fv) > 1:
            var_val = fv[1]
            # check if it is a logical variable
            if is_logical(var_name):
               val = int( float(var_val))
               if val == 0 : 
                  var_val = 'F'
               elif val == 1 : 
                  var_val = 'T'
               else:
                  print 'problem with logical : ', var_name, var_val, val
                  var_val = 'F'
            line += var_name + ' = ' + var_val + comment + '\n'
        else:
            # this is for lines that only contain comments
            line += var_name + comment + '\n'
        output.write(line)
        continue
    output.write(ll)
# done
if output != sys.stdout:
    output.close()
