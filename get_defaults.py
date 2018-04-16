#
# create namelist datafile
#

def make_type(x, type):
    if type == 'int':
        return int(x)
    if type == 'double':
        return float(x)
    if type == 'doublearray':
        return float(x)
# end

d = open('nml_defaults.data').readlines()
o = open('namelist.data','w')

found_set = False
first = True
for l in d:
    blank = (len(l.split()) == 0)
    if l.find('*') >= 0:
        if not first :
            o.write(' / \n')
        first = False
        found_set = True
        name = l.split()[1]
        print 'namelist : ', name
        lo = ' &' + name
        o.write(lo + '\n')
    elif found_set and ( not blank):
        data = l.split(',')
        type = l.split('regparm')[1].split('(')[0]
        if len(data) == 3:
            var_name = data[-2:-1][0]
            var_val = data[-1:][0].split(')')[0]
            # lc = '!   type: '+ type
            lo = '    ' + var_name + ' = '+ var_val
            # o.write(lc + '\n')
            o.write(lo + '\n')
        elif len(data) == 4:
            var_name = data[-3:-2][0]
            var_val = (data[-2:-1][0], data[-1:][0].split(')')[0])
            # lc = '!   type: '+ type
            lo = '    ' + var_name + ' = '+ var_val[0] + ', ' + var_val[1]
            # o.write(lc + '\n')
            o.write(lo + '\n')
        print 'variable : ', var_name, var_val, type
# done
o.write(' / \n')
o.close()
