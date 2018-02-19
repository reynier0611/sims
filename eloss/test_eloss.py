# test eloss calculation and correction
#
# testing the calculation of energy losses and corrections
#
import eloss as el
import numpy as np
import matplotlib.pyplot as pl

# energy
def E(p,m):
    return np.sqrt(p**2 + m**2)
#
proton = 938.5
electron = 0.511

nevents = 10000

l = 15.
dens = 0.167
z = 1.
a = 2.
ppart = 600.
M = proton
type = 1


epart = E(ppart, M)

de = []  # random

de0 = el.get_eloss(l, dens, z, a, epart , M, type) 

dem = el.eloss_results.eloss_mean # mean
demp = el.eloss_results.eloss_mp  # mostprobable

for i in range(nevents):
    de.append( el.get_eloss(l, dens, z, a, epart, M, type) )
#
de = np.array(de)
de_res = pl.hist( de, range = (2., 12), bins = 51)
de_res_mp = pl.hist( de-demp, range = (-2., 12), bins = 51)
de_res_m = pl.hist( de-dem, range = (-2., 12), bins = 51)

    
