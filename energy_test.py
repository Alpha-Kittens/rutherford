from energy_loss_2 import *

print ("beginning testing")
data = energy_data()
energies = get_energies(data, verbose = True, plot = True)
#thickness = get_thickness(energies, verbose = True)
#thiccness_dx('gold', energies['empty'].val, energies[('gold')].val)
#thiccness_dx('2gold', energies['empty'].val, energies[('2gold')].val)
thiccness_dx('4gold', energies['empty'].val, energies[('4gold')].val, 0.34490747, 0.16126968)

