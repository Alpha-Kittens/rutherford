import numpy as np
import matplotlib.pyplot as plt
from data_loader import *
from energy_fit import max_model
from data_processing import *


density_map = {
    'gold' : None,
    'iron' : None,
    'titanium' : None,
}
data = {}
recursive_read(data, "data", require = [0])
#print (data.keys())
#print (data.keys())
#print (data.keys())

e_empty = ((5.486 * 0.86 + 5.443 * 0.127) / (0.86 + 0.127))

print ("Expected beam energy: "+str(e_empty)+" MeV")

ref_data = iterationless(data)[('empty', 0)]
#print (ref)

ref_chs = [max_model(r[1], plot = False) for r in ref_data]
ref = weighted_average(ref_chs)
ref_ch, ref_unc = ref

#print (ref_e)

en = lambda ch : e_empty / ref_ch * ch
e_unc = lambda ch, ch_unc : en(ch) * np.sqrt((ref_unc / ref_ch)**2 + (ch_unc/ch)**2)
e2 = lambda eval : eval[0] + eval[1]
e = lambda eval : (en(eval[0]), e_unc(eval[0], eval[1]))
#print (ref)
#print(len(ref))
#print (e2(ref))
#print (e(ref))
ref_e = e(ref)

energies = {}
for metadata, entry in data.items():
    if metadata[0] != 'empty':
        energies[metadata[0]] = e(max_model(entry[1], plot = False))
        report (metadata[0], energies[metadata[0]], "MeV")
        #print (str(metadata[0])+": "+str(energies[metadata[0]][0])+" MeV ± "+str(energies[metadata[0]][1]))

gold_tracking = []
for foil, energy in energies.items():
    energy_loss = sum(ref_e, minus(energy))
    #e(ref_ch) - energies[foil]
    if 'gold' in foil: 
        gold_tracking.append(energy_loss)
    report (str(foil) + " energy loss", energy_loss, "MeV")
    #print (str(foil)+" energy loss: "+str(energy_loss[0])+" MeV ± "+str(energy_loss[1])+" MeV")


