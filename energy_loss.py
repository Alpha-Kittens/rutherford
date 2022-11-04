import numpy as np
import matplotlib.pyplot as plt
from data_loader import *
from energy_fit import max_model
from data_processing import *


#g/cm^3
density_map = {
    'gold' : 19.32,
    'iron' : 7.874,
    'titanium' : 4.5,
}
number = {
    'gold' : 1,
    'iron' : 1,
    'titanium' : 1,
    '2gold' : 2,
    '3gold' : 3,
}
if __name__ == '__main__':
    data = {}
    recursive_read(data, "data", require = [0])
    #print (data.keys())
    #print (data.keys())
    #print (data.keys())

    e_empty = ((5.486 * 0.86 + 5.443 * 0.127) / (0.86 + 0.127))

    print ("Expected beam energy: "+str(e_empty)+" MeV")

    ref_data = iterationless(data)[('empty', 0)]
    #print (ref)

    ref_chs = [resultify(max_model(r[1])) for r in ref_data]
    ref = resultify(weighted_average(ref_chs))
    #ref_ch, ref_unc = ref

    #print (ref_e)

    e = lambda ch : Result(e_empty) / ref * ch
    #print (ref)
    #print(len(ref))
    #print (e2(ref))
    #print (e(ref))
    ref_e = e(ref)

    energies = {}
    for metadata, entry in data.items():
        if metadata[0] != 'empty':
            energies[metadata[0]] = e(resultify(max_model(entry[1])))
            report (metadata[0], energies[metadata[0]], "MeV")
            #print (str(metadata[0])+": "+str(energies[metadata[0]][0])+" MeV ± "+str(energies[metadata[0]][1]))

    gold_tracking = []
    for foil, energy in energies.items():
        energy_loss = ref_e - energy #sum(ref_e, minus(energy))
        #e(ref_ch) - energies[foil]
        if 'gold' in foil: 
            gold_tracking.append(energy_loss / number[foil])
        report (str(foil) + " energy loss", energy_loss, "MeV")
        #print (str(foil)+" energy loss: "+str(energy_loss[0])+" MeV ± "+str(energy_loss[1])+" MeV")
    #print ([entry.val for entry in gold_tracking])
    report ("Energy loss per gold foil: ", weighted_average(gold_tracking), "MeV")


