import numpy as np
import matplotlib.pyplot as plt
from data_loader import *


density_map = {
    'gold' : None,
    'iron' : None,
    'titanium' : None,
}
data = {}
recursive_read(data, "data/", require = [0])

e_empty = ((5.486 * 0.86 + 5.443 * 0.127) / (0.86 + 0.127))

print (str(e_empty)+": "+str(e_empty)+" MeV")

ref = data[('empty', 0)]

e = lambda ch : e_empty / ref * ch

energies = {}
for metadata, entry in data.items():
    if metadata[0] != 'empty':
        energies[metadata[0]] = e(max_model(entry[1]))
        print (str(metadata[0])+": "+str(energies[metadata[0]])+" MeV")

gold_tracking = []
for foil, energy in energies.items():
    energy_loss = ref - energies[foil]
    if 'gold' in foil: 
        gold_tracking.append(energy_loss)
