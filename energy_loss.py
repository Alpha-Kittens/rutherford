import numpy as np
import matplotlib.pyplot as plt
from data_loader import *
from energy_fit import max_model
from data_processing import *
import math


#g/cm^3
density_map = {
    'gold' : 19.32,
    'iron' : 7.874,
    'titanium' : 4.5,
    '2gold' : 19.32,
    '3gold' : 19.32,
}
mass_map = {
    'gold' : 183448377e-3, # MeV
    'iron' : 520123850e-4, # average over isotopes
    #'titanium' : , # average over isotopes
    '2gold' : 183448377e-3,
    '3gold' : 183448377e-3,
}
I = {
    'gold' : 790e-6, #MeV
    '2gold' : 790e-6,
    '3gold' : 790e-6,
    'iron' : 286e-6,
    'titanium' : 233e-6
}
number = {
    'gold' : 1,
    'iron' : 1,
    'titanium' : 1,
    '2gold' : 2,
    '3gold' : 3,
}
c = 299792458 # m/s
m_a = 3.7273794066e3 #MeV
m_e = 0.51099895000 #MeV
K = 0.307075 # MeV / mol * cm^2
v = lambda b : b * c
g = lambda b : 1 / math.sqrt(1-b**2)
W_max = lambda b : (2*m_e*b**2*g(b)**2) / (1+2*g(b)*m_e/m_a + (m_e/m_a)**2)
dEdx = lambda b, foil : K * 2**2 * element_map[foil] / mass_map[foil] / b**2 * (1/2 * math.log(2*m_e*b**2*g(b)**2*W_max(b)/I[foil]**2) - b**2) * density_map[foil]
b = lambda E : math.sqrt(1 - m_a**2 / E**2)

def thiccness(foil, incident, exiting):
    E = incident
    x = 0
    dt = 1e-8
    while E > exiting:
        beta = b(E)
        x += v(beta) * dt
        E -= dEdx(b, foil) * v(beta) * dt
    return x


if __name__ == '__main__':
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium"])
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
    #Z alpha is 2

    for foil, energy in energies.items():
        energy_loss = ref_e - energy #sum(ref_e, minus(energy))
        
        thickness = energy_loss / (K * 2**2 * element_map[foil] / mass_map[foil])
        #e(ref_ch) - energies[foil]
        if 'gold' in foil: 
            gold_tracking.append(energy_loss / number[foil])
        report (str(foil) + " energy loss", energy_loss, "MeV")
        #print (str(foil)+" energy loss: "+str(energy_loss[0])+" MeV ± "+str(energy_loss[1])+" MeV")
    #print ([entry.val for entry in gold_tracking])
    report ("Energy loss per gold foil: ", weighted_average(gold_tracking), "MeV")


