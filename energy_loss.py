import numpy as np
import matplotlib.pyplot as plt
from data_loader import *
from energy_fit import max_model
from data_processing import *
import math
from plots import *


#g/cm^3
density_map = {
    'gold' : 19.32, # g/cm^3
    'iron' : 7.874,
    'titanium' : 4.5,
    '2gold' : 19.32,
    '3gold' : 19.32,
}
mass_map = {
    'gold' : 196.96657, # g/mol
    'iron' : 55.845,
    #'titanium' : ,
    '2gold' : 196.96657,
    '3gold' : 196.96657,
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
e_empty = ((5.486 * 0.86 + 5.443 * 0.127 + 5.391 * 0.014) / (0.86 + 0.127 + 0.014))
c = 299792458 # m/s
m_a = 3.7273794066e3 #MeV
m_e = 0.51099895000 #MeV
K = 0.307075 # MeV / mol * cm^2
v = lambda b : b * c
g = lambda b : 1 / math.sqrt(1-b**2)
W_max = lambda b : (2*m_e*b**2*g(b)**2) / (1+2*g(b)*m_e/m_a + (m_e/m_a)**2)
dEdx = lambda beta, foil : K * 2**2 * element_map[foil] / mass_map[foil] / beta**2 * (1/2 * math.log(2*m_e*beta**2*g(beta)**2*W_max(beta)/I[foil]**2) - beta**2) * density_map[foil]
b = lambda E : math.sqrt(1 - m_a**2 / E**2)

def thiccness(foil, incident, exiting):
    #WATCH UNITS OF X - cm
    E = incident + m_a
    x = 0
    dt = 1e-14
    print (b(E))
    while E > exiting + m_a:
        beta = b(E)
        x += v(beta) * 100 * dt # cm
        #E -= dEdx(beta, foil) * v(beta) * dt
        dE = dEdx(beta, foil) * v(beta) * dt
        E -= dE
        #print (dE)

    return x


def get_energies(foil):
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium"], condition = lambda metadata : metadata[0] == 'empty' or metadata[0] == foil)
    ref_data = iterationless(data)[('empty', 0)]
    ref_chs = [resultify(max_model(r[1])) for r in ref_data]
    ref_ch = resultify(weighted_average(ref_chs)) 
    foil_data = iterationless(data)[(foil, 0)]   
    foil_chs = [resultify(max_model(f[1])) for f in foil_data]
    foil_ch = resultify(weighted_average(foil_chs))
    e = lambda ch : Result(e_empty) / ref_ch * ch
    incident = e(ref_ch)
    exiting = e(foil_ch)
    return incident, exiting

def get_thickness(foil, do_report = False):
    incident, exiting = get_energies(foil)

    #thicc = thiccness(foil, incident.val, exiting.val)
    #thicc_high_stat = thiccness(foil, incident.val + incident.stat, exiting.val - exiting.stat)
    #thicc_low_stat = thiccness(foil, incident.val + incident.stat, exiting.val - exiting.stat)
    #thicc_high_sys = thiccness(foil, incident.val + incident.sys, exiting.val - exiting.sys)
    #thicc_low_sys = thiccness(foil, incident.val + incident.sys, exiting.val - exiting.sys)
    thicc = AsymmetricResult.asymmetric_evaluate(thiccness, foil, incident, exiting)

    if do_report:
        report("Thickness of "+foil, thicc, "cm")

    return thicc



def get_scattering_energy(foil, mode = "naive"):
    if mode == "naive":
        info = get_energies(foil)
        return Result(info[0].val, (info[0] - info[1]).val)
    else:
        #TODO
        pass



if __name__ == '__main__':
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium"])
    #print (data.keys())
    #print (data.keys())
    #print (data.keys())

    print ("Expected beam energy: "+str(e_empty)+" MeV")

    ref_data = iterationless(data)[('empty', 0)]

    empty_mean_vline = ("Estimated empty mean", max_model(ref_data[0][1])[0], 'red') # for use in plots.py
    
    plot_histogram(('empty', 0), ref_data[0][1], vlines = [empty_mean_vline]) # ref_data is a list of tuples for all empty scans. Tuples are time, hist
    
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
            entry_vline = ("Estimated "+metadata[0]+" mean", max_model(entry[1])[0], 'orange')
            plot_histogram(metadata, entry[1], [empty_mean_vline, entry_vline])
            energies[metadata[0]] = e(resultify(max_model(entry[1])))
            report (metadata[0], energies[metadata[0]], "MeV")
            #print (str(metadata[0])+": "+str(energies[metadata[0]][0])+" MeV ± "+str(energies[metadata[0]][1]))

    gold_tracking = []
    #Z alpha is 2

    for foil, energy in energies.items():
        #print (ref_e)
        #print (energy)
        energy_loss = ref_e - energy #sum(ref_e, minus(energy))
        
        thickness = energy_loss / (K * 2**2 * element_map[foil] / mass_map[foil])
        #e(ref_ch) - energies[foil]
        if 'gold' in foil: 
            gold_tracking.append(energy_loss / number[foil])
        report (str(foil) + " energy loss", energy_loss, "MeV")
        #print (str(foil)+" energy loss: "+str(energy_loss[0])+" MeV ± "+str(energy_loss[1])+" MeV")
    #print ([entry.val for entry in gold_tracking])
    report ("Energy loss per gold foil: ", weighted_average(gold_tracking), "MeV")


