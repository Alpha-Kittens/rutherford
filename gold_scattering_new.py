import os
from data_loader import read_data, recursive_read
from data_processing import *
from models import *
import matplotlib.pyplot as plt
import numpy as np
import math
from energy_loss import element_map


def scattering_unpack(cps, angle, entry, angle_err = 0.5):
    time, histogram = entry
    cps.append(Result(np.sum(histogram) / time, stat = np.sqrt(np.sum(histogram) / time)))
    angle.append(Result(angle, stat = angle_err))

def plotting_unpack(results, mode = 'tot'):
    x = []
    xerr = []
    for result in results:
        x.append(result.val)
        if mode == 'tot':
            xerr.append(result.tot)
        elif mode == 'stat':
            xerr.append(result.stat)
        elif mode == 'sys':
            xerr.append(result.sys)
    return x, xerr


def do_scattering_thing(element, min_angle, folder, emoji = ":P"):

    data = {}
    data = recursive_read(data, folder, require = [element], condition = lambda metadata : metadata[1] > min_angle)

    cps_raw = []
    angle_raw = []
    for entry in data:
        scattering_unpack(cps_raw, angle_raw, entry)
    
    angle, angle_err = plotting_unpack(angle_raw)
    cps, cps_err = plotting_unpack(cps_raw)

    plt.errorbar(angle, cps, xerr = angle_err, yerr = cps_err, marker = '.', ls = 'none')
    plt.xlabel("Angle (degrees)")
    plt.ylabel("CPS")
    plt.title(element + " scattering " + emoji)
    plt.show()

    angles = np.linspace(min_angle, 180, 1000)
    
    
    x = np.array(angle)
    Z = element_map(element)
    E = "skadoosh" #todo
    y = np.array(cps) * (E / Z) **2





    



files = os.listdir(folder)

angles = []
cpss = []
errors = []
for file_name in files:
    fp = folder + file_name

    information = read_data(fp)

    cps = information['cps']
    angle = information['angle']
    time = information['time']

    error = math.sqrt(cps)/math.sqrt(time)

    if(angle >= 10):
        angles.append(angle)
        cpss.append(cps)
        errors.append(error)



#plt.scatter(angles, cps, marker='o')
plt.errorbar(angles, cpss, yerr = np.array(errors) ,xerr = 0.5, marker='o', ls = 'none')
plt.xlabel('angle')
plt.ylabel('cps')
#plt.yscale('log')
plt.title('gold scattering :P')
plt.show()


if __name__ == '__main__':
    element = 'gold'
    min_angle = 10
    folder = 'gold_scattering/'
    do_scattering_thing(element, min_angle, folder)