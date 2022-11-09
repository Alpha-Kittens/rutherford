import os
from data_loader import read_data
import matplotlib.pyplot as plt
import numpy as np
import math
from gold_scattering_new import scattering


folder = 'gold_scattering/'

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


scattering('gold', 10, folder, emoji = ":P")