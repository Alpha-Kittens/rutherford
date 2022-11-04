import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np
from beam_profile_models import beam_profile_fit

folder = 'beam_profile/'

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

    if cps != 0:
        error = math.sqrt(cps)/math.sqrt(time)

        print('angle ' + str(angle) + ": " + str(error))
        
        angles.append(angle)
        cpss.append(cps)
        errors.append(error)

plt.errorbar(angles, cpss, yerr = errors, xerr = 0.5, marker='o', ls='none')
plt.xlabel('angle')
plt.ylabel('cps')
plt.title('beam profile :O')
plt.xticks(range(-10,11, 2))
plt.show()


beam_profile_fit(angles, cpss, errors)