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
        angles.append(angle)
        cpss.append(cps)
        errors.append(error)

plt.errorbar(angles, cpss, yerr = errors, xerr = 0.5, marker='o', ls='none')
plt.xlabel('angle')
plt.ylabel('cps')
plt.title('beam profile :O')
plt.xticks(range(-10,11, 2))
plt.show()


# Use approximations from linear fit to determine true y errors
for i in range(len(errors)):
    errors[i] = math.sqrt((errors[i])**2 + (0.5 * 30)**2) # slope is about 30 and the x error is 0.5


beam_profile_fit(angles, cpss, errors, choiceL = 'linear', choiceR= 'linear')