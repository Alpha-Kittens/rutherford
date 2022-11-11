import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np
from beam_profile_models import beam_profile_fit
from beam_profile_models import profile_sys_error
from profile import profile
plt.style.use('seaborn-colorblind')


folder = 'beam_profile/'

files = os.listdir(folder)

angles = []
cpss = []
errors = []
counts = []
for file_name in files:
    fp = folder + file_name

    information = read_data(fp)

    
    cps = information['cps']
    count = information['counts']
    angle = information['angle']
    time = information['time']

    if cps != 0:
        error = math.sqrt(cps)/math.sqrt(time)
        angles.append(angle)
        cpss.append(cps)
        errors.append(error)
        counts.append(count)

plt.errorbar(angles, cpss, yerr = np.array(errors), xerr = 0.5, marker='o', ls='none')
plt.xlabel('angle')
plt.ylabel('counts')
plt.title('beam profile :O')
plt.xticks(range(-10,11, 2))
plt.show()



#Plot fractional uncertainty
frac_unc = []
for i in range(len(errors)):
    frac_unc.append(errors[i]/cpss[i])

plt.errorbar(angles, np.array(frac_unc), marker='o', ls='none')
plt.xlabel('angle')
plt.ylabel('fractional uncertainy')
plt.title('beam profile - uncertainties :O')
plt.xticks(range(-10,11, 2))
plt.show()


# Approximate y errors from x errors with slope from linear fit. For points outside of +/- 5, the slope is 0
for i in range(len(errors)):
    if abs(angles[i]) < 7.5:
        errors[i] = math.sqrt((errors[i])**2 + (0.5 * 30)**2) # slope is about 30 and the x error is 0.5



model_comparison = {
    'linear - linear': 0,
    'exponential - linear' : 0,
    'quadratic - linear': 0,
}


choiceL = 'linear'
choiceR = 'linear'
result = beam_profile_fit(angles, cpss, errors, choiceL = choiceL, choiceR= choiceR, plot=True)


def cwrite(file, params, choiceL, choiceR):
    with open(file, 'w') as f:
        print ("now writing")
        f.write("from beam_profile_models import *\n")
        f.write("from numpy import sqrt, abs\n")
        f.write('import lmfit \n')
        f.write("params = lmfit.Parameters() \n")
        for key in params:
            f.write('params.add( \'' + str(key) + '\', ' + 'value=' + str(params[key]) + ') \n')
        f.write("def profile (x): \n")
        #f.write('\t' + 'return evaluate_beam_model(x, \'' + str(choiceL) + '\', \'' + str(choiceR) + '\', ' + 'params)')
        f.write('\t' + 'return use_beam_model(x, \'' + str(choiceL) + '\', \'' + str(choiceR) + '\', ' + 'params)')
        f.flush()


toWrite = input("Write to file: (Y/N)")

if toWrite == 'Y':
    cwrite('profile.py', result.best_values, choiceL, choiceR)

