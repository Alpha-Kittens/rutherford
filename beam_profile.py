import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np
from beam_profile_models import beam_profile_fit
from beam_profile_models import profile_sys_error
plt.style.use('seaborn-colorblind')


def get_fit():
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

<<<<<<< Updated upstream
    plt.errorbar(angles, cpss, yerr = np.array(errors)*50, xerr = 0.5, marker='o', ls='none')
    plt.xlabel('angle')
    plt.ylabel('cps')
    plt.title('beam profile :O')
    plt.xticks(range(-10,11, 2))
    plt.show()


    # Approximate y errors from x errors with slope from linear fit. For points outside of +/- 5, the slope is 0
    '''
    for i in range(len(errors)):
        if abs(angles[i]) < 5:
            errors[i] = math.sqrt((errors[i])**2 + (0.5 * 30)**2) # slope is about 30 and the x error is 0.5
    '''
=======
print(cpss)
plt.errorbar(angles, np.array(cpss), yerr = np.array(errors), xerr = 0.5, marker='o', ls='none')
plt.xlabel('angle')
plt.ylabel('cps')
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

'''
for i in range(len(errors)):
    if abs(angles[i]) < 5:
        errors[i] = math.sqrt((errors[i])**2 + (0.5 * 30)**2) # slope is about 30 and the x error is 0.5
'''
>>>>>>> Stashed changes


model_comparison = {
    'linear - linear': 0,
    'exponential - linear' : 0,
    'quadratic - linear': 0,
}


model_comparison['linear - linear'] = beam_profile_fit(angles, cpss, errors, choiceL = 'linear', choiceR= 'linear', plot=True).redchi

params_L,params_R = profile_sys_error(angles, cpss, errors, angle_error = 0.5, choiceL = 'linear', choiceR= 'linear')


print(model_comparison)

<<<<<<< Updated upstream
    return beam_profile_fit(angles, cpss, errors, choiceL = 'linear', choiceR= 'linear')

if __name__ == '__main__':
    get_fit()
=======
>>>>>>> Stashed changes
