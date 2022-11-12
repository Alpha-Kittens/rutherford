import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np
from beam_profile_models import beam_profile_fit
from beam_profile_models import profile_sys_error
from beam_profile_helpers import*
from scattering_helpers import*
plt.style.use('seaborn-colorblind')



###
'''
BEAM PROFILE
'''
###
view_omitted_data()
plot_profile_data(show=True)
plot_frac_uncertainties()

angles, cpss, errors = get_data()

triangle_fit(x=angles, y=cpss, report=True, show=True, initial_plot=True)
data_sets = generate_data_sets()
profile_sets = fit_data_sets(data_sets, show=True)


#profile_sets = [profile_sets[0]]

###
'''
RUTHERFORD SCATTERING
'''
###
convolutions, domains = get_rutherford_convolutions(profile_sets, min_angle = 1, plot=True)
'''
rutherford_models = get_rutherford_models(convolutions, domains)
#fit_to_scattering('gold', min_angle=10, folder = 'gold_scattering/', model=rutherford_models[0])
for model in rutherford_models:
    x_vals=np.linspace(11 + 0.01, 60 - 0.01, 1000)
    y_vals = []
    for x in x_vals:
        y_val = (1e-4) * model(x)
        y_vals.append(y_val)
    plt.plot(x_vals, y_vals, label = "convolution " )
plt.show()
'''




rutherford_fits('gold', min_angle=10, folder = 'gold_scattering/', convolutions=convolutions, domains=domains)