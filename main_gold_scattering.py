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
#view_omitted_data()
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

#convolutions1, domains = get_rutherford_convolutions1(profile_sets, min_angle = 1, plot=True)
convolutions2 = get_rutherford_convolutions2(profile_sets, min_angle = 10, plot=True)

data = get_scattering_data('gold', min_angle=10, folder = 'gold_scattering/', plot=True)

#rutherford_fits('gold', min_angle=10, folder = 'gold_scattering/', convolutions=convolutions1, domains=domains)
#rutherford_fits('gold', min_angle=10, folder = 'gold_scattering/', convolutions=convolutions2)

#compare_models_plot(data, rutherford_convolution=convolutions2[0])
process_scattering_data(profile=profile_sets[0], data=data)
#compare_chi2(data, rutherford_convolutions=convolutions2)