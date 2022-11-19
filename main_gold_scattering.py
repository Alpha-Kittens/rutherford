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

triangle_fit(x=angles, y=cpss, report=True, show=True, initial_plot=False)
data_sets = generate_data_sets(total=10)
profile_sets = fit_data_sets(data_sets, show=True)


#profile_sets = [profile_sets[0]]

###
'''
RUTHERFORD SCATTERING
'''
###

convolutions1, domains = get_convolutions1(profile_sets, choice = 'rutherford', min_angle = 1, plot=True)
#convolutions2 = get_convolutions2(profile_sets, choice='rutherford', min_angle = 10, plot=True)
plumconvolutions,plumdomains = get_convolutions1(profile_sets, choice = 'plum pudding', min_angle = -20, plot=True)
#plumconvolutions2 = get_convolutions2(profile_sets, choice='plum pudding', min_angle = 0, plot=True)


data = get_scattering_data('gold', min_angle=10, folder = 'gold_scattering/', plot=True, viewHist=True)


processed_data = process_scattering_data(profile=profile_sets[0], data=data, plot=True)

plot_models(processed_data, convolutions1[0], plumconvolutions[0], domain_r=domains[0], domain_p=plumdomains[0])

compare_models_plot(processed_data, rutherford_convolution=convolutions1[0], domain=domains[0])
compare_chi2(processed_data, rutherford_convolutions=convolutions1, plum_convolutions = plumconvolutions, rutherford_domains=domains, plum_domains=plumdomains)
#compare_chi2(processed_data, processed_data_plum, rutherford_convolutions=convolutions2, plum_convolutions = plumconvolutions2)
