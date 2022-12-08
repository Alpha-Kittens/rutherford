from energy_loss_2 import *
from data_loader import *
from beam_profile_helpers import*
from scattering_helpers import*

angles, cpss, errors = get_data()

triangle_fit(x=angles, y=cpss, report=True, show=True, initial_plot=True)
data_sets = generate_data_sets(total=10)
profile_sets = fit_data_sets(data_sets, show=True)


#profile_sets = [profile_sets[0]]

###
'''
RUTHERFORD SCATTERING
'''
###

convolutions1, domains = get_rutherford_convolutions1(profile_sets, min_angle = 1, plot=True)
#convolutions2 = get_rutherford_convolutions2(profile_sets, min_angle = 10, plot=True)

data_gold = get_scattering_data('gold', min_angle=10, folder = 'gold_scattering/', plot=True)
data_iron = get_scattering_data('iron', min_angle=10, folder = 'iron_scattering/', plot=True)

print (data_iron)
print (data_gold)


processed_data_gold = process_scattering_data(profile=profile_sets[0], data=data_gold, plot=True)
processed_data_iron = process_scattering_data(profile=profile_sets[0], data=data_iron, plot=True)
#compare_models_plot(processed_data, rutherford_convolution=convolutions1[0])
compare_chi2(processed_data_gold, rutherford_convolutions=convolutions1, domains=domains)
compare_chi2(processed_data_iron, rutherford_convolutions=convolutions1, domains=domains)