import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np
from beam_profile_models import beam_profile_fit
from beam_profile_models import profile_sys_error
from plots import plot_histogram
from profile import profile
plt.style.use('seaborn-colorblind')


folder = 'beam_profile/'
files = os.listdir(folder)

view_data = True

angles = []
cpss = []
errors = []
counts = []
omitted_data = []
cutoffs = []

for file_name in files:
    fp = folder + file_name
    information = read_data(fp)

    # Extracting the information
    cps = information['cps']
    count = information['counts']
    angle = information['angle']
    time = information['time']
    histogram = information['histogram']

    # Viewing all data (if desired)
    '''
    if view_data == True:
        plot_histogram(['empty', angle], histogram)

    cutoff = float(input("What is the cutoff?"))

    total_omitted_counts = 0

    for i in range(len(histogram)):
        if i<cutoff:
            total_omitted_counts += histogram[i]
            histogram[i] = 0
    
    plot_histogram(['empty', angle], histogram)

    cps = cps - total_omitted_counts
    '''
    
    # Omitting data
    '''
    if cps > 0.5:
        error = math.sqrt(cps)/math.sqrt(time)
        angles.append(angle)
        cpss.append(cps)
        errors.append(error)
        counts.append(count)
    else:
        omitted_data.append(information)
    '''
    if cps > 0:
        error = math.sqrt(cps)/math.sqrt(time)
    else:
        error = 0
    angles.append(angle)
    cpss.append(cps)
    errors.append(error)
    counts.append(count)
    #cutoffs.append(cutoff)
print(cutoffs)
    

def plot_profile_data(x = angles, y = cpss, yerr = errors, xerr = 0.5, emoji=':O', show=False):
    '''
    Plots the profile data.Can be used for any toy data set as well.
    Default is to plot the actual beam profile data
    '''
    plt.errorbar(x, y, yerr = yerr, xerr = xerr, marker='o', ls='none')
    plt.xlabel('angle')
    plt.ylabel('CPS')
    plt.title('beam profile ' + emoji)
    plt.xticks(range(-10,11, 2))
    if show:
        plt.show()


def plot_frac_uncertainties(emoji=':O'):
    '''
    Plots the fractional (poisson) errors for the beam profile data
    '''
    frac_unc = []
    for i in range(len(errors)):
        if cpss[i] > 0:
            frac_unc.append(errors[i]/cpss[i])
        else:
            frac_unc.append(0)


    plt.errorbar(angles, np.array(frac_unc), marker='o', ls='none')
    plt.xlabel('angle')
    plt.ylabel('fractional uncertainy')
    plt.title('beam profile - uncertainties :O')
    plt.xticks(range(-10,11, 2))
    plt.show()



def write_profile_to_file(file, params, choiceL, choiceR):
    '''
    Writes the profile to file (this is obsolete now)
    '''
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


def view_omitted_data():
    '''
    Plots the histograms of all the ones that were chosen to not be included
    '''
    for entry in omitted_data:
        cps = entry['cps']
        angle = entry['angle']
        histogram = entry['histogram']
        plot_histogram(['empty', angle], histogram)


def get_data():
    '''
    returns the beam profile data
    '''
    return angles, cpss, errors


def triangle_model(x, aL = 30, aR = -30, x0 = 0, y0 = 150):
    ''''
    the triangle model used for fitting to the beam profile
    '''
    a_L = aL
    a_R = aR
    x_0 = x0
    y_0 = y0

    try:
        if x < x0:
            y  = aL *(x-x0) + y0
        else:
            y = aR * (x-x0) + y0
        
        if y > 0:
            return y
        else:
            return 0
    except:
        y = []
        for i in x:
            y.append(triangle_model(i, aL = a_L, aR = a_R, x0 = x_0, y0 = y_0))
        return y

def triangle_fit(x, y, report=False, show=False, initial_plot = False):
    '''
    implements the triangle model to fit to beam profile data
    returns the fitted function
    '''
    
    model = lmfit.Model(triangle_model)

    result = model.fit(y, x=x)

    params = result.params
    init_params = result.init_params

    areaL = abs(0.5 * ((params['y0']/params['aL']) + params['x0'])*params['y0'])
    areaR = abs(0.5 * ((params['y0']/params['aR']) + params['x0'])*params['y0'])

    total_area = areaL+ areaR

    
    if report:
        plt.scatter(x, y, label='data')
    
        x_vals=np.linspace(min(x), max(x), 1000)
        y_eval=[]

        for i in x_vals:
            y_eval.append(triangle_model(i, params['aL'].value, params['aR'].value, params['x0'].value, params['y0'].value))
        plt.plot(x_vals, y_eval, label='fit')


        if initial_plot:
            initial_y= []
            for i in x_vals:
                initial_y.append(triangle_model(i, init_params['aL'].value, init_params['aR'].value, init_params['x0'].value, init_params['y0'].value))
            plt.plot(x_vals, initial_y, label = 'initial guess')
        
        plt.ylabel('CPS')
        plt.xlabel('Angle (degrees)')

        

        if show:
            plt.legend()
            plt.show()


    print(total_area)
    return lambda x : (1/(total_area)) * np.array(triangle_model(x, params['aL'].value, params['aR'].value, params['x0'].value, params['y0'].value))


def generate_data_sets(total = 10, view=True):
    '''
    generates new data sets (total= number of data sets you want to create) of the beam profile by sampling from a Gaussian
    '''
    data_sets = []
    for i in range(total):
        new_angles = []
        # Sample the angle from a gaussian distribution
        for angle in angles:
            new_angle = np.random.normal(loc=angle, scale=0.5)
            new_angles.append(new_angle)
        
        data_sets.append(new_angles)
    
    if view:
        for data_set in data_sets:
                plot_profile_data(data_set, cpss, yerr=0, xerr=0)
        plt.show()

    return data_sets

def fit_data_sets(data_sets, show = False):
    '''
    given a bunch of generated data sets, fits them to the triangle models
    returns an array of all the fits.
    '''
    fits = []
    for data_set in data_sets:
        if show:
            fit = triangle_fit(data_set, cpss, report=True, show=False)
        else:
            fit = triangle_fit(data_set, cpss, report=False, show=False)
        fits.append(fit)
    if show:
        plt.show()

    return fits
