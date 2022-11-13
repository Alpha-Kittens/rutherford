import os
from data_loader import read_data, recursive_read
from data_processing import *
from models import *
import matplotlib.pyplot as plt
import numpy as np
import math
from energy_loss import element_map
from beam_profile_models import beam_profile_fit
from profile import profile as beam_profile
import lmfit


f = lambda theta : 1/(np.sin(theta / 2 * np.pi / 180))**4    


def scattering_unpack(cps, angle, metadata, entry, time_err = 2/np.sqrt(12), angle_err = 0.5):
    """
    Given `metadata` and `entry` associated with a datapoint, appends corresponding counts per second value and angle value 
    as `Result` objects to given `cps` and `angle` arrays. 
    Arguments:
        * `cps` (array): array to which CPS `Result`s should be appended
        * `angle` (array): array to which Angle `Result`s should be appended
        * `metadata` (3-tuple): foil, angle, iteration. Same format as returned by `get_metadata`.
        * `entry` (2-tuple): time, histogram. Same format as added by `add_data` and `recursive_read`. 
        * `time_err` (float): error associated with all values of time. Since truncated, the most simplistic error is 2/sqrt(12) (default)
            since time is uniformly distributed between `time` and `time+1`. Most accurate is an asymmetric error bound, but `Result` does not 
            implement this yet. This might be my next task. If I do so, plotting_unpack results will still be compatible with plt.errorbar, so no need to change anything. (hopefully)
    Returns:
        * Nothing. Rather, appends to `cps` and `angle` arrays. 
    """
    time, histogram = entry
    cps.append(Result(np.sum(histogram), stat = np.sqrt(np.sum(histogram))) / Result(time, stat = time_err))
    #cps.append(Result(np.sum(histogram) / time, stat = np.sqrt(np.sum(histogram)) / time)) # outdated, before time_err implemented. 
    angle.append(Result(metadata[1], stat = angle_err))

def plotting_unpack(results, mode = 'tot'):
    """
    Given an array of `Result` objects, returns array of values and errors usable with `plt.errorbar`. 
    Arguments:
        * `results` (array): array of `Result` objects.
        * `mode` (string): specifying which error to use. 
            - `tot`: total error (default)
            - `stat`: statistical error
            - `sys`: systematic error
    Returns:
        * `x`: array of values
        * `xerr`: array of errors associated with `x`
    """
    x = []
    xerr = []
    for result in results:
        #print(result.val)
        x.append(result.val)
        if mode == 'tot':
            xerr.append(result.tot)
        elif mode == 'stat':
            xerr.append(result.stat)
        elif mode == 'sys':
            xerr.append(result.sys)
    return x, xerr


def scattering(element, min_angle, folder, emoji = ":P"):
    """
    Compares 1/sin^4(theta/2) angular dependence with Plum Pudding model. 
    Arguments:
        * `element`: foil to use. 
        * `min_angle`: minimum angle to use for data
        * `folder`: folder containing data
        * `emoji`: emoji to use for visualization purposes. Default: `":P"`

    Returns:
        * TBD, maybe nothing

    """

    """
    Notes: 
    -In order to ensure domain of valid convolution extends through angles as small as 10 degrees, 
    it is likely necessary to implement the "cutoff" towards the origin as recommended by Prof. Harris. 
    -I've decided that I don't need this for Z-dependence, since multiplication by a constant is an operation commutable with convolution
    -However, I will be using scattering_unpack and plotting_unpack. 
    -This will rely on energy_loss.py in order to get the scattering energy and uncertainties. I will impelemnt the following method:
    -get_energy(foil, mode) -> returns a Result object with the energy of scattering and proper statistical error. 
    modes: 
        -"naive": does it the "naive" way as discussed, i.e. energy is beam energy and uncertainty is total energy loss. Will implement asymmetric bounds soon.
        -"geometric": does it assuming geometric distribution of scattering locations inside foil. Again, will have asymmetric error bounds. 


    """

    data = {}
    recursive_read(data, folder, require = [element], condition = lambda metadata : metadata[1] > min_angle)

    cps_raw = []
    angle_raw = []
    for metadata, entry in data.items():
        scattering_unpack(cps_raw, angle_raw, metadata, entry)

    print (cps_raw)
    print (angle_raw)
    
    angle, angle_err = plotting_unpack(angle_raw)
    cps, cps_err = plotting_unpack(cps_raw)

    print ("--")
    print (angle)
    print (angle_err)
    print (cps)
    print (cps_err)
    plt.errorbar(angle, cps, xerr = angle_err, yerr = cps_err, marker = '.', ls = 'none', color = 'black')
    plt.xlabel("Angle (degrees)")
    plt.ylabel("CPS")
    plt.title(element + " scattering " + emoji)
    plt.show()

    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 180, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                             # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)
    f = lambda theta : 1/(np.sin(theta / 2 * np.pi / 180))**4    

    scattering,pdomain = convolve(beam_profile, bpdomain, f, angles)

    #print (pdomain)
    #print (scattering)

    plt.plot(pdomain, (1e-4) * scattering, label = "convolution", color = 'purple')
    plt.plot(bpdomain, beam_profile(bpdomain), label = "beam profile", color = "blue")
    plt.plot(angles, f(angles), label = "scattering expectation", color = "red")
    plt.legend()
    plt.xlabel("Angle (degrees)")
    plt.ylabel("CPS")
    plt.show()
    #model = lambda x, c : c * scattering(x)
    #result = lmfit.fit()


if __name__ == '__main__':
    element = 'gold'
    min_angle = 8
    folder = 'gold_scattering/'
    scattering(element, min_angle, folder)



'''
'''
'''
DATA METHODS
'''
'''
'''
def get_scattering_data(element, min_angle, folder, plot=True, emoji = ":P", show=True):
    '''
    gets the scattering data from a selected folder, and a given element
    '''
    data = {}
    recursive_read(data, folder, require = [element], condition = lambda metadata : metadata[1] > min_angle)

    cps_raw = []
    angle_raw = []
    for metadata, entry in data.items():
        scattering_unpack(cps_raw, angle_raw, metadata, entry)
    
    angle, angle_err = plotting_unpack(angle_raw)
    cps, cps_err = plotting_unpack(cps_raw)

    data = [angle, cps, angle_err, cps_err]

    if plot:
        plot_data(data, show=show, title = element + " scattering " + emoji)

    return data


def plot_data(data, show=True, title = None, ylabel='CPS'):
    '''
    plots any kind of data data must of be the form list: [x, y, xerr, yerr]
    ylabel is default set to CPS
    Has options to not show the plot, to add a title, or to modify the y-axis label
    '''
    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    plt.errorbar(x, y, xerr = xerr, yerr = yerr, marker = '.', ls = 'none', color = 'black')
    plt.xlabel("Angle (degrees)")
    plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    if show:
        plt.show()




def process_scattering_data(profile, data, plot=False):
    '''
    Given the result of "get scattering data", it propogates errors like the x-errors and the energy errors to give data
    ready to be fit.
    '''
    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    from energy_loss_2 import expected_E_square

    # Energy Uncertainties
    E = expected_E_square('gold')
    Evalue = E.val
    Eerror = math.sqrt(E.sys **2 + E.stat **2)
    newy = Evalue * np.array(y)

    newyerr = []
    for i in range(len(yerr)):
        new_err = math.sqrt(yerr[i]**2 + (y[i]*Eerror)**2)
        newyerr.append(new_err)

    data_new = [x, newy, xerr, newyerr]

    if plot:
        plot_data(data_new, show=True, title="New data (energy errors added)", ylabel='CPS * E^2')

    # Propogation of x-errors
    slope_function = approximate_conv_slopes(profile, data = data_new, plotConv=True)
    if plot:
        x_vals=np.linspace(min(x), max(x), 1000)
        y_eval=[]
        for i in x_vals:
            y_eval.append(slope_function(i))
        plt.plot(x_vals, y_eval, label='slope')
        plt.legend()
        plt.show()
    
    final_yerrs = []
    final_x_errs = []
    for i in range(len(newyerr)):
        final_x_errs.append(0)
        final_yerr = math.sqrt(newyerr[i]**2 + (xerr[i] * slope_function(x[i]))**2)
        final_yerrs.append(final_yerr)
    
    final_data = [x, newy, final_x_errs, final_yerrs]
    if plot:
        plot_data(final_data, show=True, title="Processed Data (All errors propogated)", ylabel='CPS * E^2')

    return final_data



def approximate_conv_slopes(profile, data, min_angle=10, plotConv=False):
    '''
    Using a basic profile fit, convolves with rutherford scattering to get a slope estimation at any point
    Note that this uses convolution1 not convolution2
    '''
    convolution, domain, angles = convolve_with_rutherford1(profile, min_angle)

    
    function, domain = interpolate(convolution, domain)
    model = lambda x, a=(1e-10), b=0 : a * np.array(function(x)) + b
    result = fit_to_scattering(data, model, plot=False)

    a = result.params['a'].value
    b = result.params['b'].value

    if plotConv:
        plot_fit(result, model, data, show=False, label='Convolution')

    return lambda x: convolution_slope(x, a * np.array(convolution), domain)







'''
'''
'''
CONVOLUTION METHODS
'''
'''
'''
def convolve_with_rutherford1(profile, min_angle = 10):
    '''
    implements the convolve method in models to convolve a beam profile with the rutherford model
    returns the resulting array (not a callable), the domain in which it is valid, and the angles 
    '''
    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 180, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                             # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)
    f = lambda theta : 1/(np.sin(theta / 2 * np.pi / 180))**4    

    scattering,pdomain = convolve(profile, bpdomain, f, angles)

    return scattering, pdomain, angles


def get_rutherford_convolutions1(profile_sets, min_angle = 10, plot=False):
    '''
    given a bunch of profile sets, implements convolve_with_rutherford for each one and returns a list of the results
    '''
    convolutions = []
    domains = []
    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 180, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                             # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)

    for i in range(len(profile_sets)):
        profile = profile_sets[i]
        scattering,pdomain = convolve(profile, bpdomain, f, angles)
        normalization = f(min(pdomain)) / scattering[0]
        convolutions.append(normalization * np.array(scattering))
        domains.append(pdomain)
        if plot:
            plt.plot(pdomain, normalization * scattering, label = "convolution " + str(i))


    if plot:
        plt.plot(pdomain, f(pdomain), label = "scattering expectation", color = "red")
        plt.legend()
        plt.show()

    return convolutions, domains


def convolve_with_rutherford2(profile, min_angle = 10):
    '''
    implements the second version of convolve, returns a callable function
    '''
    normalization = f(min_angle) / convolve2(f, profile, min_angle)
    return lambda x : normalization * convolve2(f, profile, x)

def get_rutherford_convolutions2(profiles, min_angle = 10, plot=True):
    '''
    implements convolve_with_rutherford2 for each profile and returns an array of convolutions (callable) given a list of profiles
    '''
    convolutions = []
    for profile in profiles:
        convolutions.append(convolve_with_rutherford2(profile, min_angle))    
    if plot:
        pdomain = np.linspace(min_angle, 60, 1000)

        for i in range(len(convolutions)):
            conv = convolutions[i]
            plt.plot(pdomain, conv(pdomain), label = "convolution " + str(i))

        plt.title("Convolutions of beam profile with Rutherford Scattering Expectation 2")
        plt.plot(pdomain, f(pdomain), label = "scattering expectation", color = "red")
        plt.legend()
        plt.show()
    return convolutions



'''
'''
'''
FITTING METHODS
'''
'''
'''
def scattering_model(x, a = 1e-5, b = 0):
    '''
    the model for pure (no convolution) rutherford scattering
    '''
    try:
        return a * f(x) + b
    except:
        y = []
        for i in x:
            y.append(scattering_model(i))
        return y

def fit_to_scattering(data, model, plot=True):
    ''''
    fits a given model to the given data. Again data must be of the form list: [x, y, xerr, yerr]
    '''
    #angle, cps, angle_err, cps_err = get_scattering_data(element, min_angle, folder, plot=False)

    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]


    weights = []
    for i in yerr:
        weights.append(1/i)

    themodel = lmfit.Model(model)
    result  = themodel.fit(y, x=x, weights=weights)

    print('red chi: ' + str(result.redchi))
    
    if plot:
        plot_fit(result, model, element, min_angle, folder, show=True)

    return result

def rutherford_fits(data, convolutions, domains = None, plot=True):
    '''
    given a bunch of convolutions (from either convolve method convolve1 or convolve2)
    Note that if using convolutions from convolve1 you must also pass domains
    it fits all of them and returns the distribution of chi2
    '''
    chi2 = []
    for i in range(len(convolutions)):
        print("Fitting Function: " + str(i))
        try: 
            func = convolutions[i]
            model = lambda x, a=(1e-10), b=0 : a * np.array(func(x)) + b
            result = fit_to_scattering(data, model, plot=False)

        except:
            conv = convolutions[i]
            domain = domains[i]
            function, domain = interpolate(conv, domain)
            print(':P')
        
            model = lambda x, a=(1e-10), b=0 : a * np.array(function(x)) + b

            result = fit_to_scattering(data, model, plot=False)
        
        chi = result.redchi
        chi2.append(chi)
    if plot:
        plt.hist(chi2)
        plt.xlabel('reduced chi square')
        plt.ylabel('frequency')
        plt.show()
    
    return chi2


def plot_fit(result, model, data, show=True, label=None, initial=False):
    '''
    given a result and model for scattering data, it plots the fit.
    label lets you change what is given in the legend for this fit
    intial is a boolean letting you decide if you want to plot intial fit guesses
    '''
    plot_data(data, show=False, title = None)

    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    params = result.params
    init_params = result.init_params
    x_vals=np.linspace(min(x), max(x), 1000)
    y_eval=[]
    for i in x_vals:
        y_eval.append(model(i, params['a'].value, params['b'].value))
    if label is not None:
        plt.plot(x_vals, y_eval, label='fit ' + label)
    else:
        plt.plot(x_vals, y_eval, label='fit')

    if initial:
        initial_y= []
        for i in x_vals:
            initial_y.append(model(i, init_params['a'].value, params['b'].value))
        if label is not None:
            plt.plot(x_vals, initial_y, label = 'initial guess ' + label)
        else:
            plt.plot(x_vals, initial_y, label = 'initial guess')
    if show:
        plt.legend()
        plt.show()


def rutherford_scattering_fit(data, plot=False):
    '''
    pure rutherford (no convolution) fit to  data
    again, data must be of the form list: [x, y, xerr, yerr]
    '''
    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    themodel = lmfit.Model(scattering_model)
    weights = []
    for i in yerr:
        weights.append(1/i)
    result2  = themodel.fit(y, x=x, weights=weights)

    if plot:
        plot_fit(result2, scattering_model, element, min_angle, folder, show=True, label='rutherford scattering')
    
    return result2



'''
RESULTS METHODS
'''
def compare_models_plot(data, rutherford_convolution, domain=None):
    '''
    Takes a rutherford convolution, and data specifications and prints the fitting for that in comparison to no convolution
    '''
    try:
        function, domain = interpolate(rutherford_convolution, domain)
        print(':P')
        model = lambda x, a=(1e-10), b=0 : a * np.array(function(x)) + b
        result = fit_to_scattering(data, model, plot=False)
    except:
        model = lambda x, a=(1e-10), b=0 : a * np.array(rutherford_convolution(x)) + b
        result = fit_to_scattering(data, model, plot=False)
    
    print("Fit complete. Plotting Now")
    plot_fit(result, model, data, show=False, label='beam profile x rutherford')

    result2 = rutherford_scattering_fit(data, plot=False)

    plot_fit(result2, scattering_model, data, show=False, label='rutherford scattering')

    plt.legend()
    plt.show()


def compare_chi2(data, rutherford_convolutions, domains=None):
    '''
    Plots the distribution of chi2 for rutherford convolutions and a vertical line for 
    no convolution fit chi2
    '''
    no_conv_chi = rutherford_scattering_fit(data).redchi
    print(no_conv_chi)

    plt.vlines(no_conv_chi, ymin=0, ymax=5, label='no convolution', linestyles='dashed', color='red')

    rutherford_chi = rutherford_fits(data, rutherford_convolutions, domains, plot=False)

    plt.hist(rutherford_chi, label='convolutions')
    plt.xlabel('reduced chi square')
    plt.ylabel('frequency')
    plt.legend()
    plt.show()


