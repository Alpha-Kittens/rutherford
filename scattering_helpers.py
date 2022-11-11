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

def scattering_model(x, a = 1e-5, b = 0):
    try:
        return a * f(x) + b
    except:
        y = []
        for i in x:
            y.append(scattering_model(i))
        return y


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
        print(result.val)
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





def convolve_with_rutherford(profile, min_angle = 10):
    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 180, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                             # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)
    f = lambda theta : 1/(np.sin(theta / 2 * np.pi / 180))**4    

    scattering,pdomain = convolve(profile, bpdomain, f, angles)

    return scattering, pdomain, angles


def get_rutherford_convolutions(profile_sets, min_angle = 10, plot=False):
    convolutions = []
    domains = []
    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 180, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                             # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)

    for i in range(len(profile_sets)):
        profile = profile_sets[i]
        scattering,pdomain = convolve(profile, bpdomain, f, angles)
        convolutions.append(scattering)
        domains.append(pdomain)
        if plot:
            plt.plot(pdomain, (1e-4) * scattering, label = "convolution " + str(i))


    if plot:
        plt.plot(pdomain, f(pdomain), label = "scattering expectation", color = "red")
        plt.legend()
        plt.show()

    print(type(convolutions))
    return convolutions, domains


def get_scattering_data(element, min_angle, folder, plot=True, emoji = ":P"):
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

    if plot:
        plt.errorbar(angle, cps, xerr = angle_err, yerr = cps_err, marker = '.', ls = 'none', color = 'black')
        plt.xlabel("Angle (degrees)")
        plt.ylabel("CPS")
        plt.title(element + " scattering " + emoji)
        plt.show()
    return angle, cps, angle_err, cps_err


def get_rutherford_models(rutherford_convolutions, domains, plot=False):
    models = []
    for i in range(len(rutherford_convolutions)):
        convolution = rutherford_convolutions[i]
        f_domain = domains[i]
        func, domain = interpolate(convolution, f_domain)
        model = lambda x, a=(1e-10), b=0 : a * np.array(func(x)) + b
        models.append(model)
       
        
        if plot:
            x_vals=np.linspace(min(f_domain) + 0.01, max(f_domain) - 0.01, 1000)
            y_vals = []
            for x in x_vals:
                y_val = (1e-4) * model(x)
                y_vals.append(y_val)
            plt.plot(x_vals, y_vals, label = "convolution " + str(i))
            plt.show()
    return models


def fit_to_scattering(element, min_angle, folder, model, plot=True):
    angle, cps, angle_err, cps_err = get_scattering_data(element, min_angle, folder, plot=False)

    #themodel = lmfit.Model(model)

    weights = []
    for i in cps_err:
        weights.append(1/i)

    #themodel = lmfit.Model(scattering_model)
    themodel = lmfit.Model(model)
    result  = themodel.fit(cps, x=angle, weights=weights)

    print('a: ' +  str(result.params['a'].value))
    print('b: ' +  str(result.params['b'].value))
    print('red chi: ' + str(result.redchi))
    

    if plot:
        print('plotting now')
        plt.scatter(angle, cps, label='data')
        params = result.params
        init_params = result.init_params
        x_vals=np.linspace(min(angle), max(angle), 1000)
        y_eval=[]
        for i in x_vals:
            #y_eval.append(model(i, params['a'].value, params['b'].value))
            y_eval.append(model(i, params['a'].value, params['b'].value))
        plt.plot(x_vals, y_eval, label='fit')

        initial_y= []
        for i in x_vals:
            #initial_y.append(model(i, init_params['a'].value, params['b'].value))
            initial_y.append(model(i, init_params['a'].value, params['b'].value))
        plt.plot(x_vals, initial_y, label = 'initial guess')
        plt.legend()
        plt.show()
    
    return result

def rutherford_fits(element, min_angle, folder, convolutions, domains, plot=True):
    chi2 = []
    for i in range(len(convolutions)):
        conv = convolutions[i]
        domain = domains[i]
        function, domain = interpolate(conv, domain)
        print(':P')
        model = lambda x, a=(1e-10), b=0 : a * np.array(function(x)) + b
        result = fit_to_scattering(element, min_angle, folder, model, plot=False)
        chi = result.redchi
        chi2.append(chi)
    print(chi2)
    plt.hist(chi2)
    #plt.plot(chi2)
    plt.show()