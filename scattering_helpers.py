import os
from data_loader import read_data, recursive_read
from data_processing import *
from models import *
import matplotlib.pyplot as plt
import numpy as np
import math
from energy_loss import element_map
from beam_profile_models import beam_profile_fit
#from profile import profile as beam_profile
from plots import plot_histogram, rebin_hist
import lmfit
from histogram_fitters import fit_histogram



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
        x.append(result.val)
        if mode == 'tot':
            xerr.append(result.tot)
        elif mode == 'stat':
            xerr.append(result.stat)
        elif mode == 'sys':
            xerr.append(result.sys)
    return x, xerr

'''
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
Important Functions
'''
f = lambda theta : 1/(np.sin(theta / 2 * np.pi / 180))**4    # Rutherford scattering cross section





'''
'''
'''
DATA METHODS
'''
'''
'''
def get_scattering_data(element, min_angle, folder, plot=True, emoji = ":P", show=True, viewHist = False, magnify_errors = 1):
    '''
    gets the scattering data from a selected folder, and a given element
    '''
    data = {}
    recursive_read(data, folder, require = [element], condition = lambda metadata : metadata[1] > min_angle)

    omitted = {
        20: 0.0, 
        25: 0.0, 
        30: 0.0, 
        40: 0.0, 
        45: 0.0, 
        50: 0.0}
    
    
    cps_raw = []
    cps_raw_new = []
    background_error = []
    angle_raw = []
    for metadata, entry in data.items():
        scattering_unpack(cps_raw, angle_raw, metadata, entry)
        time = entry[0]

        histogram = entry[1]
        new_bins, new_hist = rebin_hist(histogram)
        if viewHist:
            plot_histogram('gold ' + str(angle_raw[len(angle_raw) - 1]), histogram)
            plot_histogram('gold ' + str(angle_raw[len(angle_raw) - 1]), new_hist, xaxis = new_bins)


        crystal_ball = fit_histogram(new_hist, plot=True)

        print(crystal_ball.params)
        print('red chi: ' + str(crystal_ball.redchi))

        total_omitted_counts = 0
        cutoff = 750
        for i in range(len(histogram)):
            if i<cutoff:
                total_omitted_counts += histogram[i]
                histogram[i] = 0
        '''
        if viewHist:
            plot_histogram(['gold', angle_raw[len(angle_raw) - 1]], histogram)
        '''
        omitted[metadata[1]] = total_omitted_counts/time

        #cps_raw_new.append(cps_raw[len(cps_raw) - 1] - (omitted[metadata[1]]/2))
        background_error.append(omitted[metadata[1]]/2)


    print(omitted)
    
    angle, angle_err = plotting_unpack(angle_raw)
    cps, cps_err = plotting_unpack(cps_raw)
    cps_new,cps_new_err = plotting_unpack(cps_raw_new)

    total_error = []
    for i in range(len(background_error)):
        total_error.append(math.sqrt(background_error[i] **2 + cps_new_err[i]**2))

    print(total_error)

    data = [angle, cps, 0, 0]
    data_new = [angle, cps_new, 0, background_error]
    data_black = [angle, cps_new, total_error, 0]
    data_with_background = [angle, cps_new, angle_err, total_error]

    print("Angle: ")
    print(angle)
    print("CPS")
    print(cps_new)
    print("Background Error: ")
    print(background_error)
    print("Possion Error: ")
    print(cps_new_err)

    background_error_big = []
    total_error_big = []

    '''
    for i in background_error:
        background_error_big.append(i*5)
    for i in total_error:
        total_error_big.append(i*5)
    
    mag_background_data = [angle, cps_new, angle_err, total_error_big]
    mag_new_data = [angle, cps_new, 0, background_error_big]
    '''


    if plot:
        plt.text(25, 1e-3, 'Errors in CPS displayed 5 times \nlarger than true value')
        plot_data(data_with_background, show=False, mode='raw', color='b', label='total error')
        plot_data(data_new, show=show, mode='raw', color='orange', label='background error')
        #plot_data(data_black, show=show, title = element + " scattering " + emoji, mode='raw', color='black', label='new data')

        #plot_data(data, show=show, title = element + " scattering " + emoji, mode='raw', color='red', label='original data')


    return data_with_background

def plot_data(data, show=True, title = None, mode='processed', color='black', label = None):
    '''
    plots any kind of data data must of be the form list: [x, y, xerr, yerr]
    ylabel is default set to CPS
    Has options to not show the plot, to add a title, or to modify the y-axis label
    '''
    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    if label is not None:
        plt.errorbar(x, y, xerr = xerr, yerr = yerr, marker = '.', ls = 'none', color = color, label=label)
        plt.legend()
    else:
        plt.errorbar(x, y, xerr = xerr, yerr = yerr, marker = '.', ls = 'none', color = color)
    plt.xlabel("Angle (degrees)")
    if mode == 'processed':
        plt.ylabel('CPS * Energy^2')
    elif mode == 'raw':
        plt.ylabel('CPS')
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
    
    print("old: " + str(Eerror/Evalue))

    Evalue = 25.8435
    Eerror = math.sqrt(4.4640)
    newy = Evalue * np.array(y)


    yerr_rescaled = Evalue * np.array(yerr)

    print("old: " + str(Eerror/Evalue))

    energy_error = []
    for i in y:
        energy_error.append(i*Eerror)
    
    print('EnergyError: ')
    print(energy_error)

    newyerr = []
    for i in range(len(yerr_rescaled)):
        new_err = math.sqrt(yerr_rescaled[i]**2 + (energy_error[i])**2)
        newyerr.append(new_err)

    data_new = [x, newy, xerr, newyerr]

    print('Angle: ')
    print(x)
    print('CPS * Energy^2: ')
    print(newy)

    if plot:
        plot_data(data_new, show=True, title="New data (energy errors added)")

    # Propogation of x-errors
    slope_function = approximate_conv_slopes(profile, data = data_new, choice = 'rutherford', plotConv=True)
    if plot:
        x_vals=np.linspace(min(x), max(x), 1000)
        y_eval=[]
        for i in x_vals:
            y_eval.append(slope_function(i))
        plt.plot(x_vals, y_eval, label='slope')
        plt.legend()
        plt.show()

    angle_y_errors = []
    for i in range(len(x)):
        angle_y_errors.append(xerr[i] * slope_function(x[i]))
    
    print("Angle Errors: ")
    print(angle_y_errors)
    
    final_yerrs = []
    final_x_errs = []
    for i in range(len(newyerr)):
        final_x_errs.append(0)
        final_yerr = math.sqrt(newyerr[i]**2 + (angle_y_errors[i]**2))
        final_yerrs.append(final_yerr)
    
    print("Final y errors: ")
    print(final_yerrs)
    
    '''
    #Normalization
    total = 0
    for i in newy:
        total += i
    '''
    
    #final_data = [x, (1/total)* np.array(newy), final_x_errs, (1/total)* np.array(final_yerrs)]
    final_data = [x,np.array(newy), final_x_errs, np.array(final_yerrs)]
    if plot:

        plot_data(final_data, show=True, title="Processed Data (All errors propogated)")

    return final_data



def approximate_conv_slopes(profile, data, choice = 'rutherford', min_angle=10, plotConv=False):
    '''
    Using a basic profile fit, convolves with rutherford scattering to get a slope estimation at any point
    Note that this uses convolution1 not convolution2
    '''
    if choice == 'rutherford':
        convolution, domain, angles = convolve_with_rutherford1(profile, min_angle)
    elif choice == 'plum pudding':
        convolution, domain, angles = convolve_with_plum1(profile, min_angle)


    
    function, domain = interpolate(convolution, domain)
    if choice == 'rutherford':
        model = lambda x, a=(1e-10): a * np.array(function(x))
    elif choice == 'plum pudding':
        model = lambda x: np.array(function(x))
    result = fit_to_scattering(data, model, plot=False)

    a = result.params['a'].value
    print(a)

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

def convolve_with_plum1(profile, min_angle = 10):
    '''
    implements the convolve method in models to convolve a beam profile with the rutherford model
    returns the resulting array (not a callable), the domain in which it is valid, and the angles 
    '''
    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 180, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                             # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)

    scattering,pdomain = convolve(profile, bpdomain, plum_model, angles)

    return scattering, pdomain, angles


def convolve_with_rutherford2(profile):
    '''
    implements the second version of convolve, returns a callable function
    '''
    return lambda x : convolve2(f, profile, x)

def convolve_with_plum2(profile):
    '''
    implements the second version of convolve, returns a callable function
    '''
    return lambda x : convolve2(plum_model, profile, x)


def get_convolutions2(profiles, choice = 'rutherford', min_angle = 10, plot=True):
    '''
    implements convolve_with_rutherford2 for each profile and returns an array of convolutions (callable) given a list of profiles
    '''
    convolutions = []
    for profile in profiles:
        if choice == 'rutherford':
            convolutions.append(convolve_with_rutherford2(profile))
        elif choice == 'plum pudding':
            convolutions.append(convolve_with_plum2(profile))


    if plot:
        pdomain = np.linspace(min_angle, 60, 1000)

        for i in range(len(convolutions)):
            conv = convolutions[i]
            plt.plot(pdomain, conv(pdomain), label = "convolution " + str(i))

        plt.title("Convolutions of beam profile with " + choice + " Scattering Expectation 2")
        if choice == 'rutherford':
            plt.plot(pdomain, f(pdomain), label = "scattering expectation", color = "red")
        elif choice == 'plum pudding':
            plt.plot(pdomain, plum_model(pdomain), label = "scattering expectation", color = "red")

        plt.legend()
        plt.show()


        residuals = []
        for i in range(len(pdomain)):
            if choice == 'plum pudding':
                residuals.append(conv(i) - plum_model(pdomain[i]))
            elif choice == 'rutherford':
                residuals.append(conv(i) - f(pdomain[i]))

        plt.plot(pdomain, residuals)
        plt.title('residuals :)')
        plt.show()

    return convolutions



def get_convolutions1(profile_sets, choice = 'rutherford', min_angle = 10, plot=False):
    '''
    given a bunch of profile sets, convolve with the selected choice and resutnrs the results
    '''
    convolutions = []
    domains = []
    stepsize = 0.01
    bpdomain = np.arange(-10, 10, stepsize)
    angles = np.arange(min_angle, 90, stepsize) # note: to make convolution work properly, recommend we do the cutoff thing for the function `f` below 
                                    # and then change `angles` to np.arange(min_angle - 20, 180, stepsize)

    

    for i in range(len(profile_sets)):
        profile = profile_sets[i]

        if choice == 'rutherford':
            scattering,pdomain = convolve(profile, bpdomain, f, angles)
        elif choice == 'plum pudding':
            scattering,pdomain = convolve(profile, bpdomain, plum_model, angles)
        else:
            raise Exception('not a valid choice')
        
        if choice == 'rutherford':
            eval = f(pdomain)
        elif choice == 'plum pudding':
            eval = plum_model(pdomain)

        convolution = np.array(scattering)* (np.sum(eval)/np.sum(scattering))
        convolutions.append(convolution) # Treat the original model as "normalized" Normalize the convolution w/ respect to it
        domains.append(pdomain)
        if plot:
            plt.plot(pdomain, convolution, label = "convolution " + str(i))

    

    if plot:
        '''
        scattering_domain = np.linspace(3, 60, 1000)
        for i in scattering_domain:
            if choice == 'rutherford':
                eval.append(f(i))
            elif choice == 'plum pudding':
                eval.append(plum_model(i))
            else:
                raise Exception('not a valid choice')
        '''
        plt.plot(pdomain, eval, label = choice + " scattering expectation", color = "black", linestyle='dashed')
        plt.ylabel('(Unnormalized) Scattering Cross Section')
        plt.xlabel('Angle (degrees)')
        plt.legend()
        plt.show()

        residuals = []
        
        for i in range(len(pdomain)):
            if choice == 'plum pudding':
                residuals.append(convolution[i] - eval[i])
            elif choice == 'rutherford':
                residuals.append(convolution[i] - eval[i])

        plt.plot(pdomain, residuals)
        plt.title('residuals (Convolved Function - Original Function)')
        plt.ylabel('(Unnormalized) Scattering Cross Section')
        plt.xlabel('Angle (degrees)')
        plt.show()

    return convolutions, domains

'''
'''
'''
FITTING METHODS
'''
'''
'''
#def scattering_model(x, a = 1e-5, b = 0):
def scattering_model(x, a = 1e-5):
    '''
    the model for pure (no convolution) rutherford scattering
    '''
    try:
        #return a * f(x) + b
        return a * f(x)
    except:
        y = []
        for i in x:
            y.append(scattering_model(i))
        return y

def plum_model(x):
    return (2/math.sqrt(math.pi)) * np.exp(-(np.power(x,2)))


def plot_models(data, rutherford_conv, plum_conv, domain_r = None, domain_p = None):
    # Fit the rutherford model
    function, domain = interpolate(rutherford_conv, domain_r)
    model = lambda x, a=(1e-10): a * np.array(function(x))
    result = fit_to_scattering(data, model, plot=False)

    plum_func, plum_domain = interpolate(plum_conv, domain_p)

    a = result.params['a'].value

    # Plot the data
    plot_data(data, show=False, label='data')


    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]
    #Plot the unconvolved functions
    x_vals = np.linspace(0.2, 55, 1000)
    plum_evals = []
    rutherford_evals = []

    for i in x_vals:
        plum_evals.append(plum_model(i))
    for i in x_vals:
        rutherford_evals.append(a * f(i))
    
    plt.plot(x_vals, plum_evals, label='plum scattering')


    index = 0
    for i in range(len(rutherford_evals)):
        if rutherford_evals[i] > 1:
            index = i
    
    x_vals_r = x_vals[index:]
    rutherford_evals_new = rutherford_evals[index:]

    plt.plot(x_vals_r, rutherford_evals_new, label='rutherford scattering')


    # Plot the convolutions for each
    x_r = np.linspace(min(domain_r) + 0.01, 55, 1000)
    x_p = np.linspace(min(domain_p) + 0.01, 55, 1000)

    r_conv_eval = []
    p_conv_eval = []

    for i in x_r:
        r_conv_eval.append(a * function(i))
    for i in x_p:
        p_conv_eval.append(plum_func(i))

    index = 0
    for i in range(len(r_conv_eval)):
        if r_conv_eval[i] > 1:
            index = i
    
    r_conv_eval_new = r_conv_eval[index:]
    x_r_new = x_r[index:]
    
    plt.plot(x_r_new, r_conv_eval_new, label='rutherford convolution')

    index = 0
    for i in range(len(p_conv_eval)):
        if p_conv_eval[i] > 1:
            index = i
    
    p_conv_eval_new = p_conv_eval[index:]
    x_p_new = x_p[index:]

    plt.plot(x_p_new, p_conv_eval_new, label='plum convolution')


    #plt.yticks((1/10)* np.array(range(-1, 10)))
    plt.legend()
    plt.show()


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
    
    if plot:
        plot_fit(result, model, element, min_angle, folder, show=True)

    return result

def fits_chi2(data, convolutions, choice = 'rutherford', domains = None, plot=True, plotFit = True):
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
            if choice == 'rutherford':
                model = lambda x, a=(1e-5): a * np.array(func(x))
            elif choice == 'plum pudding':
                model = lambda x: np.array(func(x))
            result = fit_to_scattering(data, model, plot=False)

        except:
            conv = convolutions[i]
            domain = domains[i]
            function, domain = interpolate(conv, domain)        
            #model = lambda x, a=(1e-10), b=0 : a * np.array(function(x)) + b
            model = lambda x, a=1: a * np.array(function(x))
            result = fit_to_scattering(data, model, plot=False)
        
        chi = result.redchi
        chi2.append(chi)
    if plot:
        plt.hist(chi2)
        plt.xlabel('reduced chi square')
        plt.ylabel('frequency')
        plt.show()
    
    if plotFit:
        plot_fit(result, model, data, show=True, residuals=True, label=choice)
    
    return chi2


def plot_fit(result, model, data, show=True, label=None, initial=False, residuals=False):
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
        y_eval.append(model(i, params['a'].value))
    if label is not None:
        plt.plot(x_vals, y_eval, label='fit ' + label)
    else:
        plt.plot(x_vals, y_eval, label='fit')

    if initial:
        initial_y= []
        for i in x_vals:
            initial_y.append(model(i, init_params['a'].value))
        if label is not None:
            plt.plot(x_vals, initial_y, label = 'initial guess ' + label)
        else:
            plt.plot(x_vals, initial_y, label = 'initial guess')


    if show:
        plt.legend()
        plt.show()

        if residuals:
            residuals = []
            for i in range(len(y)):
                residuals.append((y[i] - model(x[i], params['a'].value))/yerr[i])
            plt.errorbar(x, residuals, np.array(yerr)/np.array(yerr), fmt='o')
            plt.hlines(0, 18, 50, colors='red', linestyles='dashed')
            plt.xlabel('Angle (degrees)')
            plt.ylabel('Residuals (CPS * Energy^2)')
            plt.title('Residuals Between Scattering Data and Fitted Convolution')
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


def plum_scattering_fit(data, plot=False):
    '''
    pure plum (no convolution) fit to  data
    again, data must be of the form list: [x, y, xerr, yerr]
    '''
    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    themodel = lmfit.Model(plum_model)
    weights = []
    for i in yerr:
        weights.append(1/i)
    result2  = themodel.fit(y, x=x, weights=weights)

    if plot:
        plot_fit(result2, plum_model, element, min_angle, folder, show=True, label='rutherford scattering')
    
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
        #model = lambda x, a=(1e-10), b=0 : a * np.array(function(x)) + b
        model = lambda x, a=(1e-10): a * np.array(function(x))
        result = fit_to_scattering(data, model, plot=False)
    except:
        #model = lambda x, a=(1e-10), b=0 : a * np.array(rutherford_convolution(x)) + b
        model = lambda x, a=(1e-10): a * np.array(rutherford_convolution(x))
        result = fit_to_scattering(data, model, plot=False)
    
    print("Fit complete. Plotting Now")
    plot_fit(result, model, data, show=True, label='beam profile x rutherford', residuals=True)

    
    result2 = rutherford_scattering_fit(data, plot=False)
    plot_fit(result2, scattering_model, data, show=False, label='rutherford scattering')
    
    
    x = data[0]
    y = data[1]
    xerr = data[2]
    yerr = data[3]

    a = result.params['a']
    #b = result.params['b']
    a2 = result2.params['a']
    x_vals=np.linspace(min(x), max(x), 1000)
    y_eval=[]
    for i in x_vals:
        y_eval.append(scattering_model(i, a))
    plt.plot(x_vals, y_eval, label='unconvolved')
    

    plt.legend()
    plt.show()

    


def compare_chi2(data, rutherford_convolutions, plum_convolutions, rutherford_domains=None, plum_domains=None):
    '''
    Plots the distribution of chi2 for rutherford convolutions and a vertical line for 
    no convolution fit chi2
    '''
    no_conv_chi = rutherford_scattering_fit(data).redchi

    #plt.vlines(no_conv_chi, ymin=0, ymax=len(rutherford_convolutions)/4, label='no convolution', linestyles='dashed', color='red')

    '''
    rutherford_chi = fits_chi2(data_r, rutherford_convolutions, rutherford_domains, choice = 'rutherford', plot=False)
    plum_chi = fits_chi2(data_p, plum_convolutions, plum_domains, choice = 'plum pudding', plot=False)
    '''

    rutherford_chi = fits_chi2(data, rutherford_convolutions, choice = 'rutherford', plot=False, domains=rutherford_domains)
    #plum_chi = fits_chi2(data, plum_convolutions, choice = 'plum pudding', plot=False, domains=plum_domains)

    #print(plum_chi[0])

    plt.hist(rutherford_chi, label='rutherford convolutions')
    #plt.hist(plum_chi, label = 'plum convolutions')
    plt.xlabel('reduced chi square')
    plt.ylabel('frequency')
    plt.legend()
    plt.show()

    sum = 0
    for i in rutherford_chi:
        sum += i
    mean = sum/len(rutherford_chi)
    print('mean: ' + str(mean))

    sum = 0
    for i in rutherford_chi:
        sum += (i - mean)**2
    std_dev = math.sqrt(sum/len(rutherford_chi))
    print('standard deviation: ' + str(std_dev))

    


