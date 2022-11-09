import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np

'''
Contains methods useful for fitting the beam profile
The useful method to call is beam_profile_fit
All other methods are helpers

To add another model choice:
 * add the function to evaluate_function_choice
 * add parameters to get_params (make sure not to give any parameters names with underscores)
'''


def evaluate_function_choice(x, choice, params, x0, y0):
    '''
    evaluates a chosen function, given its parameters

    Arguments: 
        * `x` (float) : the point at which the evaluate the function
        * `choice` (string) : what type of function you want to evaluate (e.g. exponential, linear, etc)
        * `params` (dictionary) : a dictionary witht the parameters needed in order to evaluate this function
        * `x0` the "origin" where you split between left, and right models
        * `y0` the value of the function at x0
    Returns: 
        * (float) result of evaluating the function f(x)
    '''
    if choice == 'exponential':
        a = params['a']
        b= params['b']
        return y0 * math.exp(a*((x-x0)-b))
    if choice == 'linear':
        a = params['a']
        return a*(x-x0) + y0
    if choice == 'quadratic':
        a = params['a']
        b= params['b']
        return a*((x-x0)**2) + b*(x-x0) + y0

def get_params(choice, side):
    '''
    returns a parameters objects for the given choice, also sets initial values of parameters

    Arguments: 
        * `choice` (string) : what type of function you want to parameters for (e.g. exponential, linear, etc) 
        * `side` (string) : L or R indicatin wether you want parameters for the left side or right side model
    
    Returns: 
        * (Parameters) parameters object for the selected choice and side
    '''
    factor = 1
    if(side == "L"):
        factor = -1

    params = lmfit.Parameters()
    if choice == 'exponential':
        params.add('a' + "_" + str(side), value=-0.5*factor)
        params.add('b' + "_" + str(side), value=0)
    if choice == 'linear':
        params.add('a' + "_" + str(side), value=-20*factor)
    if choice == 'quadratic':
        params.add('a' + "_" + str(side), value=-0.5*factor)
        params.add('b' + "_" + str(side), value=-20*factor)

    # Since the right side is constrained by the left, dont't add this parameter for right side model
    if(side == "L"):
        params.add('y0' + "_" + str(side), value=150)


    return params

def get_init_params(choice_L, choice_R):
    '''
    creates the parameters object to pass to lmfit.fit
    combines left and right parameters

    Arguments: 
        * `choice_L` (string) : choice for left side function
        * `choice_R` (string) : choice for right side function

    Returns: 
        * (Parameters) parameters object containing all parameters for the full beam profile model
    '''

    params = lmfit.Parameters()

    params.add('x0', value=0)
    left_params = get_params(choice_L, side='L')
    right_params = get_params(choice_R, side='R')

    for key in left_params:
        params.add(key, value=left_params[key])
    for key in right_params:
        params.add(key, value=right_params[key])

    return params



def evaluate_beam_model(x, choice_L, choice_R, params):
    '''
    evaluates the full beam_profile model

    Arguments: 
        * `x` (float) :  the point at which to evaluate the model
        * `choice_L` (string) : choice for left side function
        * `choice_R` (string) : choice for right side function
        * `params` (dictionary) : the dictionary of parameters and their value

    Returns: 
        * (float) the result of evaluating the beam profile model
    '''

    left_params = lmfit.Parameters()
    right_params = lmfit.Parameters()

    for key in params:
        if key != 'x0':
            if 'L' not in key:
                split = key.split('_')
                right_params.add(split[0], value=params[key])
            else:
                split = key.split('_')
                left_params.add(split[0], value=params[key])

    x0 = params['x0']

    if x < x0:
        answer = evaluate_function_choice(x, choice=choice_L, params=left_params, x0=x0, y0=left_params['y0'])
        if answer > 0:
            return answer
        else:
            return 0
    else:
        h = evaluate_function_choice(x0, choice=choice_L, params=left_params, x0=x0, y0=left_params['y0']) #height at the center
        answer = evaluate_function_choice(x, choice=choice_R, params=right_params, x0=x0, y0=h)
        if answer > 0:
            return answer
        else:
            return 0


def beam_model(x, choice_L, choice_R, **params):
    '''
    the actual beam model, implements evaluate_beam_model

    Arguments: 
        * `x` (float or list of floats) :  the point at which to evaluate the model or a list of points
        * `choice_L` (string) : choice for left side function
        * `choice_R` (string) : choice for right side function
        * `params` (floats) : the neccesary parameters for this model

    Returns: 
        * (float or list of floats) the result of evaluating the beam profile model
    '''
    try: 
        return evaluate_beam_model(x, choice_L, choice_R, params)
    except:
        y_values = []
        for i in range(0, len(x)):
            y_values.append(evaluate_beam_model(x[i], choice_L, choice_R, params))
        return y_values


#The actual beam model for selected choices
the_beam_model = lambda choiceL, choiceR : lambda x, **params : beam_model(x, choiceL, choiceR, **params)


def fit_report(result, choiceL, choiceR):
    '''
    given the results of fit, it prints important features nicely

    Arguments: 
        * `result` (ModelResult) :  the resut of fittings
        * `choiceL` (string) : choice for left side function
        * `choiceR` (string) : choice for right side function

    '''

    params = result.params
    init_params = result.init_params


    print("-----------")
    print("FIT REPORT")
    print("-----------")

    print("x0 (init = " + str(init_params['x0'].value) + ') : ' + str(params['x0'].value) + ' +/- ' + str(params['x0'].stderr))
    print("-----------")
    print("Left side model: " + choiceL)
    for param in params:
        if 'R' not in param and param !='x0':
            split = param.split('_')
            name = split[0]
            print(name + ' (init = ' + str(float(init_params[param].value)) + ') : ' + str(params[param].value) + ' +/- ' + str(params[param].stderr))
    print("-----------")
    print("Right side model: " + choiceR)
    for param in params:
        if 'L' not in param and param !='x0':
            split = param.split('_')
            name = split[0]
            print(name + ' (init = ' + str(float(init_params[param].value)) + ') : ' + str(params[param].value) + ' +/- ' + str(params[param].stderr))
    print("-----------")
    print('chi2: ' + str(result.chisqr))
    print('reduced chi2: ' + str(result.redchi))
    

def plot_fit(x, y, yerr, result, choiceL, choiceR):
    '''
    given the results of fit, plots the raw data, the fit, and initial guess

    Arguments: 
        * `x` the x coordinates of the raw data
        * `y` the y coordinates of the raw data
        * `yerr` the errors on the y coordinates
        * `result` (ModelResult) :  the resut of fittings
        * `choiceL` (string) : choice for left side function
        * `choiceR` (string) : choice for right side function
    '''

    plt.errorbar(x, y, yerr=yerr, marker='o', ls='none', label='data')

    params = result.params
    init_params = result.init_params

    x_vals=np.linspace(min(x), max(x), 1000)
    y_eval=[]
    initial_y= []

    for x in x_vals:
        y_eval.append(evaluate_beam_model(x, choiceL, choiceR, params))
        initial_y.append(evaluate_beam_model(x, choiceL, choiceR, params=init_params))

    plt.plot(x_vals, y_eval, label='fit')
    plt.plot(x_vals, initial_y, label='initial guess')

    plt.legend()
    plt.show()

def unpack_params(raw_params):
    params = {}
    for name, param in raw_params.items():
        params[name] = param.value
    return params

def beam_profile_fit(x, y, yerr, choiceL, choiceR, plot=False):
    '''
    implements lmfit.fit for the selected model to fit the beam profile data to a function

    Arguments: 
        * `x` the x coordinates of the raw data
        * `y` the y coordinates of the raw data
        * `yerr` the errors on the y coordinates
        * `result` (ModelResult) :  the resut of fittings
        * `choiceL` (string) : choice for left side function
        * `choiceR` (string) : choice for right side function
    '''

    model = lmfit.Model(the_beam_model(choiceL, choiceR))

    weights = []
    for i in yerr:
        if i==0:
            weights.append(float(0))
        else:
            weights.append(1/i)

    weights = np.array(weights)

    params = get_init_params(choiceL, choiceR)
    
    result = model.fit(y, x=np.array(x), params = params, weights=weights)

    fit_report(result, choiceL, choiceR)

<<<<<<< Updated upstream
    return lambda x : the_beam_model(choiceL, choiceR)(x, **unpack_params(result.params))


    
=======
    if plot:
        plot_fit(x, y, yerr, result, choiceL, choiceR)

    return result

    
def profile_sys_error(angles, cpss, errors, angle_error, choiceL = 'linear', choiceR= 'linear'):
    left_angles = angles
    right_angles = angles
    for i in range(len(angles)):
        left_angles[i] -= angle_error
        right_angles[i] += angle_error
        
    result_L = beam_profile_fit(left_angles, cpss, errors, choiceL = 'linear', choiceR= 'linear', plot=True)
    result_R = beam_profile_fit(right_angles, cpss, errors, choiceL = 'linear', choiceR= 'linear', plot=True)

    return result_L.params, result_R.params
>>>>>>> Stashed changes
