import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np


def evaluate_beam_model(x, choice, params):
    if choice == 'exponential':
        return c * math.exp(a*(x-b))
    if choice == 'linear':
        return a*x + b
    if choice == 'quadratic':
        return a*(x**2) + b*x + c

def get_params(choice, side):
    '''
    returns a parameters objects for the given choice

    Arguements:
    `choice` : a string designating the model choice
    `side` : integer, either 1 or -1. 1 means right side, -1 means left side (which side of the data these paramaters are for)
    '''
    params = lmfit.Parameters()
    if choice == 'exponential':
        params.add('a' + " " + str(side), value=-0.5*side)
        params.add('b' + " " + str(side), value=4 * side)
    if choice == 'linear':
        params.add('a' + " " + str(side), value=-20*side)
    
    if(side == -1):
        params.add('y_int' + " " + str(side), value=50)


    return params

def beam_model(x, choice_L, choice_R, **params):
    choice = 'exponential' #model choice for the left side

    left_params = lmfit.Parameters()
    right_params = lmfit.Parameters()

    for key in params:
        if key != 'x_0':
            if '-' not in key:
                right_params.add(key, value=params[key])
            else:
                left_params.add(key, value=params[key])

    x_0 = params['x_0']

    h = left_beam_model(x_0, a, b, c, choice=choice) #height at the center


    try:  
        if x < x_0:
            return evaluate_beam_model(x, choice=choice_L, params=left_params)
        else:
            return evaluate_beam_model(x, choice=choice_R, params=right_params)
    except:
        y_values = []
        for i in range(0, len(x)):
            y_values.append(beam_model(x[i], x_0, a, b, c, d))
        return y_values



'''
def beam_model(x, x_0 = 0, a=0.5, b=-4, c = 20, d=-20):
    choice = 'exponential' #model choice for the left side
    h = left_beam_model(x_0, a, b, c, choice=choice) #height at the center

    try:  
        if x < x_0:
            return left_beam_model(x, a, b, c, choice=choice)
        else:
            return h + (d*x)
    except:
        y_values = []
        for i in range(0, len(x)):
            y_values.append(beam_model(x[i], x_0, a, b, c, d))
        return y_values
'''

def beam_profile_fit(x, y, yerr):
    model = lmfit.Model(beam_model)

    weights = []
    for i in yerr:
        if i==0:
            weights.append(float(0))
        else:
            weights.append(1/i)

    weights = np.array(weights)
    
    result = model.fit(y, x=np.array(x), weights=weights)
    print(lmfit.fit_report(result))
    
    plt.errorbar(x, y, yerr=yerr, xerr = 0.5, marker='o', ls='none', label='data')
    x_vals=np.linspace(min(x), max(x), 1000)
    y_eval=[]
    intial_y= []
    for i in x_vals:
        y_eval.append(beam_model(i, result.best_values['x_0'], result.best_values['a'], result.best_values['b'], result.best_values['c'], result.best_values['d']))
        intial_y.append(beam_model(i))
    plt.plot(x_vals, y_eval, label='fit')
    plt.plot(x_vals, intial_y, label='initial guess')

    plt.ylabel('cps')
    plt.xlabel('angle')
    plt.legend()
    plt.show()

    residuals = []
    for i in range(0, len(x)):
        predicted_y = beam_model(x[i], result.best_values['x_0'], result.best_values['a'], result.best_values['b'], result.best_values['c'], result.best_values['d'])
        residuals.append(y[i] - predicted_y)
    plt.errorbar(x, residuals, yerr=yerr, xerr = 0.5, marker='o', ls='none')
    plt.show()