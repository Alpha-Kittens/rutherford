import os
from data_loader import read_data
import matplotlib.pyplot as plt
import math
import lmfit
import numpy as np



def left_beam_model(x, a, b, c, choice = 'linear'):
    if choice == 'exponential':
        return c * math.exp(a*(x-b))
    if choice == 'linear':
        return a*x + b
    if choice == 'quadratic':
        return a*(x**2) + b*x + c



def beam_model(x, x_0 = 0, a=0.5, b=-4, c = 20, d=-20):
    '''
    Test beam profiles
    choice = exponential should could have initial values
    
    '''
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