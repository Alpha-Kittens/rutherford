import os
from data_loader import read_data
import matplotlib.pyplot as plt
import numpy as np
import math


folder = 'gold_scattering/'
file = 'gold_0.Spe'

fp = folder+file

information = read_data(fp)

histogram = information['histogram']



def max_model(histogram, plot = False, target = None, threshold = 1/2, give_weights = False):

    data = np.zeros(shape=(len(histogram), 2))
    weights = []

    for i in range(len(histogram)):
        data[i, 0] = i
        data[i, 1] = histogram[i]

        if(histogram[i] == 0):
            weights.append(math.inf)
        else:
            weights.append(1/math.sqrt(histogram[i]))
    

    data_half = (max(data[:,1]) + min(data[:,1])) * threshold

    weights = np.array(weights)

    data_above = data[:,1] > data_half
    data_below = data[:,1] <= data_half
    cut_wavelengths = data[data_above,0]
    cut_weights = weights[data_above]
    
    averaging_weights = data[data_above,1]
    
    w_max = np.dot(cut_wavelengths, averaging_weights) / np.sum(averaging_weights)
    w_err = np.sqrt(np.dot((cut_wavelengths - w_max)**2, averaging_weights) / np.sum(averaging_weights))

    if plot:
        plt.title("max model on scan for "+str(target))
        plt.xlabel("Channel Number")
        plt.ylabel("CPS")
        plt.errorbar(data[data_below,0], data[data_below,1], yerr = weights[data_below], label = "excluded data", marker = '.', c = 'b', ls='none')
        plt.errorbar(cut_wavelengths, data[data_above,1], yerr = cut_weights, label = "included data", marker = '.', c = 'orange', ls='none')
        plt.axvline(x = w_max, label = "wavelength estimate", c = 'r', linestyle = '--')
        plt.axvline(x = w_max - w_err, label = "error bounds", c = 'magenta', linestyle = '--')
        plt.axvline(x = w_max + w_err, c = 'magenta', linestyle = '--')
        plt.legend()
        plt.show()

    if give_weights:
        return w_max, w_err, averaging_weights
    else:
        return w_max, w_err


if __name__ == '__main__':
    max_model(histogram, plot=True, target=information['target'])