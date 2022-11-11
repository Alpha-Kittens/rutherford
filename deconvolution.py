import numpy as np
from scipy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt 

def deconvolve (output, input):
    H = fft(kill_noise(output))
    G = fft(kill_noise(input))
    F = H / G
    return ifft(F)

def reduce_noise (histogram, damping_constant = 1/10, plot = False, title=None):
    """
    Given data, uses a fourier bandpass filter to reduce noise and estimate uncertainty of each counts per second measurement. 
    
    Arguments:
        * `data` (nx2 numpy array): 2 columns. First is wavelength in angstroms, second is measured counts per second. 
    Keyword Arguments:
        * `plot` (Boolean): Whether the results should be displayed in a plot. Defaults to false. 

    Returns:
        * As a 2-tuple:
            - `new_data` (nx2 numpy array): numpy array with the same format as `data`, hopefully with reduced noise
            - `weights` (numpy array length n): weights associated with each new data point.  
    """

    

    amplitudes = fft(histogram)

    new_amplitudes, removed_amplitudes = sophisticated_bandpass(amplitudes, damping_constant = damping_constant)

    new_data = ifft(new_amplitudes)

    if plot:
        psr = primary_signal_region(histogram)
        #plt.scatter(x[psr], y[psr], marker = '.')
        #plt.plot(x[psr], new_data[psr])
        #plt.scatter(x, y, marker = '.', c = 'b')
        plt.scatter(x, y, marker = '.', label="raw data")
        #plt.plot(x, new_data, c = 'r')
        plt.plot(x, new_data, c = 'red', label="noise reduced")
        plt.xlabel('Monochromator Step')
        plt.ylabel('Counts per second')        
        plt.legend()
        if title is not None:
            plt.title(title)
        plt.show()

    return np.real(new_data)

def sophisticated_bandpass(amplitudes, damping_constant = 1/10):
    """
    Given frequencies and amplitudes, uses the cutoff method from max_model to retain the dominant frequencies and damp the 
    higher frequencies with less amplitude using an exponential model. 

    Arguments:
        * `amplitudes` (numpy array): Amplitudes of fourier transform of data. 
    Keyword Arguments: 
        * `damping_constant` (float): Positive number for the damping factor e^(-`damping_constant` * distance), where distance
            is the distance (as measured by index in the amplitudes array, which is ordered by fequency) between the point and the
            nearest "signal region" as determined by the cutoff algorithm. 

    Returns:
        * As a 2-tuple:
            `new_amplitudes` (numpy array of same dimension as `amplitudes`): Array of new amplitudes. 
            `removed_amplitudes` (numpy array of same dimension as `amplitudes`): Array of removed amplitudes. These two should sum to the original. 
    """

    cutoff = get_cutoff(np.abs(amplitudes))
    signals, backgrounds = get_regions(np.abs(amplitudes), cutoff, reduce = False)

    #new_amplitudes = np.where(np.abs(amplitudes) > cutoff, amplitudes, 0)

    new_amplitudes = np.array([amplitudes[i] if abs(amplitudes[i]) > cutoff else amplitudes[i] * np.exp(-damping_constant * get_distance_to_nearest(signals, i)) for i in range(len(amplitudes))])

    return new_amplitudes, amplitudes - new_amplitudes

def get_distance_to_nearest(signals, i):
    """
    Used by sophisticated bandpass.
    """

    min = -1

    for signal in signals:

        if min == -1 or abs(signal[0] - i) < min:
            min = abs(signal[0] - i)
        if min == -1 or abs(i - signal[1]) < min:
            min = abs(signal[1] - i)

    
    
    return min

def regions(histogram, cutoff = 2, reduce = True):
    """
    Given an array of positive values, determines a cutoff to distinguish signal from noise, then determines regions of indices
    which correspond to background and signal data. 
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
    Keyword Arguments:
        * `reduce` (Boolean): True if regions containing only one point should be removed; False otherwise. 
            Generally should be set to True for wide signals.  

    Returns:
        * As a 2-array
            - `backgrounds` (Array of tuples): Ranges of indices of `cps` classified as background
            - `signals` (Array of tuples):  Ranges of indices of `cps` classified as signal
    """


    return get_regions(histogram, cutoff, reduce = reduce)




def get_regions(cps, cutoff, reduce = True):
    """
    Given an array of positive values and a cutoff, determines signal and background regions of the data. 
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
        * `cutoff` (float): Classifying value. values of `cps` above `cutoff` are signals; otherwise background.
    Keyword Arguments:
        * `reduce` (Boolean): True if regions containing only one point should be removed; False otherwise. 
            Generally should be set to True for wide signals.  

    Returns:
        * As a 2-array 
            - `backgrounds` (Array of tuples): Ranges of indices of `cps` classified as background
            - `signals` (Array of tuples):  Ranges of indices of `cps` classified as signal
    """

    background_regions = []
    signal_regions = []

    data_below = cps <= cutoff
    
    tracking_background = True
    start = 0

    for i in range(len(cps)):

        if data_below[i] != tracking_background:
            if tracking_background:
                if (not reduce or (i-1) - start > 1) and i-1 >= 0:
                    background_regions.append((start, i-1))
            else:
                if not reduce or (i-1) - start > 1:
                    signal_regions.append((start, i-1))
            tracking_background = not tracking_background
            start = i

    if tracking_background:
        background_regions.append((start, len(cps) - 1))
    else:
        signal_regions.append((start, len(cps) - 1)) 


    return background_regions, signal_regions




# max model, HWHM cutoff
def get_cutoff(cps):
    """
    Given an array of positive values, determines a cutoff do distinguish signal from noise. This algorithm is not particularly refined.
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
    Returns:
        * `cutoff` (float): signal/noise cutoff for `cps` array. 
    """

    hist, binedges = np.histogram(cps, int(max(cps) // min(cps)))

    cutoff = 1 + int(max(cps) // min(cps) / 50)

    #print ("The cutoff is ", binedges[cutoff], ". ", np.sum(hist[cutoff:]), "/", np.sum(hist), " of the data is above the cutoff.")

    return binedges[cutoff]

def primary_signal_region (histogram, cutoff = 2):

    backgrounds, signals = regions(histogram, cutoff)

    maxval = 0
    max_signal = (0,1)
    for signal in signals:
        smax = max(histogram[signal[0]:signal[1]+1])
        if smax > maxval:
            maxval = smax
            max_signal = signal
    
    i = np.array(range(len(histogram)))
    return np.logical_and(i >= max_signal[0], i <= max_signal[1])

def smear (arr, width = 5):
    return [smear_i(arr, i, width) for i in range(len(arr))]
    
def smear_i(arr, i, width):
    min_i = max(i - width, 0)
    max_i = min(i + width, len(arr) - 1) + 1 # endpoints don't matter at all
    dists = np.array(range(min_i - i, max_i - i))
    return np.dot(arr[min_i : max_i], (1 - (dists / width)**2)) / np.sum(1 - (dists / width)**2)

def kill_noise(histogram):
    h = np.array(histogram)[:]
    psr = primary_signal_region(h)
    h_nr = reduce_noise(h[psr])
    h[psr] = h_nr
    return smear(np.where(psr, h, 0))

if __name__ == '__main__':
    from data_loader import recursive_read
    r = {}
    recursive_read(r, "data", require = [0], reject = ["titanium", "2gold", "3gold", "iron"])
    output = np.array(r[("gold", 0, 1)][1]) / r[("gold", 0, 1)][0]
    input = np.array(r[("empty", 0, 1)][1]) / r[("empty", 0, 1)][0]
    x = np.array(range(len(output)))

    f = deconvolve(output, input)

    print (f)
    plt.plot(x, output, label = "output", color = 'r')
    plt.plot(x, input, label = "input", color = 'b')
    plt.plot(x, f, label = "f", color = 'purple')
    plt.legend()
    plt.show()