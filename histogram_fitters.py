import lmfit
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import moyal
from data_processing import *


"""
def convolve(x,f1,f2,iMin=-10,iMax=10,iN=2000):
    convolution = []
    for xi in x:
        step=(iMax-iMin)/iN
        pInt=0
        for i0 in range(iN):
                pX   = i0*step+iMin
                pVal = f1(xi-pX)*f2(pX)
                pInt += pVal*step
                #print (pInt)
        convolution.append(pInt)
    return convolution"""
def convolve(f, fdomain, g, gdomain):

    # everything needs to be sorted
    # works better if the shorter of fdomain and gdomain is of odd length, I think?
    M, N = (gdomain, fdomain) if len(fdomain) > len(gdomain) else (fdomain, gdomain)
    a = math.floor((len(M)-1) / 2) #3: a=b=1. 4: a=1,b=2
    b = math.ceil((len(M)-1) / 2)
    return np.convolve(f(fdomain), g(gdomain))[len(M):len(N)+1], 1/2 * (N[a:-b] + N[b:-a])

def simple_convolve(f, fdomain, g, gdomain):
    return np.convolve(f(fdomain), g(gdomain), mode = 'same')


def moyal_convolution(x, f_histogram, loc, scale, test_bounds = 100):
    test_domain = loc + np.arange(-1*test_bounds, test_bounds, 1)
    test_vals = moyal.pdf(test_domain, loc = loc, scale = scale)
    #print (test_vals)
    good_domain = test_domain[test_vals > 1e-4]
    #print (len(good_domain))
    if len(good_domain) % 2 == 0:
        good_domain = np.append(good_domain, good_domain[-1]+1)
    #print (len(good_domain))
    convolution = convolve(f_histogram, x, lambda x : moyal.pdf(x, loc = loc, scale = scale), good_domain)
    padding = int((len(x) - len(convolution[0])) /2)
    return [0 for i in range(padding)] + list(convolution[0]) + [0 for i in range(padding)]

def moyal_convolution_2(x, f_histogram, loc, scale):
    return simple_convolve(f_histogram, x, lambda x : moyal.pdf(x, loc = loc, scale = scale), x)

def init_moyal_convolution(x, f_histogram, loc, scale, test_bounds = 100):
    test_domain = loc + np.arange(-1*test_bounds, 2*test_bounds, 1)
    test_vals = moyal.pdf(test_domain, loc = loc, scale = scale)
    #plt.plot(test_domain, test_vals, label = 'moyal')
    #plt.legend()
    #plt.show()
    #plt.plot(x, np.convolve(f_histogram(x), moyal.pdf(x, loc = loc, scale = scale), mode = 'same'), label = 'convolution')
    #plt.legend()
    #plt.show()
    #print (test_vals)
    good_domain = test_domain[test_vals > 1e-4]
    #print (len(good_domain))
    if len(good_domain) % 2 == 0:
        good_domain = np.append(good_domain, good_domain[-1]+1)
    #print (len(good_domain))
    convolution = simple_convolve(f_histogram, x, lambda x : moyal.pdf(x, loc = loc, scale = scale), x)
    #print (convolution)
    padding = int((len(x) - len(convolution)) /2)
    return np.array([0 for i in range(padding)] + list(convolution) + [0 for i in range(padding)])
    
    #return convolve(x, f_histogram, lambda x : moyal.pdf(x, loc = loc, scale = scale), iMin = min(good_domain), iMax = max(good_domain))

r_moyal_convolution = lambda f_histogram, results : lambda x : moyal_convolution_2(x, f_histogram, results.params['loc'].value, results.params['scale'].value)

def fit_moyal_convolution(ihistogram, fhistogram, plot = False, foil = None, include_params = False):
    x = np.array(range(len(ihistogram)))
    f_histogram = optimal_function_fit(ihistogram)
    #e_difference = (optimal_energy_fit(ihistogram) - optimal_energy_fit(fhistogram)).val
    #model = lmfit.Model(lambda x, c, loc, scale : c*np.array(moyal_convolution_2(x, f_histogram, loc, scale)))
    model = lmfit.Model(lambda x, loc, scale : np.array(moyal_convolution_2(x, f_histogram, loc, scale)))
    params = lmfit.Parameters()
    loc_init = 0
    max_overlap = 0
    for loc in np.arange(-1000, 2000, 100):
        overlap = np.dot(init_moyal_convolution(x, f_histogram, loc, 20), fhistogram)
        if overlap > max_overlap:
            loc_init = loc
            max_overlap = overlap
    #while np.dot(init_moyal_convolution(x, f_histogram, loc_init, 20), fhistogram) < 1:
    #    loc_init += 100
    params.add('loc', value = loc_init, min = -500 + loc_init, max = 500 + loc_init)
    params.add('scale', value = 30, min = 0.5, max = 100)
    #params.add('c', value = np.sum(fhistogram)/np.sum(ihistogram), min = 0)
    nfhistogram = np.array(fhistogram)
    for e in [loc_init]:
        init = init_moyal_convolution(x, f_histogram, e, 20)
        #plt.errorbar(x, ihistogram/np.sum(ihistogram), yerr = np.sqrt(ihistogram) / np.sum(ihistogram), color = 'r', label = "Initial histogram", ls = 'none')
        #plt.errorbar(x, fhistogram/np.sum(fhistogram), yerr = np.sqrt(fhistogram) / np.sum(fhistogram), color = 'b', label = "Final (convolved) histogram", ls = 'none')
        #plt.plot(x, init, color = 'cyan', label = "first guess")
        #plt.legend()
        #plt.show()
    result = model.fit(nfhistogram / np.sum(nfhistogram), x = x, params = params, weights = np.sum(nfhistogram) / np.where(nfhistogram > 0, np.sqrt(nfhistogram), 1e5))
    #result = model.fit(nfhistogram / np.sum(nfhistogram), x = x, params = params, weights = np.sum(nfhistogram) / np.where(nfhistogram > 0, np.sqrt(nfhistogram), 1))
    print(lmfit.fit_report(result))
    rfunc = r_moyal_convolution(f_histogram, result)
    if plot:
        plt.errorbar(x, ihistogram/np.sum(ihistogram), yerr = np.sqrt(ihistogram) / np.sum(ihistogram), color = 'r', label = "Initial histogram", ls = 'none')
        plt.errorbar(x, fhistogram/np.sum(fhistogram), yerr = np.sqrt(fhistogram) / np.sum(fhistogram), color = 'b', label = "Final (convolved) histogram", ls = 'none')
        plt.plot(x, rfunc(x), color = 'magenta', label = "Convolved fit. Reduced χ² = %.3f" % result.redchi)
        plt.plot(x, init, color = 'cyan', label = "first guess")
        plt.legend()
        plt.xlabel("Energy channel")
        plt.ylabel("Normalized counts")
        if foil is not None:
            plt.title("Fitting moyal convolution from empty histogram to histogram of "+foil)
        plt.show()
    rmoyal = lambda x : moyal.pdf(x, result.params['loc'].value, result.params['scale'].value)
    if include_params:
        return rmoyal, Result(result.params['loc'].value, result.params['loc'].stderr), Result(result.params['scale'].value, result.params['scale'].stderr)
    return rmoyal



## Crystal ball stuff ##

A = lambda n, α : (n / abs(α))**n * math.exp(-abs(α)**2/2)
B = lambda n, α : n / abs(α) - abs(α)
N = lambda n, α, σ : 1 / (σ * (C(n, α) + D(n, α)))
C = lambda n, α : n / abs(α) * 1 / (n - 1) * math.exp(-abs(α)**2/2)
D = lambda n, α : math.sqrt(math.pi / 2) * (1 + math.erf(abs(α) / math.sqrt(2)))

def crystal_ball(x, α, n, xbar, σ):
    return np.where((x - xbar)/σ > -α, np.exp(-(x - xbar)**2/(2*σ**2)), A(n, α) * (B(n, α) - (x - xbar)/σ) ** -n)
    
def normalized_crystal_ball(x, α, n, xbar, σ):
    return N(n, α, σ) * crystal_ball(x, α, n, xbar, σ)

def ncb_plus_const(x, α, n, xbar, σ, c):
    #print ("-----")
    #print (α)
    #print (n)
    #print (xbar)
    #print (σ)
    #print (c)
    return c + normalized_crystal_ball(x, α, n, xbar, σ) * (1 - c*len(x))

def soft_cutoff(x, xbar, σ, c, z_cutoff):
    z = (x - xbar) / σ
    return np.where(z > z_cutoff, c*np.exp(-(z - z_cutoff)), c)

def sigmoid_soft_cutoff(x, xbar, σ, c, b):
    z = (x - xbar) / σ
    return c*np.exp(-z)/(1+z) + b

def ncb_plus_soft_const(x, α, n, xbar, σ, c, z_cutoff):
    cutoff_integral = np.dot(np.ones_like(x), soft_cutoff(x, xbar, σ, c, z_cutoff))
    return soft_cutoff(x, xbar,  σ, c, z_cutoff) + normalized_crystal_ball(x, α, n, xbar, σ) * (1 - cutoff_integral)

this_crystal_ball = lambda points : lambda x, α, n, xbar, σ, c : ncb_plus_const(x, α, n, xbar, σ, c, points)

def crystal_ball_params(x, y):
    xbar = x[np.argmax(y)]
    n = 2
    α = 1.5
    above_half = y > max(y)/2
    xmin = min(x[above_half])
    xmax = max(x[above_half])
    σ = (xmax - xmin) / 2 / np.sqrt(2 * np.log(2)) # sigma from hwhm
    params = lmfit.Parameters()
    params.add('xbar', value = xbar, min = xmin, max = xmax)
    params.add('α', value = α, min = 1e-8)
    params.add('n', value = n, min = 1 + 1e-8, max = 50)
    params.add('σ', value = σ, min = 1e-8)
    return params

def result_ball(params):
    return lambda x : normalized_crystal_ball(x, params['α'].value, params['n'].value, params['xbar'].value, params['σ'].value)

def constant_result_ball(params, soft = False):
    if not soft:
        return lambda x : ncb_plus_const(x, params['α'].value, params['n'].value, params['xbar'].value, params['σ'].value, params['c'].value)
    else:
        return lambda x : ncb_plus_soft_const(x, params['α'].value, params['n'].value, params['xbar'].value, params['σ'].value, params['c'].value, params['z_cutoff'].value)

def fit_histogram(histogram, const = False, soft = False, plot = False, plot_initial = False, foil = None):

    npch = np.array(range(len(histogram)))
    nphist = np.array(histogram)
    errs = np.sqrt(histogram)
    nonzero = nphist > 0
    x = npch[nonzero]
    y = nphist[nonzero] / np.sum(nphist)
    weights = 1/(errs / np.sum(nphist))[nonzero]
    params = crystal_ball_params(x, y)
    if const:
        params.add('c', value = min(y), min = 0, max = 50 / np.sum(histogram))
        if soft:
            params.add('z_cutoff', value = -1, min = -5, max = 5)
            model = lmfit.Model(ncb_plus_soft_const)
        else:
            model = lmfit.Model(ncb_plus_const)
    else:
        model = lmfit.Model(normalized_crystal_ball)
    if plot_initial:
        plt.errorbar(x, y, yerr = 1/weights, color = 'b', label = "Normalized energy channel histogram", ls = 'none', marker = '.')
        plt.plot(npch, result_ball(params)(npch), color = 'r', label = "Normalized crystal ball guess")
        plt.axvline(x = params['xbar'].value, label = "Mean", color = 'red', ls = '--')
        if foil is not None:
            plt.title("Initial guess for "+foil)
        plt.xlabel("Channel #")
        plt.ylabel("Normalized counts")
        plt.show()

    result = model.fit(y, x = x, params = params, weights = weights)

    #print (lmfit.fit_report(result))
    
    if plot:
        text_x = 0 if result.params['xbar'] > 800 else 1200
        if const:
            if not soft:
                plt.axhline(result.params['c'].value, color = 'cyan', label = "Background level", ls = '--')
            else:
                plt.plot(x, soft_cutoff(x, result.params['xbar'].value, result.params['σ'].value, result.params['c'].value, result.params['z_cutoff'].value), color = 'cyan', label = "Background level", ls = '--')
                #plt.show()
            plt.text(text_x, np.max(y) / 2 - np.max(y) / 10, "BG level = %.3fx10⁻⁵" % (result.params['c'].value*1e5)) 
        plt.plot(npch, result_ball(result.params)(npch), color = 'r', label = "Crystal ball fit")
        plt.errorbar(x, y, yerr = 1/weights, color = 'black', label = "Energy channel histogram", ls = 'none', marker = '.')
        plt.text(text_x, np.max(y) / 2, "Reduced χ² = %.3f" % result.redchi)
        plt.axvline(x = result.params['xbar'].value, label = "Estimate mean", color = 'red', ls = '--')
        if result.params['xbar'].stderr is not None:
            plt.axvline(x = result.params['xbar'].value + result.params['xbar'].stderr, label = "Estimate error bounds", color = 'orange', ls = '--')
            plt.axvline(x = result.params['xbar'].value - result.params['xbar'].stderr, color = 'orange', ls = '--')
         
        if foil is not None:
            mtitle = "crystal Ball"
            if const:
                mtitle += " plus"
                if soft:
                    mtitle += " soft"
                mtitle += " constant"
            plt.title("Energy channel fit for "+foil+". Model: "+mtitle)
        plt.xlabel("Channel #")
        plt.ylabel("Normalized counts")
        plt.legend()
        plt.show()
    return result

def optimal_energy_fit(histogram):
    
    results = {}
    minresult = None
    minchi = 0
    for const in (False, True):
        for soft in (False, True):
            if (const, soft) != (False, True):
                results[(const, soft)] = (result:= fit_histogram(histogram, const = const, soft = soft))
                if not soft and result.params['xbar'].stderr is not None and (minchi == 0 or result.redchi < minchi): # do not accept soft cutoff as estimate
                    minresult = result
                    minchi = result.redchi
    sys = 0
    for key, result in results.items():
        if result.redchi < minchi and result.params['xbar'].stderr is not None:
            if abs(result.params['xbar'].value - minresult.params['xbar'].value) > sys:
                sys = abs(result.params['xbar'].value - minresult.params['xbar'].value)

    return Result(minresult.params['xbar'].value, stat = minresult.params['xbar'].stderr, sys = sys)

def optimal_function_fit(histogram):
    
    results = {}
    minresult = None
    minchi = 0
    for const in (False, True):
        for soft in (False, True):
            if (const, soft) != (False, True):
                results[(const, soft)] = (result:= fit_histogram(histogram, const = const, soft = soft))
                if (minchi == 0 or result.redchi < minchi): # now, just seek a best fit for convolution purposes or whatever
                    minresult = (const, soft)
                    minchi = result.redchi

    const, soft = minresult
    if const:
        return constant_result_ball(results[minresult].params, soft = soft)
    return result_ball(results[minresult].params)
    




if __name__ == '__main__':
    padding = [0] * 1000
    from data_loader import *
    data = {}
    recursive_read(data, "data", require = [0], reject = ['titanium'])    
    #foils = ('gold', '2gold', '4gold', 'iron')
    foils = ['gold', '2gold', '4gold']
    x = np.array(range(2048))
    numbers = []
    xs = []
    ys = []
    locs = []
    scales = []
    for foil in foils:
        func, loc, scale = fit_moyal_convolution(data[('empty', 0, 1)][1] + padding, data[(foil, 0, 1)][1] + padding, plot = False, foil = foil, include_params = True)
        numbers.append(1 if foil[0] == 'g' else int(foil[0]))
        xs.append(x)
        ys.append(func(x))
        locs.append(loc)
        scales.append(scale)
    for i, foil in enumerate(foils):
        plt.plot(xs[i], ys[i], label = foil)
    plt.legend()
    plt.show()
    from scattering_helpers import plotting_unpack
    ly, lyerr = plotting_unpack(locs)
    sy, syerr = plotting_unpack(scales)
    plt.title("Locations")
    plt.errorbar(numbers, ly, yerr = lyerr, ls = 'none', marker = 'o')
    plt.show()
    plt.title("Scales")
    plt.errorbar(numbers, sy, yerr = syerr, ls = 'none', marker = 'o')
    plt.show()
    
    
    """x = np.array(range(2048))
    # x = np.arange(0, 2047, 1)
    f_histogram = lambda x : x+1
    c = moyal_convolution(x, f_histogram, 0, 3)
    #print (c)
    print (len(x))
    print (len(c))
    #moyal_convolution(x, f_histogram, 0, 1)
    """
"""
if __name__ == '__main__':
    from data_loader import *
    data = {}
    recursive_read(data, "data", require = [0, 'empty'], reject = ['titanium'])
    #recursive_read(data, "data", require = [0], reject = ['titanium'])
    for metadata, entry in data.items():
        foil = metadata[0]
        if foil == 'empty':
            histogram = [0] * 1000 + entry[1][1000:]
        else:
            histogram = entry[1]
        for const in (True, False):
            if const:
                for soft in (True, False):
                    fit_histogram(histogram, const = const, soft = soft, plot = True, foil = foil)
            else:
                fit_histogram(histogram, const = const, plot = True, foil = foil)"""