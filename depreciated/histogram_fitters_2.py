import lmfit
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import moyal
from data_processing import *


#### Fit results ####
"""a_loc = Result(-1031375.697618693, stat = 26240.65526430214)
b_loc = Result(1507.6611799242887, stat = 6.86358751455451)
a_scale = Result(55233.98556503942, stat = 391.05899088989685)
b_scale = Result(0.981144735446772, stat = 0.09993857071754861)

loc = lambda x : a_loc * x + b_loc
scale = lambda x : a_loc * x + b_loc"""
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

loc_quadratic = lambda x : -342518953.9609103*x**2 + -1116484.5605326646*x + 1526.356081390396
scale_quadratic = lambda x : 115775685.27995577*x**2 + 6377.829708029055*x + 14.703687598788186

loc_quadratic = lambda x : -400825937.7634321*x**2 + -1070175.2906392978*x + 1519.0544474418773
scale_quadratic = lambda x : 115775685.27995577*x**2 + 6377.829708029055*x + 14.703687598788186

nmoyal = lambda x, loc, scale : moyal.pdf(-x, -loc, scale)

def make_predictive_plot():
    padding = [0] * 1000
    from data_loader import recursive_read
    data = {}
    recursive_read(data, "data", require = [0], reject = ['titanium'])    
    ihistogram = data[('empty', 0, 1)][1] + padding
    ch = np.array(range(len(ihistogram)))
    x = np.linspace(0, 0.0006, 12)
    Es = []
    mus = []
    plt.plot(ch, ihistogram/np.sum(ihistogram), label = "initial")
    f_histogram = optimal_function_fit(ihistogram)
    plt.plot(ch, f_histogram(ch))
    plt.show()
    from energy_loss_3 import get_energies
    #get_energies(data)
    for xi in x:
        print (xi)
        #if xi != x[2]:
        if True:
            print (xi)
            fhistogram = moyal_convolution_2(ch, f_histogram, loc_quadratic(xi), scale_quadratic(xi))
            plt.plot(ch, fhistogram/np.sum(fhistogram), label = "final")
            plt.title("x = "+str(xi))
            plt.legend()
            #plt.show()
            Ei = optimal_energy_fit(fhistogram * np.sum(ihistogram), plot = False)
            mus.append(np.dot(fhistogram, ch))
            print (Ei)
            Es.append(Ei.val)
    plt.show()
    print (x)
    print (Es)
    print ("mus")
    print (mus)
    plt.plot(x, Es, ls = 'none', marker = 'o')
    plt.xlabel("x position")
    plt.ylabel("Fitted crystal ball mean of moyal prediction")
    plt.show()

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
    test_vals = nmoyal(test_domain, loc = loc, scale = scale)
    #print (test_vals)
    good_domain = test_domain[test_vals > 1e-4]
    #print (len(good_domain))
    if len(good_domain) % 2 == 0:
        good_domain = np.append(good_domain, good_domain[-1]+1)
    #print (len(good_domain))
    convolution = convolve(f_histogram, x, lambda x : nmoyal(x, loc = loc, scale = scale), good_domain)
    padding = int((len(x) - len(convolution[0])) /2)
    return [0 for i in range(padding)] + list(convolution[0]) + [0 for i in range(padding)]

def moyal_convolution_2(x, f_histogram, loc, scale):
    return simple_convolve(f_histogram, x, lambda x : nmoyal(x, loc = loc, scale = scale), x)

def init_moyal_convolution(x, f_histogram, loc, scale, test_bounds = 100):
    """    
    test_domain = loc + np.arange(1*test_bounds, 2*test_bounds, 1)
    print (test_domain)

    test_vals = nmoyal(test_domain, loc = loc, scale = scale) 
    test_vals = moyal.pdf(test_domain, loc = loc, scale = scale)
    plt.plot(test_domain, test_vals, label = 'moyal')
    plt.legend()
    plt.show()
    plt.plot(x, np.convolve(f_histogram(x), nmoyal(x, loc = loc, scale = scale), mode = 'same'), label = 'convolution')
    plt.legend()
    plt.show()
    #print (test_vals)
    good_domain = test_domain[test_vals > 1e-4]
    #print (len(good_domain))
    if len(good_domain) % 2 == 0:
        good_domain = np.append(good_domain, good_domain[-1]+1)"""
    #print (len(good_domain))
    convolution = simple_convolve(f_histogram, x, lambda x : nmoyal(x, loc = loc, scale = scale), x)
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
    for loc in np.arange(-2000, 2000, 100):
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
        plt.plot(x, init, color = 'cyan', label = "first guess")
        plt.plot(x, rfunc(x), color = 'magenta', label = "Convolved fit. Reduced χ² = %.3f" % result.redchi)
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
    #print (α, n, xbar, σ, "-", sep = '\n')
    if False:
    #if xbar <= 1400:
        plt.plot(x, N(n, α, σ) * crystal_ball(x, α, n, xbar, σ))
        plt.show()
    if True in np.isnan(crystal_ball(x, α, n, xbar, σ)):
        plt.plot(x, N(n, α, σ) * crystal_ball(x, α, n, xbar, σ))
        plt.show()
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
    params.add('α', value = α, min = 1e-8, max = 50)
    params.add('n', value = n, min = 1 + 1e-8, max = 50)
    params.add('σ', value = σ, min = 1e-8, max = 200)
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
    print (nphist)
    errs = np.sqrt(histogram)
    nonzero = nphist > 0
    x = npch[nonzero]
    y = nphist[nonzero] / np.sum(nphist)
    print (len(x))
    print (len(y))
    weights = 1/(errs[nonzero] / np.sum(nphist))
    print (x)
    print (y)
    print (weights)
    params = crystal_ball_params(x, y)
    if const:
        params.add('c', value = min(y[x < params['xbar'].value]), min = 0, max = 50 / np.sum(histogram))
        if soft:
            params.add('z_cutoff', value = -1, min = -5, max = 5)
            model = lmfit.Model(ncb_plus_soft_const)
        else:
            model = lmfit.Model(ncb_plus_const)
    else:
        model = lmfit.Model(normalized_crystal_ball)
    if plot_initial:
        plt.errorbar(x, y, yerr = 1/weights, color = 'b', label = "Normalized energy channel histogram", marker = '.')#, ls = 'none')
        plt.plot(npch, result_ball(params)(npch), color = 'r', label = "Normalized crystal ball guess")
        plt.title("dat boi")
        
        plt.axvline(x = params['xbar'].value, label = "Mean", color = 'red', ls = '--')
        if foil is not None:
            plt.title("Initial guess for "+foil)
        plt.xlabel("Channel #")
        plt.ylabel("Normalized counts")
        plt.legend()
        plt.show()

    result = model.fit(y, x = x, params = params, weights = weights)
    if plot_initial:
        print(lmfit.fit_report(result))
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

def optimal_energy_fit(histogram, plot = False):
    
    results = {}
    minresult = None
    minchi = 0
    for const in (False, True):
        for soft in (False, True):
            if (const, soft) != (False, True):
                print (const, soft)
                result = fit_histogram(histogram, const = const, soft = soft, plot = plot, plot_initial = plot)
                results[(const, soft)] = result
                if not soft and result.params['xbar'].stderr is not None and (minchi == 0 or result.redchi < minchi): # do not accept soft cutoff as estimate
                    minresult = result
                    minchi = result.redchi
    sys = 0
    for key, result in results.items():
        #print (result.redchi)
        #print (result.params['xbar'].stderr)
        #print (minchi)
        if minchi == 0 or result.redchi < minchi and result.params['xbar'].stderr is not None:
            if abs(result.params['xbar'].value - minresult.params['xbar'].value) > sys:
                sys = abs(result.params['xbar'].value - minresult.params['xbar'].value)
        #print (minchi)

    return Result(minresult.params['xbar'].value, stat = minresult.params['xbar'].stderr, sys = sys)

def optimal_function_fit(histogram):
    
    results = {}
    minresult = None
    minchi = 0
    for const in (False, True):
        for soft in (False, True):
            if (const, soft) != (False, True):
                print (const, soft)
                result = fit_histogram(histogram, const = const, soft = soft)
                results[(const, soft)] = result
                if (minchi == 0 or result.redchi < minchi): # now, just seek a best fit for convolution purposes or whatever
                    minresult = (const, soft)
                    minchi = result.redchi

    const, soft = minresult
    if const:
        return constant_result_ball(results[minresult].params, soft = soft)
    return result_ball(results[minresult].params)
    



def fitting_moyals(plot = False, do_quadratic = False):
    padding = [0] * 1000
    from data_loader import recursive_read
    data = {}
    recursive_read(data, "data", require = [0], reject = ['titanium'])    
    #foils = ('gold', '2gold', '4gold', 'iron')
    foils = ['empty', 'gold', '2gold', '4gold']
    x = np.array(range(2048))
    numbers = []
    xs = []
    ys = []
    locs = []
    scales = []
    for foil in foils:
        func, loc, scale = fit_moyal_convolution(data[('empty', 0, 1)][1] + padding, data[(foil, 0, 1)][1] + padding, plot = False, foil = foil, include_params = True)
        numbers.append(0 if foil == 'empty' else 1 if foil[0] == 'g' else int(foil[0]))
        xs.append(x)
        ys.append(func(x))
        locs.append(loc)
        scales.append(scale)
    from scattering_helpers import plotting_unpack
    ly, lyerr = plotting_unpack(locs)
    sy, syerr = plotting_unpack(scales)
    from energy_loss_2 import get_energies, get_thickness
    thicknesses = get_thickness(get_energies(data))
    thicknii = []
    for foil in foils:
        if foil == 'empty':
            thicknii.append(Result(0))
        else:
            thicknii.append(thicknesses[foil])
    thickx, thickxerr = plotting_unpack(thicknii)
    quadratic = lambda x, a, b, c : a*x**2 + b*x + c
    linear = lambda x, a, b : a*x + b
    test_thickx = np.linspace(min(thickx) - min(thickx) / 2, max(thickx) + min(thickx) / 2, 100)
    if do_quadratic:
        model = lmfit.Model(quadratic)
    else:
        model = lmfit.Model(linear)
    slope = lambda yarr, xarr : (yarr[1] - yarr[0]) / (xarr[1] - xarr[0])
    adjusted_lyerr = np.sqrt([lyerr[i]**2 + (slope(ly, thickx) * thickxerr[i])**2 for i in range(len(ly))])
    adjusted_syerr = np.sqrt([syerr[i]**2 + (slope(sy, thickx) * thickxerr[i])**2 for i in range(len(sy))])
    if do_quadratic:
        loc_result = model.fit(ly, a = -10000, b = -10000, c = 0, x=thickx, weights = 1/adjusted_lyerr)
        scale_result = model.fit(sy[1:], a = 0.5, b = 1, c = 0, x=thickx[1:], weights = 1/adjusted_syerr[1:])
    else:
        loc_result = model.fit(ly, a = -10000, b = 0, x=thickx, weights = 1/adjusted_lyerr)
        scale_result = model.fit(sy[1:], a = 1, b = 0, x=thickx[1:], weights = 1/adjusted_syerr[1:])
    print ("Loc")
    print (loc_result.params['a'].value, loc_result.params['a'].stderr)
    a_loc = loc_result.params['a'].value
    print (loc_result.params['b'].value, loc_result.params['b'].stderr)
    b_loc = loc_result.params['b'].value
    if do_quadratic:
        c_loc = loc_result.params['c'].value
        print (loc_result.params['c'].value, loc_result.params['c'].stderr)
    print (loc_result.redchi)
    print ("Scale")
    print (scale_result.params['a'].value, scale_result.params['a'].stderr)
    a_scale = scale_result.params['a'].value
    print (scale_result.params['b'].value, scale_result.params['b'].stderr)
    b_scale = scale_result.params['b'].value
    if do_quadratic:
        c_scale = scale_result.params['c'].value
        print (scale_result.params['c'].value, scale_result.params['c'].stderr)
    print (scale_result.redchi)
    if plot:
        for i, foil in enumerate(foils):
            plt.plot(xs[i], ys[i], label = foil)
        plt.xlabel("Channel number")
        plt.ylabel("Moyal evaluation")
        plt.title("Moyals modeling energy spreading for each foil")
        plt.legend()
        plt.show()
        plt.title("Locations")
        plt.xlabel("Foil Thickness (cm)")
        plt.ylabel("'Loc' of moyal for convolution")
        plt.errorbar(thickx, ly, yerr = lyerr, xerr = thickxerr, ls = 'none', marker = 'o', label = "data")
        if do_quadratic:
            plt.plot(test_thickx, quadratic(test_thickx, a_loc, b_loc, c_loc), label = "Quadratic fit. Redchi: %.3f" % loc_result.redchi)
        else:
            plt.plot(test_thickx, linear(test_thickx, a_loc, b_loc), label = "linear fit. Redchi: %.3f" % loc_result.redchi)
        plt.legend()
        plt.savefig(("Quadratic" if do_quadratic else "Linear") + "fit_for_Location.png")
        plt.show()
        plt.title("Scales")
        plt.xlabel("Foil Thickness (cm)")
        plt.ylabel("'Scale' of moyal for convolution")
        plt.errorbar(thickx, sy, yerr = syerr, xerr = thickxerr, ls = 'none', marker = 'o', label = "data")
        if do_quadratic:
            plt.plot(test_thickx, quadratic(test_thickx, a_scale, b_scale, c_scale), label = "Quadratic fit. Redchi: %.3f" % scale_result.redchi)
        else:
            plt.plot(test_thickx, linear(test_thickx, a_scale, b_scale), label = "linear fit. Redchi: %.3f" % scale_result.redchi)
        plt.legend()
        plt.savefig(("Quadratic" if do_quadratic else "Linear") + "fit_for_Scale.png")
        plt.show()
    #return loc_result


    
    
    
def moyal_convolution_pdf(ihistogram, x, padding = [0] * 1000):
    convolution_histogram = (ihistogram + padding) / np.sum(ihistogram)
    moyal_domain = np.array(range(len(convolution_histogram)))
    return np.convolve(convolution_histogram, moyal.pdf(moyal_domain, loc = loc_quadratic(x).val, scale = scale_quadratic(x).val), mode = 'same')

def plot_inverted_moyal(x, loc = 0, scale = 1):
    y = moyal.pdf(-x, loc = loc, scale = scale)
    plt.plot(x, y)
    plt.show()

if __name__ == '__main__':
    make_predictive_plot()
    #x = np.linspace(-5, 5, 100)
    #plot_inverted_moyal(x)
    #fitting_moyals(plot = False, do_quadratic = True)
    #fitting_moyals(plot = True)
    
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