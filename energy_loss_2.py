import numpy as np
import matplotlib.pyplot as plt
from data_loader import *
from data_processing import *
import math
from plots import *
import lmfit
from scattering_helpers import plotting_unpack

### COPY ME! ###
"""
α
σ
χ²

"""
### CONSTANTS ###
Z_map = {
    'gold' : 79,
    'iron' : 26,
    'titanium' : 22,
    '2gold' : 79,
    '4gold' : 79,
}
density_map = {
    'gold' : 19.32, # g/cm^3
    'iron' : 7.874,
    'titanium' : 4.5,
    '2gold' : 19.32,
    '4gold' : 19.32,
}
mass_map = {
    'gold' : 196.96657, # g/mol
    'iron' : 55.845,
    #'titanium' : ,
    '2gold' : 196.96657,
    '4gold' : 196.96657,
}
I = {
    'gold' : 790e-6, #MeV
    '2gold' : 790e-6,
    '4gold' : 790e-6,
    'iron' : 286e-6,
    'titanium' : 233e-6
}
number = {
    'gold' : 1,
    'iron' : 1,
    'titanium' : 1,
    '2gold' : 2,
    '4gold' : 3,
}

e_empty = Result((5.486 * 0.86 + 5.443 * 0.127 + 5.391 * 0.014) / (0.86 + 0.127 + 0.014))
c = 299792458 # m/s
m_a = 3.7273794066e3 #MeV
m_e = 0.51099895000 #MeV
K = 0.307075 # MeV / mol * cm^2

### DEFINITIONS ###
v = lambda b : b * c
g = lambda b : 1 / math.sqrt(1-b**2)
W_max = lambda b : (2*m_e*b**2*g(b)**2) / (1+2*g(b)*m_e/m_a + (m_e/m_a)**2)
dEdx = lambda beta, foil : K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * (1/2 * math.log(2*m_e*beta**2*g(beta)**2*W_max(beta)/I[foil]**2) - beta**2) * density_map[foil]
b = lambda E : math.sqrt(1 - m_a**2 / E**2)

### EXTRA HELPERS ###


### INTEGRATION ###
def thiccness_dx(foil, incident, exiting):
    E = incident + m_a
    x = 0
    dt = 1e-14
    dx = dt/2 * v(b(E)) * 100
    Es = [incident]
    xs = [x]
    while E > exiting + m_a:
        beta = b(E)
        #print (beta)
        x += dx # cm
        dE = dEdx(beta, foil) * dx
        E -= dE
        Es.append(E - m_a)
        xs.append(x)
    return xs, Es  

def thiccness(foil, incident, exiting):
    E = incident + m_a
    x = 0
    dt = 1e-14
    #print (b(E))
    while E > exiting + m_a:
        beta = b(E)
        x += v(beta) * 100 * dt # cm
        dE = dEdx(beta, foil) * v(beta) * dt
        E -= dE
    return x

### FITTERS ###
A = lambda n, α : (n / abs(α))**n * math.exp(-abs(α)**2/2)
B = lambda n, α : n / abs(α) - abs(α)
N = lambda n, α, σ : 1 / (σ * (C(n, α) + D(n, α)))
C = lambda n, α : n / abs(α) * 1 / (n - 1) * math.exp(-abs(α)**2/2)
D = lambda n, α : math.sqrt(math.pi / 2) * (1 + math.erf(abs(α) / math.sqrt(2)))

def crystal_ball(x, α, n, xbar, σ):
    return np.where((x - xbar)/σ > -α, np.exp(-(x - xbar)**2/(2*σ**2)), A(n, α) * (B(n, α) - (x - xbar)/σ) ** -n)
    
def normalized_crystal_ball(x, α, n, xbar, σ):
    return N(n, α, σ) * crystal_ball(x, α, n, xbar, σ)

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
    params.add('n', value = n, min = 1 + 1e-8)
    params.add('σ', value = σ, min = 1e-8)
    return params

def result_ball(params):
    return lambda x : normalized_crystal_ball(x, params['α'].value, params['n'].value, params['xbar'].value, params['σ'].value)

def fit_histogram(histogram, plot = False, foil = None):

    npch = np.array(range(len(histogram)))
    nphist = np.array(histogram)
    errs = np.sqrt(nphist)
    nonzero = nphist > 0
    x = npch[nonzero]
    y = nphist[nonzero] / np.sum(nphist)
    weights = 1/(errs / np.sum(nphist))[nonzero]
    model = lmfit.Model(normalized_crystal_ball)
    params = crystal_ball_params(x, y)
    if False:
        plt.errorbar(x, y, yerr = 1/weights, color = 'b', label = "Normalized energy channel histogram", ls = 'none', marker = '.')
        plt.plot(npch, result_ball(params)(npch), color = 'r', label = "Normalized crystal ball guess")
        plt.axvline(x = params['xbar'].value, label = "Mean", color = 'red', ls = '--')
        if foil is not None:
            plt.title("Initial guess for "+foil)
        plt.xlabel("Channel #")
        plt.ylabel("Normalized counts")
        plt.show()

    result = model.fit(y, x = x, params = params, weights = weights)

    #print(lmfit.fit_report(result))

    if plot:
        plt.errorbar(x, y, yerr = 1/weights, color = 'b', label = "Energy channel histogram", ls = 'none', marker = '.')
        plt.plot(npch, result_ball(result.params)(npch), color = 'r', label = "Crystal ball fit")
        plt.text(0, np.max(y) / 2, "Reduced χ² = %.3f" % result.redchi)
        plt.axvline(x = result.params['xbar'].value, label = "Estimate mean", color = 'red', ls = '--')
        if result.params['xbar'].stderr is not None:
            plt.axvline(x = result.params['xbar'].value + result.params['xbar'].stderr, label = "Estimate error bounds", color = 'orange', ls = '--')
            plt.axvline(x = result.params['xbar'].value - result.params['xbar'].stderr, color = 'orange', ls = '--')
        
        if foil is not None:
            plt.title("Energy channel fit for "+foil+" (NORMALIZED)")
        plt.xlabel("Channel #")
        plt.ylabel("Normalized counts")
        plt.legend()
        plt.show()

    return Result(result.params['xbar'].value, stat = result.params['xbar'].stderr)

def test_fits():
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium"])

    for metadata, entry in data.items():
        fit_histogram(entry[1], foil = metadata[0])

### ENERGY STUFF ###

def get_energies(data, verbose = False, plot = False):

    channels = {}
    for (foil, angle, iteration), entry in data.items():
        if foil in channels.keys():
            channels[foil] = weighted_average([channels[foil], fit_histogram(entry[1], foil = foil)])
        else:
            channels[foil] = fit_histogram(entry[1], foil = foil)
    energies = {}
    e = lambda ch : e_empty / channels["empty"] * ch
    for foil, channel in channels.items():
        energies[foil] = e(channel) if foil != "empty" else e_empty
        if verbose:
            report(foil + " energy", energies[foil], "MeV")
    return energies

def get_thickness(energies, verbose = False):

    thickness = {}
    incident = e_empty
    for foil, exiting in energies.items():
        thickness[foil] = AsymmetricResult.asymmetric_evaluate(thiccness, foil, incident, exiting)
        if verbose:
            report(foil + " thickness", thickness[foil], "cm")
    return thickness

def do_gold_thickness(thickness):
    x = []
    thicknii = []
    for foil, thick in thickness.items():
        if "gold" in foil:
            try:
                x.append(int(foil[0]))
            except:
                x.append(1)
            thicknii.append(thick)
    y, y_err = plotting_unpack(thicknii)
    plt.errorbar(x, y, yerr = y_err, label = "estimations", color = "black", ls = 'none', marker = '.')
    estithicc = resultify(weighted_average(thicknii))
    report("Gold thickness per foil", estithicc, "cm")
    estix = np.linspace(min(x) - 0.5, max(x) + 0.5)
    plt.plot(estix, [estithicc.val * esti for esti in estix], label = "Average thickness curve", color = "red", ls = '--')
    plt.xticks(x)
    plt.xlabel("# gold foils")
    plt.ylabel("Thickness")
    plt.title("Thickness of gold foils")
    plt.legend()
    plt.show()

### attempting probability stuff ###

def poly_coeffs(energies):
    E2 = 1/np.array(energies)**2
    c = [E2]
    while len(c) < len(E2):
        c.append(c[-1][:-1] * E2[len(c):])
    return [0] + [np.sum(ci) for ci in c]

#p_tot = lambda energies, N, p : N * np.polyval(poly_coeffs(energies), p)
scatter = lambda cps, empty_cps : (empty_cps - cps) / empty_cps
p_tot = lambda coeffs, N, p : N * np.polyval(coeffs, p)
model_p_tot = lambda coeffs : lambda foil, N, p : p_tot(coeffs[f], N, p) if type(foil) == int else [p_tot(coeffs[f], N, p) for f in foil] 

def fit_poly_gold(gold_data):
    
    p_scatter = []
    coeffs = []
    empty_cps = 0
    for time, histogram in iterationless(gold_data)[('empty', 0)]:
        counts = Result(np.sum(histogram), np.sqrt(np.sum(histogram)))
        seconds = Result(time + 0.5, 1/np.sqrt(12))
        if empty_cps == 0:
            empty_cps = counts/seconds
        else:
            empty_cps = weighted_average(empty_cps, counts/seconds)
    foil_energies = get_energies(gold_data)
    incident = foil_energies['empty'].val
    # temp, I think. Need to fit multiple times for more uncertainties on p and N
    for metadata, (time, histogram) in gold_data.items():
        foil = metadata[0]
        if foil != 'empty':
            # add probability of scattering
            counts = Result(np.sum(histogram), np.sqrt(np.sum(histogram)))
            seconds = Result(time + 0.5, 1/np.sqrt(12))
            p_scatter.append(scatter(counts/seconds, empty_cps))

            # get array of energy coefficients
            exiting = foil_energies[foil].val
            x, energies = thiccness_dx(foil, incident, exiting)
            dx = x[1] - x[0]
            coeffs.append(poly_coeffs(energies))

    #p_scatter_1 = np.dot(p_scatter, [2] + [1] * (len(p_scatter) - 2) + [0]) / len(p_scatter)
    #p_scatter_2 = np.dot(p_scatter, [0] + [1] * (len(p_scatter) - 2) + [2]) / len(p_scatter)
    #print (p_scatter_1)
    #print (p_scatter_2)
    
    """
    model = lmfit.Model(model_p_tot(coeffs))
    params = lmfit.Parameters()
    params.add('p', value = 0.5 * dx, min = 0, max = 1)
    params.add('N', value = 1/energies[0], min = 0)
    result = model.fit([p.val for p in p_scatter], foil = range(len(coeffs)), weights = [1/p.tot for p in p_scatter], params = params)
    print (lmfit.fit_report(result))
    """

#from sympy import symbols, Eq, solve
def solve_simultaneous(func):
    pass #TODO





### Testing ###

def energy_data():
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium"])
    #print (data.keys())
    return data

def gold_data(empty = True):
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium", "iron"] if empty else ["titanium", "iron", "empty"])
    #print (data.keys())
    return data

if __name__ == '__main__':

    data = gold_data()
    #energies = get_energies(data)
    fit_poly_gold(data)
    
    #data = energy_data()
    #energies = get_energies(data, verbose = True)
    #thickness = get_thickness(energies, verbose = True)
    #do_gold_thickness(thickness)

    ## Gold thickness ##
