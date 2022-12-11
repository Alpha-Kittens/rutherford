import numpy as np
import matplotlib.pyplot as plt
from data_loader import *
from data_processing import *
import math
from plots import *
import lmfit
from scattering_helpers import plotting_unpack
from histogram_fitters import *

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
    '3gold' : 79,
}
density_map = {
    'gold' : 19.32, # g/cm^3
    'iron' : 7.874,
    'titanium' : 4.5,
    '2gold' : 19.32,
    '3gold' : 19.32,
}
mass_map = {
    'gold' : 196.96657, # g/mol
    'iron' : 55.845,
    #'titanium' : ,
    '2gold' : 196.96657,
    '3gold' : 196.96657,
}
I = {
    'gold' : 790e-6, #MeV
    '2gold' : 790e-6,
    '3gold' : 790e-6,
    'iron' : 286e-6,
    'titanium' : 233e-6
}
number = {
    'gold' : 1,
    'iron' : 1,
    'titanium' : 1,
    '2gold' : 2,
    '3gold' : 4,
}
b_est = {
    'gold' : 3,
    'iron' : 2,
    #'titanium' : 1,
    '2gold' : 3,
    '3gold' : 3,
}

mode = 2
consts = { # jackson + mccarthy with chi squred is 2. 
    'gold' : (0.477, 0.1385) if mode == 2 else (0.607, 0.175),
    '2gold' : (0.477, 0.1385) if mode == 2 else (0.607, 0.175),
    '3gold' : (0.477, 0.1385) if mode == 2 else (0.607, 0.175),
    'iron' : (0.477, 0.1385),
    
}

c_a = {
    'gold' : Result(0.3310636, sys = 0.00793844),
    '2gold' : Result(0.3310636, sys = 0.00793844),
    '3gold' : Result(0.3310636, sys = 0.00793844),
    #'iron' : Result(0.344490747, sys = 0.01103380) unimplemented
}
c_b = {
    'gold' : Result(0.15372616, sys = 0.00397555),
    '2gold' : Result(0.15372616, sys = 0.00397555),
    '3gold' : Result(0.15372616, sys = 0.00397555)
}


#p = Result(0.09742704033851624, stat = 0.002604316920042038, sys = 0.0024188943207263947)
                                                                # error is actually total. "sys" used for asymmetric evaluate
#N = Result(0.8846554831893968, stat= 0.002250923948972613)

p = Result(0.0870)

#e_empty = Result((5.486 * 0.86 + 5.443 * 0.127 + 5.391 * 0.014) / (0.86 + 0.127 + 0.014))
e_empty = Result(4.85527376, sys = 0.00304907)
c = 299792458 # m/s
m_a = 3.7273794066e3 #MeV
m_e = 0.51099895000 #MeV
K = 0.307075 # MeV / mol * cm^2
a_fs = 1/137

def expected_E_square(foil, use_moyal = False, simplistic = False):
    if foil == 'iron':
        print ("Iron is unimplemented")
        raise Exception
    data = {}
    recursive_read(data, "data", require = [0], reject = ['iron', 'titanium'], condition = lambda metadata : metadata[0] == 'empty' or metadata[0] == foil)
    energies, e_map = get_energies(data, verbose = False, channel_map = True)
    if simplistic:
        return Result(energies['empty'].val**2, stat = (energies['empty'] - energies[foil]).val)
    xs, Es = thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = 1e-6)
    print (Es[0])
    E2_x = []
    E2_base = [Result(0)]*120 + [e_map(ch+120)**-2 for ch in range(2048 + 1000 - 120)]
    #print (E2_base[-1])
    #print (E2_base[5])
    #print (E2_base)
    if use_moyal:
        for x in xs:
            moyal_pdf = moyal_convolution_pdf(data[('empty', 0, 1)][1], x)
            #print (np.sum(moyal_pdf))
            E2_x.append(np.dot(E2_base, moyal_pdf))
            #print (E2_x[-1])
    else:
        E2_x = np.array(Es)**2
    #squares = np.array(Es)**2
    pdf = [poly_x(i+1, Es, p.val) for i in range(len(xs))]
    print (np.sum(pdf))
    #plt.title("pdf of scattering occuring at any point in gold foil")
    #plt.plot(xs, np.array(pdf) / np.sum(pdf), label = "Scattering pdf", color = 'blue')
    ##plt.ylim(0, 1.2*max(pdf)/ np.sum(pdf))
    #plt.xlabel("x")
    #plt.ylabel("Probability of scattering")
    #plt.show()
    
    #for x in xs:
    #    print (moyal_convolution_pdf(data[('empty', 0, 1)][1], x))
    #print (E2_x)
    E = np.dot(E2_x, pdf) / np.sum(pdf)
    var = np.dot(np.array(E2_x)**2, pdf) / np.sum(pdf) - E**2
    print (E)
    print (var)
    #plt.plot(xs, [px.val for px in pdf])
    #plt.axvline(E.val)
    #plt.show()
    #print ("-")
    #print ("incident E:", energies['empty'])
    #print ("exiting E:", energies[foil])
    #print ("Expected E:", np.dot(pdf, Es)/np.sum(pdf))
    
    return Result(E, stat = var.val**(1/2))


def expected_E_inverse_square(foil):
    if foil == 'iron':
        print ("Iron is unimplemented")
        raise Exception
    data = {}
    recursive_read(data, "data", require = [0], condition = lambda metadata : metadata[0] == 'empty' or metadata[0] == 'foil')
    energies = get_energies(data)
    xs, Es = thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = 1e-6)
    inverse_square = 1/np.array(Es)**2
    pdf = [poly_x(i+1, Es, p) for i in range(len(xs))]
    E = np.dot(pdf, inverse_square) / np.sum(pdf)
    var = np.dot(pdf, inverse_square**2) / np.sum(pdf) - E**2
    E.stat = (var **(1/2)).val
    #plt.plot(xs, [px.val for px in pdf])
    #plt.axvline(E.val)
    #plt.show()
    print ("-")
    print ("incident E:", energies['empty'])
    print ("exiting E:", energies[foil])
    print ("Expected E:", np.dot(pdf, Es)/np.sum(pdf))
    
    return E



"""
OLD
### DEFINITIONS ###
v = lambda b : b * c
g = lambda b : 1 / math.sqrt(1-b**2)
W_max = lambda b : (2*m_e*b**2*g(b)**2) / (1+2*g(b)*m_e/m_a + (m_e/m_a)**2)
dEdx = lambda beta, foil, c_a, c_b : K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * loss(beta, foil, c_a, c_b) * density_map[foil]
dEdx_simple = lambda beta, foil : K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * L_a(beta, foil) * density_map[foil]
loss = lambda beta, foil, c_a, c_b : L_a(beta, foil) - C_Z(beta, foil, c_a, c_b) + 2*L1(beta, foil) + 2**2 * L2(beta)
L_a = lambda beta, foil : (1/2 * math.log(2*m_e*beta**2*g(beta)**2*W_max(beta)/I[foil]**2) - beta**2)
x = lambda beta, foil : beta**2/(a_fs**2 * Z_map[foil])
#C_Z = lambda beta, foil : 2*(Z_map[foil]**0.4 / (2*x(beta, foil)) + (b_est[foil] + 2) / x(beta, foil)**2 * Z_map[foil])
C_Z = lambda beta, foil, c_a, c_b : c_a * math.log(x(beta, foil)) + c_b
V = lambda beta, foil : beta*c / math.sqrt(Z_map[foil]* a_fs * c)
L1 = lambda beta, foil : L_a(beta, foil) / V(beta, foil)**2 * (consts[foil][0] - consts[foil][1] * math.log(V(beta, foil) + 2)) / math.sqrt(Z_map[foil])
y = lambda beta : 2 * a_fs / beta
L2 = lambda beta : -y(beta)**2 * np.dot([1/(l+1) for l in range(100)], [1/((l+1)**2 + y(beta)**2) for l in range(100)])
b = lambda E : math.sqrt(1 - m_a**2 / E**2)
"""

### DEFINITIONS ###
v = lambda b : b * c
g = lambda b : 1 / math.sqrt(1-b**2)

b = lambda E : math.sqrt(1 - m_a**2 / E**2)
x = lambda beta, foil : beta**2/(a_fs**2 * Z_map[foil])
y = lambda beta : 2 * a_fs / beta
V = lambda beta, foil : beta / (math.sqrt(Z_map[foil]) * a_fs)

dEdx = lambda beta, foil, c_a, c_b : K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * loss(beta, foil, c_a, c_b) * density_map[foil]

W_max = lambda b : (2*m_e*b**2*g(b)**2) / (1+2*g(b)*m_e/m_a + (m_e/m_a)**2)

#loss = lambda beta, foil, c_a, c_b : L_a(beta, foil) - C_Z(beta, foil, c_a, c_b) + 2*L1(beta, foil, c_a, c_b) + 2**2 * L2(beta)
loss = lambda beta, foil, c_a, c_b : L_a(beta, foil) - C_Z(beta, foil, c_a, c_b) + 2*L1(beta, foil, c_a, c_b) + L2(beta)
L_a = lambda beta, foil : (1/2 * math.log(2*m_e*beta**2*g(beta)**2*W_max(beta)/I[foil]**2) - beta**2)
C_Z = lambda beta, foil, c_a, c_b : c_a * math.log(x(beta, foil)) + c_b
L1 = lambda beta, foil, c_a, c_b : (L_a(beta, foil) - C_Z(beta, foil, c_a, c_b)) / V(beta, foil)**2 * (consts[foil][0] - consts[foil][1] * math.log(V(beta, foil) + 2)) / math.sqrt(Z_map[foil])
L2 = lambda beta : -y(beta)**2 * np.dot([1/(l+1) for l in range(100)], [1/((l+1)**2 + y(beta)**2) for l in range(100)])


dEdx_simple = lambda beta, foil : K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * L_a(beta, foil) * density_map[foil]
### EXTRA HELPERS ###

def energy_data():
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium"])
    #print (data.keys())
    return data

def read_gold_data(empty = True):
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium", "iron"] if empty else ["titanium", "iron", "empty"])
    #print (data.keys())
    return data

def all_data():
    data = {}
    recursive_read(data, "data", reject = ["titanium", "iron"], condition = lambda metadata : abs(metadata[1]) <= 10)
    return data

def get_cps(data, time_err = 0.5):
    cps = {}
    for metadata, (time, histogram) in data.items():
        cps[metadata[0]] = Result(np.sum(histogram), stat = np.sqrt(np.sum(histogram))) / Result(time, stat = time_err / math.sqrt(12))
    return cps

### INTEGRATION ###

def thiccness_float(foil, incident, exiting, c_a, c_b):
    if incident <= exiting:
        raise ValueError
    E = incident + m_a
    x = 0
    dx = 1e-16 * v(b(E)) * 100
    #steps = 0
    while E > exiting + m_a:
        beta = b(E)
        x += dx
        E -= dEdx(beta, foil, c_a, c_b) * dx
        #steps += 1
        #if foil == '3gold' and steps % 2000 == 0:
        #    print (steps)
        #    print (E)
        #    print (beta)
        #    print (dEdx(beta, foil, c_a, c_b) * dx)
    return x

def thiccness_dx(foil, incident, exiting, c_a, c_b, dx = None):
    if incident <= exiting:
        raise ValueError
    E = incident + m_a
    x = 0
    xs, Es = [x], [E - m_a]
    if dx == None:
        dx = 1e-16 * v(b(E)) * 100
    while E > exiting + m_a:
        beta = b(E)
        x += dx
        E -= dEdx(beta, foil, c_a, c_b) * dx
        xs.append(x)
        Es.append(E - m_a)
    return xs, Es


### ENERGY STUFF ###

def get_energies(data, verbose = False, plot = False, channel_map = False):

    channels = {}
    for (foil, angle, iteration), entry in data.items():
        #print (foil, angle, iteration)
        if foil in channels.keys():
            #print ("averaging")
            channels[foil] = resultify(weighted_average([channels[foil], e := optimal_energy_fit(entry[1], plot = plot)]))
            #print (e.report())
        else:
            channels[foil] = optimal_energy_fit(entry[1], plot = plot)
        #print(channels[foil].report())
    energies = {}
    e = lambda ch : e_empty / channels["empty"] * ch
    for foil, channel in channels.items():
        energies[foil] = e(channel) if foil != "empty" else e_empty
        if verbose:
            report(foil + " energy", energies[foil], "MeV")
    if channel_map:
        return energies, e
    return energies

def get_thickness(energies, verbose = False):
    # iron thickness is not working yet--would need to estimate shift. 
    thickness = {}
    incident = e_empty
    for foil, exiting in energies.items():
        if foil != 'iron' and foil != 'empty':
            thickness[foil] = AsymmetricResult.asymmetric_evaluate(thiccness_float, foil, incident, exiting, c_a[foil], c_b[foil])
            # consider bootstrapping
            if verbose:
                report(foil + " thickness", thickness[foil], "cm")
    return thickness

def do_gold_thickness(thickness):
    x = []
    thicknii = []
    for foil, thick in thickness.items():
        if "gold" in foil:
            try:
                #x.append(int(foil[0]))
                if int(foil[0]) == 4:
                    x.append(3)
                else:
                    x.append(int(foil[0]))
            except:
                x.append(1)
            thicknii.append(thick)
    y, y_err = plotting_unpack(thicknii)
    plt.errorbar(x, y, yerr = y_err, label = "estimations", color = "black", ls = 'none', marker = '.')
    estithicc = resultify(weighted_average(thicknii))
    report("Gold thickness per foil by weighted average:", estithicc, "cm")
    estix = np.linspace(min(x) - 0.5, max(x) + 0.5)
    model = lmfit.Model(lambda x, a, b : a*x + b)
    result = model.fit(y, x=x, a=estithicc.val, b = 0, weights = 1/np.array(y_err))
    plt.plot(estix, [estithicc.val * esti for esti in estix], label = "Average thickness curve", color = "orange", ls = '--')
    plt.plot(estix, result.params['a'].value * estix + result.params['b'].value, label = "Fit curve", color = "red", ls= '--')
    plt.text(0.5, 0.0007, "Reduced χ²: %.3f" % result.redchi)
    report ("Gold thickness per foil by lmfit", Result(result.params['a'].value, tot = result.params['a'].stderr), "cm")
    report ("Offset by lmfit", Result(result.params['b'].value, tot = result.params['b'].stderr), "cm")
    plt.xticks(x)
    plt.xlabel("# gold foils")
    plt.ylabel("Thickness")
    plt.title("Thickness of gold foils")
    plt.legend(loc = "upper left")
    plt.show()

    return thicknii, estithicc

### attempting probability stuff ###

def plot_polys_v2(e1, e2, ps1, ps2):
    pvals = np.linspace(0, 1, 100)
    pevals1 = np.array([evaluate_poly(e1, p) for p in pvals])
    pevals2 = np.array([evaluate_poly(e2, p) for p in pvals])
    plt.plot(pvals, pevals1, label = "p1")
    plt.plot(pvals, pevals2, label = "p2")
    plt.legend()
    plt.show()
    plt.plot(pvals, pevals2/pevals1, label = "ratio")
    plt.axhline(ps2 / ps1, ls = '--', label = "target ratio")
    plt.legend()
    plt.show()

p_scatter = lambda empty_cps, cps : (empty_cps - cps) / empty_cps

def evaluate_poly(energies, p, resultit = True):
    evaluation = 0 if not resultit else Result(0)
    pp = 1 if not resultit else Result(1)
    for e in energies:
        evaluation += pp * p/e**2
        pp = (p/e**2*-1 + 1) * pp
    #print (evaluation)
    return evaluation

def poly_x(index, energies, p):
    evaluation = 0
    pp = 1
    for e in energies[:index]:
        ret = p/e**2 * pp
        evaluation = ret + evaluation
        pp = (p/e**2*-1 + 1) * pp
    return ret

def poly_bisection(f, g, ratio):
    #print (ratio)
    #print ("beginning search")
    pmin = 0
    pmax = 1
    pguess = (pmax + pmin)/2
    threshold = 1e-8
    sign = 1 if ratio > 1 else -1
    rguess =  g(pguess) / f(pguess)
    while abs((rguess) - ratio) > threshold:
        if sign*rguess > sign*ratio:
            pmin = pguess
        else:
            pmax = pguess
        #print (rguess)
        #print (pguess)
        pguess = (pmax + pmin)/2
        rguess =  g(pguess) / f(pguess)
        if pguess < 1e-5:
            raise Exception
    return pguess

def single_poly_bisection(f, p_scatter):
    #print (ratio)
    #print ("beginning search")
    pmin = 0
    pmax = 1
    pguess = (pmax + pmin)/2
    threshold = 1e-8
    rguess =  f(pguess)
    #print (pguess)
    #print (rguess)
    while abs((rguess) - p_scatter) > threshold:
        if rguess < p_scatter:
            pmin = pguess
        else:
            pmax = pguess
        #print (rguess)
        #print (pguess)
        pguess = (pmax + pmin)/2
        rguess =  f(pguess)
        #print (pguess)
        #print (rguess)
    return pguess

def do_poly_fit(incident, exiting1, exiting2, c_a, c_b, ps1, ps2):
    #print ("Now doing: ", incident, exiting1, exiting2, c_a, c_b, ps1, ps2)
    foil = 'gold'
    E1 = thiccness_dx(foil, incident, exiting1, c_a, c_b, dx = 1e-6)[1]
    E2 = thiccness_dx(foil, incident, exiting2, c_a, c_b, dx = 1e-6)[1]
    p1 = lambda p : evaluate_poly(E1, p, resultit = False)
    p2 = lambda p : evaluate_poly(E2, p, resultit = False)
    #plot_polys_v2(E1, E2, ps1, ps2)
    return poly_bisection(p1, p2, ps2 / ps1)

def calculate_p():
    data = read_gold_data()
    
    cps = get_cps(data)
    empty_cps = cps['empty']
    ps = {}
    for metadata, cps_i in cps.items():
        if metadata != 'empty':
            ps[metadata] = p_scatter(empty_cps, cps_i)
    for foil1 in ps.keys():
        print (foil1, ps[foil1])
        for foil2 in ps.keys():
            if foil1 != foil2:
                print (foil1, foil2, ps[foil1]/ps[foil2])
    energies = get_energies(data, verbose = False)
    thicknii = get_thickness(energies)
    p_tot = None
    E1 = {}
    p1 = {}
    p = {}
    for foil1 in ('gold', '2gold', '3gold'):
    #foil1 = 'gold'
        E1[foil1] = thiccness_dx(foil1, e_empty.val, energies[foil1].val, c_a[foil1].val, c_b[foil1].val, dx = 1e-6)[1]
        p1[foil1] = lambda p : evaluate_poly(E1[foil1], p, resultit = False)
        #print (E1[foil1])
        p[foil1] = single_poly_bisection(p1[foil1], ps[foil1].val)
        print (foil1+": %.4f" % p[foil1])
    for foil1, foil2 in (('gold', '2gold'), ('2gold', 'gold')):
        print (evaluate_poly(E1[foil2], p[foil1], resultit = False))
        print (ps[foil2])
    #print( do_poly_fit(energies['empty'].val, energies[foil1].val, energies[foil2].val, c_a[foil1].val, c_b[foil1].val, ps[foil1].val, ps[foil2].val))
    
    for i, foil1 in enumerate(energies.keys()):
        for j in range(i, len(energies.keys())):
            foil2 = list(energies.keys())[j]
            if foil2 != foil1 and foil2 != 'empty' and foil1 != 'empty': #and (foil1, foil2) != ('2gold', '3gold'):
                #print (foil1, "x", foil2)
                #p = AsymmetricResult.asymmetric_evaluate(do_poly_fit, energies['empty'], energies[foil1], energies[foil2], c_a[foil1], c_b[foil1], ps[foil1], ps[foil2], progress = False)
                p = Result(do_poly_fit(energies['empty'].val, energies[foil1].val, energies[foil2].val, c_a[foil1].val, c_b[foil1].val, ps[foil1].val, ps[foil2].val))
                if p_tot is None:
                    p_tot = p
                else:
                    p_tot = resultify(weighted_average((p_tot, p)))
                print (foil1, "and", foil2+": p = %.4f" % p.val)
                E1 = thiccness_dx(foil1, energies['empty'].val, energies[foil1].val, c_a[foil1].val, c_b[foil1].val, dx = 1e-6)[1]
                N1 = ps[foil1] / evaluate_poly(E1, p, resultit = True)
                E2 = thiccness_dx(foil2, energies['empty'].val, energies[foil2].val, c_a[foil2].val, c_b[foil2].val, dx = 1e-6)[1]
                N2 = ps[foil2] / evaluate_poly(E2, p, resultit = True)
                print (foil1, N1)
                print (foil2, N2)
    return
    """
    print ("Weighted average of p:", p_tot)"""
    foil1 = 'gold'
    foil2 = '2gold'
    for key, val in ps.items():
        print (key)
        print (val)
    #print (do_poly_fit(energies['empty'].val, energies['gold'].val, 1/2*(energies['gold'].val + energies['2gold'].val), c_a['gold'].val, c_b['gold'].val, ps['gold'].val, 1/2*(ps['gold'].val + ps['2gold'].val)))
    #print( do_poly_fit(energies['empty'].val, energies[foil1].val, energies[foil2].val, c_a[foil1].val, c_b[foil1].val, ps[foil1].val, ps[foil2].val))
    #p = AsymmetricResult.asymmetric_evaluate(do_poly_fit, energies['empty'], energies[foil1], energies[foil2], c_a[foil1], c_b[foil1], ps[foil1], ps[foil2], progress = False)
    p = Result(0.09742704033851624, stat = 0.002604316920042038, sys = 0.0024188943207263947)
    #print (p)
    def get_N(p):
        print (p)
        foils = ['gold', '2gold']#, '3gold']
        Ns = []
        for foil in foils:
            E = thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = 1e-6)[1]
            Ns.append(ps[foil] / evaluate_poly(E, p, resultit = True))
        wa = resultify(weighted_average(Ns))
        #print (wa)
        return wa
    print (get_N(p))
    tdx = lambda energy, c_a, c_b : thiccness_dx('3gold', e_empty.val, energy, c_a, c_b, dx = 1e-6)[-1]
    super_evaluate = lambda energy, c_a, c_b, p : evaluate_poly(tdx(energy, c_a, c_b), p, resultit = False)
    print(get_N(p) * AsymmetricResult.asymmetric_evaluate(super_evaluate, energies['3gold'], c_a['3gold'], c_b['3gold'], p, progress = True))
    print(ps['3gold'])

                

    #N = AsymmetricResult.asymmetric_evaluate(get_N, p)
    #print (N)






### PLOTS ###



def plot_dEdx():
    b0 = b(2 + m_a)
    betas = np.linspace(b0, 0.06)
    foil = 'gold'
    const = K * 2**2 * Z_map[foil] / mass_map[foil] * density_map[foil] 
    plt.plot(betas, [const/beta**2 * L_a(beta, foil) for beta in betas], label = "Bethe")
    plt.plot(betas, [dEdx(beta, foil, c_a[foil].val, c_b[foil].val) for beta in betas], color = 'blue', label = "total")
    plt.plot(betas, [const/beta**2 * C_Z(beta, foil, c_a[foil].val, c_b[foil].val) for beta in betas], color = 'orange', label = "C/Z")
    plt.plot(betas, [const/beta**2 * L1(beta, foil, c_a[foil].val, c_b[foil].val) for beta in betas], label = "L1")
    plt.plot(betas, [const/beta**2 * L2(beta) for beta in betas], label = "L2")
    plt.legend()
    plt.axhline(0, color = 'r', ls = '--')
    plt.show()



def poly_plot_test(dx = 1e-5):
    E1 = np.arange(5, 4, -0.1)
    E2 = np.arange(5, 3, -0.1)
    
    data = energy_data()
    energies = get_energies(data, verbose = False)
    empty_cps = (np.sum(data[('empty', 0, 1)][1]) / data[('empty', 0, 1)][0])
    foil = 'gold'
    E1 = thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = dx)[1]
    ps1 = (empty_cps - (np.sum(data[(foil, 0, 1)][1]) / data[(foil, 0, 1)][0]))/empty_cps
    foil = '2gold'
    E2 = thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = dx)[1]
    ps2 = (empty_cps - (np.sum(data[(foil, 0, 1)][1]) / data[(foil, 0, 1)][0]))/empty_cps
    plot_polys_v2(E1, E2, ps1, ps2)

def gold_data_2():
    data = {}
    recursive_read(data, "data", require = [0], reject = ["titanium", "iron"], condition = lambda metadata : metadata[2] != 3 and ("gold" in metadata[0] or "empty" in metadata[0]))
    return data

if __name__ == '__main__':
    #data = all_data()
    data = {}
    #recursive_read(data, "data", require = ["2gold", 0], reject = ["titanium", "iron"], condition = lambda metadata : abs(metadata[1]) <= 10)
    recursive_read(data, "data", require = [0], reject = ["titanium", "iron", "3gold"], condition = lambda metadata : metadata[2] != 3 and ("gold" in metadata[0] or "empty" in metadata[0]))
    #print (data.keys())
    #recursive_read(data, "data", require = ["empty", 0])
    #print (data.keys())
    cps = get_cps(data)
    for metadata, (time, histogram) in data.items():
        print (metadata, time, np.sum(histogram), Result(np.sum(histogram), stat = np.sqrt(np.sum(histogram))) / Result(time, stat = 0.5 / math.sqrt(12)))
    energies = get_energies(data, verbose = True, plot = False)
    thiccness = get_thickness(energies, verbose = True)
    do_gold_thickness(thiccness)
    #print (data.keys())
    #for metadata, val in iterationless(data).items():
    #    print (metadata, ":", len(val))
    #expected_E_square('gold', use_moyal = True)
    #calculate_p()
    #data = read_gold_data()
    #energies = get_energies(data)
    #thickness = get_thickness(energies, verbose = True)
    #do_gold_thickness(thickness)   
    """
    for foil in ('gold', '2gold', '3gold'):
        print ("Foil:", foil)
        print ("With moyals:")
        print (expected_E_square(foil))
        print ("Without:")
        print (expected_E_square(foil, use_moyal = False))
        print ("-")
    """





#p_tot = lambda energies, N, p : N * np.polyval(poly_coeffs(energies), p)
scatter = lambda cps, empty_cps : (empty_cps - cps) / empty_cps
p_tot = lambda coeffs, N, p : N * np.polyval(coeffs, p)
model_p_tot = lambda coeffs : lambda foil, N, p : p_tot(coeffs[f], N, p) if type(foil) == int else [p_tot(coeffs[f], N, p) for f in foil] 

### Testing ###


def gold_thickness():
    data = energy_data()
    energies = get_energies(data, verbose = True)
    thickness = get_thickness(energies, verbose = True)
    return do_gold_thickness(thickness)
