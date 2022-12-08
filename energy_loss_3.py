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
    '4gold' : 4,
}
b_est = {
    'gold' : 3,
    'iron' : 2,
    #'titanium' : 1,
    '2gold' : 3,
    '4gold' : 3,
}
consts = { # jackson + mccarthy with chi squred is 2. 
    'gold' : (0.477, 0.1385),
    '2gold' : (0.477, 0.1385),
    '4gold' : (0.477, 0.1385),
    'iron' : (0.477, 0.1385),
    
}

c_a = {
    'gold' : Result(0.3310636, sys = 0.00793844),
    '2gold' : Result(0.3310636, sys = 0.00793844),
    '4gold' : Result(0.3310636, sys = 0.00793844),
    #'iron' : Result(0.344490747, sys = 0.01103380) unimplemented
}
c_b = {
    'gold' : Result(0.15372616, sys = 0.00397555),
    '2gold' : Result(0.15372616, sys = 0.00397555),
    '4gold' : Result(0.15372616, sys = 0.00397555)
}


#p = Result(0.09742704033851624, stat = 0.002604316920042038, sys = 0.0024188943207263947)
                                                                # error is actually total. "sys" used for asymmetric evaluate
#N = Result(0.8846554831893968, stat= 0.002250923948972613)

p = Result(0.0870)

e_empty = Result((5.486 * 0.86 + 5.443 * 0.127 + 5.391 * 0.014) / (0.86 + 0.127 + 0.014))
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
        #if foil == '4gold' and steps % 2000 == 0:
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
        if foil in channels.keys():
            channels[foil] = weighted_average([channels[foil], optimal_energy_fit(entry[1])])
        else:
            channels[foil] = optimal_energy_fit(entry[1])
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
    for foil1 in ('gold', '2gold', '4gold'):
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
            if foil2 != foil1 and foil2 != 'empty' and foil1 != 'empty': #and (foil1, foil2) != ('2gold', '4gold'):
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
        foils = ['gold', '2gold']#, '4gold']
        Ns = []
        for foil in foils:
            E = thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = 1e-6)[1]
            Ns.append(ps[foil] / evaluate_poly(E, p, resultit = True))
        wa = resultify(weighted_average(Ns))
        #print (wa)
        return wa
    print (get_N(p))
    tdx = lambda energy, c_a, c_b : thiccness_dx('4gold', e_empty.val, energy, c_a, c_b, dx = 1e-6)[-1]
    super_evaluate = lambda energy, c_a, c_b, p : evaluate_poly(tdx(energy, c_a, c_b), p, resultit = False)
    print(get_N(p) * AsymmetricResult.asymmetric_evaluate(super_evaluate, energies['4gold'], c_a['4gold'], c_b['4gold'], p, progress = True))
    print(ps['4gold'])

                

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


if __name__ == '__main__':
    expected_E_square('gold', use_moyal = True)
    #calculate_p()
    #data = read_gold_data()
    #energies = get_energies(data)
    #thickness = get_thickness(energies, verbose = True)
    #do_gold_thickness(thickness)   
    """
    for foil in ('gold', '2gold', '4gold'):
        print ("Foil:", foil)
        print ("With moyals:")
        print (expected_E_square(foil))
        print ("Without:")
        print (expected_E_square(foil, use_moyal = False))
        print ("-")
    """
def test_solve():
    data = energy_data()
    energies = get_energies(data, verbose = True)
    dx = 100 * 1e-16 * v(b(m_a + energies['empty'].val)) 
    energy_arrays = []
    p_scatters = []
    empty_cps = (np.sum(data[('empty', 0, 1)][1]) / data[('empty', 0, 1)][0])
    for foil in ('gold', '2gold', '4gold'):
        p_scatters.append((empty_cps - (np.sum(data[(foil, 0, 1)][1]) / data[(foil, 0, 1)][0]))/empty_cps)
        energy_arrays.append(thiccness_dx(foil, energies['empty'].val, energies[foil].val, c_a[foil].val, c_b[foil].val, dx = 1e-5)[1])
    solve_probability_distribution(energy_arrays[0:2], p_scatters[0:2], 30, 1e-5)


    
#if __name__ == '__main__':
#    calculate_p()
    #poly_plot_test(dx = 1e-8)
    #test_poly_coeffs()
    #test_solve()
    #plot_dEdx()
    #test_poly_coeffs()
    

    #print (dx)

    #thickness = get_thickness(energies, verbose = True)
    #do_gold_thickness(thickness)    
"""
def poly_coeffs(energies):
    E2 = 1/np.array(energies)**2
    c = [E2]
    while len(c) < len(E2):
        c.append(c[-1][:-1] * E2[len(c):])
    return [0] + [np.sum(ci) for ci in c]
"""

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

    print (gold_thickness())
    dx = gold_thickness().val / 1000
    for metadata, (time, histogram) in gold_data.items():
        foil = metadata[0]
        if foil != 'empty':
            # add probability of scattering
            counts = Result(np.sum(histogram), np.sqrt(np.sum(histogram)))
            seconds = Result(time + 0.5, 1/np.sqrt(12))
            p_scatter.append(scatter(counts/seconds, empty_cps))

            # get array of energy coefficients
            exiting = foil_energies[foil].val
            #x, energies = thiccness_dx(foil, incident, exiting)
            print ("--")
            print (foil)
            #print (incident)
            #print (exiting)
            thickness, energy_func = thiccness_dx(foil, incident, exiting)
            x = np.arange(0, thickness, dx)
            energies = energy_func(x)
            print (x)
            print (energies)
            coeffs.append(poly_coeffs(energies))

    N, positive_roots = cheater_poly(coeffs[0], coeffs[1], p_scatter[0].val, p_scatter[1].val)
    N2, positive_roots2 = cheater_poly(coeffs[2], coeffs[1], p_scatter[2].val, p_scatter[1].val)
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

def cheater_poly(coeffs1, coeffs2, y1, y2):
    c1 = coeffs1[::-1]
    c2 = coeffs2[::-1]
    factor = np.polyval(c1, 0.5) 
    #print ("p = 0.5:", factor)
    c1_prior =np.array([0] * (max(len(coeffs1), len(coeffs2)) - len(coeffs1)) + c1)
    c2_prior = np.array([0] * (max(len(coeffs1), len(coeffs2)) - len(coeffs2)) + c2) 
    #print ("c1_prior:", c1_prior)
    #print ("Prior:", np.polyval(c1_prior, 0.5))
    c1 = c1_prior / factor
    c2 = c2_prior / factor
    #print (coeffs1)
    #print (c1)
    #print (c2)
    print ("scaling factor", factor)

    pvals = np.linspace(0, 1, 100)
    pevals1 = np.polyval(c1, pvals)
    pevals2 = np.polyval(c2, pvals)
    plt.plot(pvals, pevals1)
    plt.plot(pvals, pevals2)
    plt.plot(pvals, pevals2/pevals1)
    plt.axhline(y2 / y1, ls = '--')
    plt.show()
    #print ("After scaling:", np.polyval(c1, 0.5), np.polyval(c2, 0.5))
    #print (np.polyval(c1_prior, 0.5) / factor)
    # coefficients seem to be working so far. 
    """
    new_coeffs = (y2 * c1 - y1 * c2)
    #print (new_coeffs)
    roots = np.roots(new_coeffs)
    print (roots)
    complex_threshold = 1e-8
    real_roots = np.real(roots)[np.abs(roots - np.real(roots)) < complex_threshold]
    print (real_roots)
    #print (real_roots)
    positive_roots = real_roots[real_roots > 0]
    print (positive_roots)
    valid_roots = positive_roots[positive_roots < 1]
    print (valid_roots)
    #print (real_roots)
    #print ("roots:", positive_roots)
    evals = np.polyval(c1, valid_roots), np.polyval(c2, valid_roots)
    #print ("Polynomial evaluations:", evals)
    #print ("check:", y2 * np.polyval(c1, valid_roots) - y1 * np.polyval(c2, valid_roots))
    print (evals)
    N1 = y1 / (evals[0] * factor)
    print (N1)
    N2 = y2 / (evals[1] * factor)
    print (N2)
    if abs(N1 - N2) > 1e-8:
        print ("Something went wrong")
        raise Exception
    """    
    #print (np.polyval(new_coeffs, valid_roots))
    #print (y1 * evals[1])
    #print (y2 * evals[0])
    #print (y2/y1)
    #print (evals[1] / evals[0])
    print ("--")
    
    #print (N2)
    #print ("y1 =",y1)
    #print (N1 * np.polyval(coeffs1, positive_roots))
    #print ("y2 =",y2)
    #print (N2 * np.polyval(coeffs2, positive_roots))
    
    #return N, positive_roots

#from sympy import symbols, Eq, solve
def solve_simultaneous(poly1, y1, poly2, y2):
    p = 1
    N = y1 / poly1(p)
    bigboi = poly2(p)/poly1(p) 





### Testing ###


def gold_thickness():
    data = energy_data()
    energies = get_energies(data, verbose = True)
    thickness = get_thickness(energies, verbose = True)
    return do_gold_thickness(thickness)
"""
if __name__ == '__main__':

    #data = read_gold_data()
    #energies = get_energies(data)
    #fit_poly_gold(data)
    
    data = energy_data()
    energies = get_energies(data, verbose = True)
    thickness = get_thickness(energies, verbose = True)
    do_gold_thickness(thickness)

    ## Gold thickness ##"""


""" 
depreciated
from models import interpolate
def thiccness_dx(foil, incident, exiting, c_a, c_b):
    E = incident + m_a
    x = 0
    dt = 2e-16
    Es = [incident]
    xs = [x]
    dEs = [dEdx_simple(b(E), foil) * v(b(E)) * 100 * dt]
    dE_complexs = [dEdx(b(E), foil, c_a, c_b) * v(b(E)) * 100 * dt]
    Cs = []
    L1s = []
    L2s = []
    i = 0
    while E > exiting + m_a:
        beta = b(E)
        #print (beta)
        dx = v(beta) * 100 * dt # cm
        x += dx
        dE = dEdx_simple(beta, foil) * dx
        dE_complexs.append(dEdx(b(E), foil, c_a, c_b) * v(b(E)) * 100 * dt)
        if i % 50 == 0:
            print (C_Z(beta, foil, c_a, c_b))
            print ("x =",beta**2/(a_fs**2 * Z_map[foil]))
        dEs.append(dE)
        E -= dE
        Es.append(E - m_a)
        i += 1
        Cs.append(C_Z(beta, foil, c_a, c_b))
        L1s.append(L1(beta, foil))
        L2s.append(L2(beta))
        xs.append(x)
        #if i % 50 == 0:
        #    plt.plot(xs, Es, label = "Es")
        #    plt.plot(xs, [d*20 for d in dEs], label = "dEs")
        #    plt.axhline(exiting, color = 'r', ls = '--')
        #    plt.legend()
        #    plt.show()
    plt.plot(xs, Es, label = "Es")
    plt.plot(xs, [d*20 for d in dEs], label = "dEs*20")
    plt.plot(xs, dE_complexs, label = "With 'small' corrections")
    labels = ["Cs", "L1s", "L2s"]
    cols = ["magenta", "orange", "green"]
    for i,y in enumerate((Cs, L1s, L2s)):
        plt.plot(xs[1:], y, label = labels[i], color = cols[i])
    plt.axhline(exiting, color = 'r', ls = '--')
    plt.legend()
    plt.show()
    #print (xs)
    print (x)
    print (len(xs))
    #print (dEs)
    return x, interpolate(Es, xs)[0]
    #return xs, Es  

def thiccness(foil, incident, exiting):
    E = incident + m_a
    x = 0
    dt = 1e-16
    #print (b(E))
    #print (E)
    #print ("---")
    #print (exiting + m_a)
    #print (E)
    if exiting > incident:
        raise ValueError
    i = 0
    print ("--")
    print (exiting)
    while E > exiting + m_a:
        beta = b(E)
        dx = v(beta) * dt * 100 # cm
        #print (dx)
        x += dx
        dE = dEdx(beta, foil) * v(beta) * dx
        E -= dE
        #print (dE)
        #print (E)
        #print ("-")
        #print (E)
        i += 1
        if i == 82 or i == 81 or i == 80:
            print (x)
            print (E)
    print ("Timesteps:", i)
    return x
"""

"""
Moved to histogram_fitters.py
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

def fit_histogram(histogram, min_counts = 1, plot = False, foil = None):

    npch = np.array(range(len(histogram)))
    nphist = np.array(histogram) - min_counts
    errs = np.sqrt(histogram)
    nonzero = nphist >= 0
    x = npch[nonzero]
    y = nphist[nonzero] / np.sum(nphist[nonzero])
    weights = 1/(errs / np.sum(nphist[nonzero]))[nonzero]
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
        print (len(x))
        removed = np.logical_not(nonzero)
        plt.errorbar(x, y, yerr = 1/weights, color = 'black', label = "Signal histogram", ls = 'none', marker = '.')
        plt.errorbar(npch[removed], nphist[removed]/np.sum(nphist), yerr = (errs / np.sum(nphist))[removed], color = 'orange', label = "Removed points", ls = 'none', marker = '.')
        plt.axhline(0, color = 'cyan', label = 'zero', ls = '--')
        #plt.plot(npch, result_ball(result.params)(npch), color = 'r', label = "Crystal ball fit")
        #plt.text(0, np.max(y) / 2, "Reduced χ² = %.3f" % result.redchi)
        plt.text(0, np.max(y) / 2 - np.max(y) / 10, "Minimum signal counts = " +str(min_counts))
        #plt.axvline(x = result.params['xbar'].value, label = "Estimate mean", color = 'red', ls = '--')
        #if result.params['xbar'].stderr is not None:
        #    plt.axvline(x = result.params['xbar'].value + result.params['xbar'].stderr, label = "Estimate error bounds", color = 'orange', ls = '--')
        #    plt.axvline(x = result.params['xbar'].value - result.params['xbar'].stderr, color = 'orange', ls = '--')
        
        if foil is not None:
            plt.title("Channel histogram for "+foil+" (NORMALIZED)")
        plt.xlabel("Channel #")
        plt.ylabel("Normalized counts")
        plt.legend()
        plt.show()


        ##########
        plt.errorbar(x, y, yerr = 1/weights, color = 'black', label = "Energy channel histogram", ls = 'none', marker = '.')
        plt.plot(npch, result_ball(result.params)(npch), color = 'r', label = "Crystal ball fit")
        plt.text(0, np.max(y) / 2, "Reduced χ² = %.3f" % result.redchi)
        plt.text(0, np.max(y) / 2 - np.max(y) / 10, "Minimum signal counts = " +str(min_counts))
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

    return Result(result.params['xbar'].value, stat = result.params['xbar'].stderr)#, result.redchi
"""

"""
#probably need to limit degree.
def poly_coeffs(energies, pp, cs, max_degree, dx, debug = False):
    # fucked due to floating point errors :(
    product_previous = pp[:]
    coefficient_sum = cs[:]
    coefficient = [0]
    test_p = 0.01
    threshold = 1e-8
    for E in energies:

        print ("-")
        #print (np.polyval(product_previous[::-1], 0))
        

        contribution = np.array(product_previous)/E**2
        coefficient = [0] + list(contribution)
        anticoefficient = [1] + list(-contribution)
        print ("Expected pp:", expected:= (np.polyval(product_previous[::-1], test_p) * (1 - np.polyval(coefficient[::-1], test_p))))
        #print (np.polyval(coefficient[::-1], 0))
        #print ("contribution:",coefficient)
        product_previous = np.convolve(product_previous, [1] + list(-contribution))
        product_previous = product_previous[:min(len(product_previous), max_degree)]
        print ("pp:", got:= (np.polyval(product_previous[::-1], test_p)))
        if abs(expected - got) > threshold:
            raise Exception
        print ("Normalization:")
        print (np.polyval(coefficient[::-1], test_p))
        print ("1-:", np.polyval([1] + list(-contribution), test_p))
        print (np.polyval(coefficient[::-1], test_p) + np.polyval(anticoefficient, test_p))
        total = np.array(coefficient[::-1]) + np.array(anticoefficient[::-1])
        print (np.polyval(total, test_p))
        #print ("product:",product_previous)
        coefficient_sum = list(np.array(coefficient_sum + [0] * (len(coefficient) - len(coefficient_sum))) + np.array(coefficient))
        coefficient_sum = coefficient_sum[:min(len(coefficient_sum), max_degree)]
        #print ("sum:", coefficient_sum)
        #print (np.polyval(coefficient_sum[::-1], 0))
        print ("Total:")
        print (np.polyval(coefficient_sum[::-1], test_p))
    return coefficient_sum#, product_previous

def test_poly_coeffs():
    energies = np.ones(20) #np.arange(5, 3, -0.1)
    pp = [1]
    cs = [0]
    print(poly_coeffs(energies, pp, cs, 50, 1))
"""


def solve_probability_distribution(Elists, p_scatters, max_degree, dx, plot = False):
    #make sure you order them by length

    p_coeff_list = []
    pp = [1]
    cs = [0]
    
    """    
    for i in range(len(Elists)):
        Es = Elists[i]
        Esprev = [] if i == 0 else Elists[i-1]
        new_Es = Es[len(Esprev):]
        print (len(Esprev))
        cs, pp = poly_coeffs(new_Es, pp, cs, max_degree)
        p_coeff_list.append(cs)
        print (p_coeff_list[-1][1])
        print ("-----")
        print (len(Es))"""
    for Elist in Elists:
        #print (pp)
        #print (cs)
        p_coeff_list.append(poly_coeffs(Elist, pp, cs, max_degree, dx, debug = len(p_coeff_list) == 0)[::-1])
    #going to try a pair just for now
    y1 = p_scatters[0]
    y2 = p_scatters[1]
    c1 = np.array(p_coeff_list[0])
    c2 = np.array(p_coeff_list[1])
    #c1 = c1[0:3]
    #c2 = c2[0:3]
    print (c1)
    print (c2)
    #print (c1 - c2)
    final_coeffs = y2*c1 - y1*c2
    #print (final_coeffs)
    if plot:
        plot_polys(c1, c2, y1, y2)

    roots = np.roots(final_coeffs)
    complex_threshold = 1e-8
    real_roots = np.real(roots)[np.abs(roots - np.real(roots)) < complex_threshold]
    positive_roots = real_roots[real_roots > 0]
    valid_roots = positive_roots[positive_roots < 1]
    print (valid_roots)


# Plots
def plot_polys(c1, c2, y1, y2):
    pvals = np.linspace(0, 1, 150)
    pevals1 = np.polyval(c1, pvals)
    pevals2 = np.polyval(c2, pvals)
    #print (max(pevals1))
    #print (max(pevals2))
    #pevals1 = np.polyval(c1, pvals)
    #pevals2 = np.polyval(c1 - c2, pvals)
    plt.plot(pvals, pevals1, label = "p1")
    plt.plot(pvals, pevals2, label = "p2")
    plt.legend()
    plt.show()
    plt.plot(pvals[1:], pevals2[1:]/pevals1[1:], label = "ratio")
    plt.axhline(y2 / y1, ls = '--', label = "target ratio")
    plt.legend()
    plt.show()