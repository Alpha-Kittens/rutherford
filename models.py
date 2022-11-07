import numpy as np
import math

def const_cross_section(c, theta):
    return c / (np.sin(theta))**4
    
def cross_section(Z, E, theta, Z_alpha = 2, e = 1):
    return (Z * Z_alpha * e**2 / 4*E)**2 / (np.sin(theta))**4

def convolve(f, fdomain, g, gdomain):
    # everything needs to be sorted
    # works better if the shorter of fdomain and gdomain is of odd length, I think?
    M, N = (gdomain, fdomain) if len(fdomain) > len(gdomain) else (fdomain, gdomain)
    a = math.floor((len(M)-1) / 2) #3: a=b=1. 4: a=1,b=2
    b = math.ceil((len(M)-1) / 2)
    return np.convolve(f(fdomain), g(gdomain))[len(M):len(N)+1], 1/2 * (N[a:-b] + N[b:-a])

def interpolate(info):
    p, pdomain = info
    #p, pdomain are np array
    def pfunc(x):
        eval = []
        for xi in x:
            if xi >= max(pdomain) or xi <= min(pdomain):
                print ("can only interpolate between valid convolution bounds:  must have x >", min(pdomain), "| x <", max(pdomain))
                raise ValueError
            lows = pdomain <= x
            highs = pdomain >= x
            u = (x - pdomain[lows][-1])/(pdomain[highs][0] - pdomain[lows][-1]) # between 0 and 1. 0 if at xlow, 1 if at xhigh
            eval.append(p[lows][-1] * (1 - u) + p[highs][0] * (u))
        return eval
    return pfunc, pdomain
        



if __name__ == '__main__':
    f = lambda x : x + 1
    g = lambda x : 2*x
    fdomain = np.array([-1, 0, 1, 2])
    gdomain = np.array([2, 3, 4, 5, 6])
    print (convolve(f, fdomain, g, gdomain))


    #import matplotlib.pyplot as plt
    #def bp(x):
    #    return np.where(abs(x) > 6, 0, 6 - x**2 / np.abs(x))
    #cs = lambda x : 1 / (np.sin(x/2 * np.pi / 180))**4
    #angles = np.linspace(-10, 70, 160)
    #bp_domain = angles[angles < 10]
    #cs_domain = angles[angles > 25]
    #prediction = convolve(bp, bp_domain, cs, cs_domain)
    #M = len(bp_domain)
    #N = len(cs_domain)
    #prediction = np.convolve(bp(cs_domain), cs(bp_domain), mode = 'full')
    #print (len(prediction))
    #print (len(angles))
    #print (len(cs_domain))
    #print (len(bp_domain))
    #valid = prediction[M:N+1]
    #a = np.linspace(30, 90, 1000)
    #print (l := len(valid))
    #print (l + 2*(M-1))
    #plt.plot(a, cs(a))
    #plt.plot(bp_domain, bp(bp_domain), label = "beam profile", color = 'r')
    #plt.plot(cs_domain[cs_domain > 10], cs(cs_domain[cs_domain > 10]), label = "cross section", color = 'b')
    #plt.plot(cs_domain[cs_domain > 30], prediction[cs_domain > 30], label = "convolution", color = 'purple')
    #print (prediction[np.logical_and(cs_domain > 30, prediction > 0)])
    #print (cs(cs_domain))
    #print (bp(bp_domain))
    #print (prediction)
    #plt.legend()
    #plt.show()
    
