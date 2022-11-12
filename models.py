import numpy as np
import math

def const_cross_section(c, theta):
    # as yet unused
    return c / (np.sin(theta))**4
    
def cross_section(Z, E, theta, Z_alpha = 2, e = 1):
    # as yet unused
    return (Z * Z_alpha * e**2 / 4*E)**2 / (np.sin(theta))**4

def convolve(f, fdomain, g, gdomain):
    """
    Executes a convolution of two functions on given domains. 
    **Ensure the spacing of points in the domains are identical**
    The domain of valid convolution might be shorter than you expect. Consider extrapolating out the longer of the two arrays 
        to expand the domain of the convolution. 
    Phil sent us a notebook that may have a better version of this? Will check later. 
    Arguments:
        * `f` (function: float->float): First function for convolution
        * `fdomain` (np array of floats): Domain of first function. 
        * `g` (function: float->float): Second function for convolution
        * `gdomain` (np array of floats): Domain of second function.
    Returns:
        * `c` (np array of floats): Evaluation of convolution on convolution domain
        * `cdomain` (np array of floats): Domain of convolution function. 
    """
    # everything needs to be sorted
    # works better if the shorter of fdomain and gdomain is of odd length, I think?
    M, N = (gdomain, fdomain) if len(fdomain) > len(gdomain) else (fdomain, gdomain)
    a = math.floor((len(M)-1) / 2) #3: a=b=1. 4: a=1,b=2
    b = math.ceil((len(M)-1) / 2)
    return np.convolve(f(fdomain), g(gdomain))[len(M):len(N)+1], 1/2 * (N[a:-b] + N[b:-a])

def interpolate(p, pdomain):
    """
    Given the output of `convolve` above, returns a function that approximates the convolved function. 
    Valid only on convolution domain. 
    Presently unused.
    Arguments:
        `p` (function): same as output of `convolve`. i.e., evaluation of function 
        `pdomain` corresponding domain.
    Returns:
        `pfunc` (function: float -> float): Function which uses linear interpolation to approximate convolution.
        Compatible with np arrays. 
    """
    #p, pdomain = info
    #p, pdomain are np array
    def pfunc(x):
        try:
            eval = []
            for xi in x:
                if xi >= max(pdomain) or xi <= min(pdomain):
                    print(xi)
                    print ("can only interpolate between valid convolution bounds: must have x >", min(pdomain), "| x <", max(pdomain))
                    raise ValueError
                lows = pdomain <= xi
                highs = pdomain >= xi
                u = (xi - pdomain[lows][-1])/(pdomain[highs][0] - pdomain[lows][-1]) # between 0 and 1. 0 if at xlow, 1 if at xhigh
                eval.append(p[lows][-1] * (1 - u) + p[highs][0] * (u))
        except:
            if x >= max(pdomain) or x <= min(pdomain):
                print ("can only interpolate between valid convolution bounds: must have x >", min(pdomain), "| x <", max(pdomain))
                raise ValueError
            lows = pdomain <= x
            highs = pdomain >= x
            u = (x - pdomain[lows][-1])/(pdomain[highs][0] - pdomain[lows][-1]) # between 0 and 1. 0 if at xlow, 1 if at xhigh
            eval = p[lows][-1] * (1 - u) + p[highs][0] * (u)
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
    
