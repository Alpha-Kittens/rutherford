import numpy as np
from scipy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt 

#def deconvolve (output, input, width = 5):
#    H = fft(smear(output, width))
##    G = fft(smear(input, width))
 #   F = H / G
#    return ifft(F)

def convolve(f1,f2,x,iMin=-10,iMax=10,iN=2000):
    step=(iMax-iMin)/iN
    pInt=0
    for i0 in range(iN):
            pX   = i0*step+iMin
            pVal = f1(x-pX)*f2(pX)
            pInt += pVal*step
    return pInt

def smear (arr, width):
    return [smear_i(arr, i, width) for i in range(len(arr))]

def smear_i(arr, i, width):
    min_i = max(i - width, 0)
    max_i = min(i + width, len(arr) - 1) + 1 # endpoints don't matter at all
    dists = np.array(range(min_i - i, max_i - i))
    return np.dot(arr[min_i : max_i], (1 - (dists / width)**2)) / np.sum(1 - (dists / width)**2)


if __name__ == '__main__':
    from data_loader import recursive_read
    r = {}
    recursive_read(r, "data", require = [0], reject = ["titanium", "2gold", "3gold", "iron"])
    output = np.array(r[("gold", 0, 1)][1]) / r[("gold", 0, 1)][0]
    input = np.array(r[("empty", 0, 1)][1]) / r[("empty", 0, 1)][0]
    x = np.array(range(len(output)))

    f = deconvolve(output, input, width = 10)

    print (f)
    plt.plot(x, smear(output, 10), label = "smeared output", color = 'r')
    plt.plot(x, smear(input, 10), label = "smeared input", color = 'b')
    plt.plot(x, f, label = "f", color = 'purple')
    plt.plot(x, np.convolve(input, f, mode = 'same'), label = "input*f", color = 'lime')
    plt.legend()
    plt.show()