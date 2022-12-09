import matplotlib.pyplot as plt
import numpy as np


def plot_histogram(title, histogram, xaxis=None, vlines = [], show=True):
    print(title)
    plt.title("Energy channel histogram for " + title)
    plt.xlabel("MCA channel")
    plt.ylabel("Particle counts")
    plt.plot([], [], color = 'blue', label = "Poisson errors on counts")
    if xaxis is not None:
        width = xaxis[1] - xaxis[0]
        plt.bar(xaxis, histogram, width = width, color = 'cyan', ecolor = 'blue', label = "Histogram data", yerr = np.sqrt(histogram))
    else:
        plt.bar(range(len(histogram)), histogram, width = 1, color = 'cyan', ecolor = 'blue', label = "Histogram data", yerr = np.sqrt(histogram))
    for label, x, color in vlines:
        plt.axvline(x, color = color, label = label, ls = '--')
    plt.legend()
    if show:
        plt.show()


def rebin_hist(histogram, factor = 20):
    nbins = 2048/factor
    bin_centers = []
    counts = []
    for i in range(len(histogram)):
        bin_centers.append(i)
        counts.append(histogram[i])
    # To rebin to a smaller number of bins
    bin_centers_merged=np.asarray([max(a) for a in np.array_split(bin_centers,nbins)])
    counts_merged=np.asarray([sum(a) for a in np.array_split(counts,nbins)])

    return bin_centers_merged, counts_merged