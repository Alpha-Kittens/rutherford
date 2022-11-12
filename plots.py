import matplotlib.pyplot as plt
import numpy as np

def plot_histogram(metadata, histogram, vlines = []):
    plt.title("Energy channel histogram for "+metadata[0]+" target, "+str(metadata[1])+"° angle")
    plt.xlabel("MCA channel")
    plt.ylabel("Particle counts")
    plt.plot([], [], color = 'blue', label = "Poisson errors on counts")
    plt.bar(range(len(histogram)), histogram, width = 1, color = 'cyan', ecolor = 'blue', label = "Histogram data", yerr = np.sqrt(histogram))
    for label, x, color in vlines:
        plt.axvline(x, color = color, label = label, ls = '--')
    plt.legend()
    plt.show()
