import csv
from histogram_fitters import *
import lmfit

mode = 'e'
if mode == 'e':
    data_set_path = "digitized/e_howitzer.csv"
elif mode == 'ch':
    data_set_path = "digitized/howitzer.csv"


xs = []
ys = []
yup = 6.06767e-3
with open(data_set_path) as datafile:
    reader = csv.reader(datafile, delimiter = ',')
    for row in reader:
        #xs.append(round(float(row[0])))
        xs.append(float(row[0]))
        ys.append(float(row[1]))
rel = (yup - max(ys))/max(ys)
peak = 1/rel**2
scaling = peak/max(ys)

#histogram = [0] * 2048
#for i, x in enumerate(xs):
#    histogram[x] = scaling*ys[i]

#print (xs)
#print (ys)
scaled_crystal_ball = lambda x, α, n, xbar, σ, k : k*normalized_crystal_ball(x, α, n, xbar, σ)
model = lmfit.Model(scaled_crystal_ball)
params = crystal_ball_params(np.array(xs), np.array(ys))
params.add("k", value = scaling, min = 0)
#print (params)
result = model.fit([y*scaling for y in ys], x = np.array(xs), params = params, weights = [1/np.sqrt(y*scaling) for y in ys])
print(lmfit.fit_report(result))

import matplotlib.pyplot as plt


if mode == 'e':
    testx = np.linspace(4, 5.5, 200)
    plt.errorbar(xs, [y*scaling for y in ys], yerr = np.sqrt([y*scaling for y in ys]), color = 'black', ls = 'none', marker = '.')
    plt.plot(testx, result.params['k'].value * result_ball(result.params)(testx), color = 'r')
elif mode == 'ch':
    testx = np.array(range(2048))
    plt.errorbar(xs, ys, yerr = np.sqrt(ys)/np.sqrt(scaling), color = 'black', ls = 'none', marker = '.')
    plt.plot(testx, result_ball(result.params)(testx), color = 'r')
plt.show()
#print (max(histogram))
#optimal_energy_fit(histogram, plot = True)





