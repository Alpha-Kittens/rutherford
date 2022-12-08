import lmfit
import math
import numpy as np
import matplotlib.pyplot as plt
#Old data
#x = [math.log(i/10) for i in range(4, 10)]
#y = [-0.1-10/16/10, -0.1+(4/16+1/32)/10, 0-2/16/10, 0+(7/16+1/32)/10, 0.1-3/16/10, 0.1+3/16/10]
#yerr = [1/32/math.sqrt(12)] * len(x)
x = [0.5525328183974436, 0.6824347533073096, 0.7855910129686542, 0.9605782549513739, 0.421166818858211, 0.6048663001275308, 1.174543202064665, 0.46105789009622866, 0.4996776867506683, 0.7546216443429123, 0.886335722493827]
y = [-0.03750000000000003, 0.028571428571428525, 0.08124999999999999, 0.1419642857142857, -0.14464285714285713, -0.006250000000000033, 0.1910714285714286, -0.10714285714285715, -0.07500000000000001, 0.06607142857142856, 0.11785714285714288]
yerr = [1/math.sqrt(12) * 0.0008928571428571397  for yi in y]
# upper bound:
xerr = [1/math.sqrt(12) * 0.016876130131426192 for xi in x]


print (x)
print (y)
print (yerr)
def linear(x, a, b): 
    print (x, a, b)
    return a*np.log(x)+b
model = lmfit.Model(linear)
result = model.fit(y, x=x, weights = [1/err for err in yerr], a = 1, b = -1)
print (lmfit.fit_report(result))
slope = result.params['a'].value
new_yerr = np.sqrt(np.array(yerr)**2 + slope * (np.array(xerr))**2)
new_result = model.fit(y, x=x, weights = [1/err for err in new_yerr], a = 1, b = -1)
print(lmfit.fit_report(new_result))
plt.errorbar(x, y, xerr = xerr, yerr = yerr, ls = 'none')
plt.show()
plt.errorbar(np.log(x), y, yerr = new_yerr, ls = 'none')
plt.show()