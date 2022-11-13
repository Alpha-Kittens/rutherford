import lmfit
import math
x = [math.log(i/10) for i in range(4, 10)]
y = [-0.1-10/16/10, -0.1+(4/16+1/32)/10, 0-2/16/10, 0+(7/16+1/32)/10, 0.1-3/16/10, 0.1+3/16/10]
yerr = [1/32/math.sqrt(12)] * len(x)

print (x)
print (y)
print (yerr)
def linear(x, a, b): 
    print (x, a, b)
    return a*x+b
model = lmfit.Model(linear)
result = model.fit(y, x=x, weights = [1/err for err in yerr], a = 1, b = -1)
print (lmfit.fit_report(result))