import numpy as np
import math
import matplotlib.pyplot as plt

Angle =  [20, 25, 30, 40, 45, 50]
CountsPerSecond = [0.030165912518853696, 0.018148820326678763, 0.012191103789126852, 0.0034184802484364253, 0.002088424409435172, 0.0007995507867335604]
BackgroundError = [0.001005530417295123, 0.0009074410163339383, 0.0004942339373970346, 0.00012087101513057168, 9.532227594932321e-05, 0.00010375849904175974]
PoissonError = [0.0027992818674142913, 0.0029404492002642517, 0.0014456289034221312, 0.00014811724102963635, 0.00019454361827982264, 7.425156597072258e-05]
Angle2 = [20, 25, 30, 40, 45, 50]

EnergyError = [0.06373508366211407, 0.03834515468963849, 0.025757583810796045, 0.007222626681530248, 0.004412460732761214, 0.0016893053128335812]

CPS_Energy2 = [0.77959276, 0.46902904, 0.31506079, 0.08834549, 0.0539722,  0.02066319]
AngleErrors = [-0.12007401547169562, -0.03553648105480852, -0.013554149955260972, -0.003077374639828499, -0.0016952035340902715, -0.0009991759398191439]
Finalyerrors = [0.15616903948161573, 0.09517284737019262, 0.0490518982255622, 0.009276143232303677, 0.007327333254954428, 0.003837279001518818]



# Check 1:
print(np.array(CPS_Energy2)/np.array(CountsPerSecond))

energy = 25.8435
Eerror = math.sqrt(4.4640)



# Modify the earlier functions
BackgroundError = energy * np.array(BackgroundError)
PoissonError = energy * np.array(PoissonError)

#Fractional Uncertainties

b = abs(np.array(BackgroundError)/np.array(CPS_Energy2))
p = abs(np.array(PoissonError)/np.array(CPS_Energy2))
a = abs(np.array(AngleErrors)/np.array(CPS_Energy2))
e = abs(np.array(EnergyError)/np.array(CPS_Energy2))
f = abs(np.array(Finalyerrors)/np.array(CPS_Energy2))

errors = {
    'background' : b,
    'poisson' : p,
    'energy' : e,
    'angle' : a,
    'total' : f
}

for key in errors:
    plt.errorbar(Angle, errors[key], label=key + ' error')
plt.legend()
plt.ylabel('Fractional uncertainty')
plt.xlabel('Angle (degrees)')
plt.show()