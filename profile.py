from beam_profile_models import *
from numpy import sqrt, abs
import lmfit 
params = lmfit.Parameters() 
params.add( 'x0', value=0.5962950765376037) 
params.add( 'a_L', value=21.856218228626407) 
params.add( 'y0_L', value=153.30557430489304) 
params.add( 'a_R', value=-28.36507313403734) 
def profile (x): 
	return use_beam_model(x, 'linear', 'linear', params)