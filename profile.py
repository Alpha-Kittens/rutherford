from beam_profile_models import *
from numpy import sqrt, abs
import lmfit 
params = lmfit.Parameters() 
params.add( 'x0', value=0.2914487912491559) 
params.add( 'a_L', value=30.145550872464096) 
params.add( 'y0_L', value=161.8128845721189) 
params.add( 'a_R', value=-28.289148951801547) 
def profile (x): 
	return use_beam_model(x, 'linear', 'linear', params)