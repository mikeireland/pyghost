"""Given (x,wave,matrices, slit_profile), extract the flux from each order. For 
readability, we keep this separate from the simulator.... but the"""

from __future__ import division, print_function
import ghostsim

class Extractor():
    """A class for each arm of the spectrograph. The initialisation function takes a 
    single string representing the configuration. For GHOST, it can be "red" or "blue"."""
    
    def __init__(self,arm):
        sim = ghostsim.Arm(arm)
        self.x_map,self.w_map,self.blaze,self.matrices = sim.spectral_format_with_matrix()
        
    def define_lenslets(self,positions,weights):
        print("Not implemented yet!")