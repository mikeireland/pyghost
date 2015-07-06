"""Given (x,wave,matrices, slit_profile), extract the flux from each order. For 
readability, we keep this separate from the simulator.... but the simulator is
required in order to run this."""

from __future__ import division, print_function
import ghostsim
import numpy as np
import matplotlib.pyplot as plt
import pdb

class Extractor():
    """A class for each arm of the spectrograph. The initialisation function takes a 
    single string representing the configuration. For GHOST, it can be "red" or "blue".
    
    The extraction is defined by 3 key parameters: an "x_map", which is equivalent to
    2dFDR's tramlines and contains a physical x-coordinate for every y (dispersion direction)
    coordinate and order, and a "w_map", which is the wavelength corresponding to every y
    (dispersion direction) coordinate and order. """
    
    def __init__(self,arm, mode):
        self.sim = ghostsim.Arm(arm)
        self.x_map,self.w_map,self.blaze,self.matrices = self.sim.spectral_format_with_matrix()
        #Fill in the slit dimensions in "simulator pixel"s. based on if we are in the 
        #high or standard resolution mode.
        if mode == 'high':
            self.mode = mode
            self.lenslet_width = self.sim.lenslet_high_size
            self.nl = 28
            ## Set default profiles - object, sky and reference
            fluxes = np.zeros( (self.nl,3) )
            fluxes[2:2+19,0]=1.0
            fluxes[2+19:,1]=1.0
            fluxes[0,2]=1.0
            self.define_profile(fluxes)
            
        elif mode == 'std':
            self.mode = mode
            self.lenslet_width = self.sim.lenslet_std_size
            self.nl = 17
            ## Set default profiles - object 1, sky and object 2
            fluxes = np.zeros( (self.nl,3) )
            fluxes[0:7,0]  = 1.0
            fluxes[7:10,1] = 1.0
            fluxes[10:,2] = 1.0
            self.define_profile(fluxes)
            
        #Set some default pixel offsets for each lenslet, as used for a square lenslet profile
        ny = self.x_map.shape[1]
        nm = self.x_map.shape[0]
        pix_offset_ix = np.append(np.append([0],np.arange(1,self.nl).repeat(2)),self.nl)
        self.square_offsets = np.empty( (2*self.nl,nm) )
        # The [0,0] component of "matrices" measures the size of a detector pixel in the 
        # simulated slit image space. i.e. slitmicrons/detpix.
        for i in range(self.nl):
            self.square_offsets[:,i] = (pix_offset_ix - self.nl/2.0) * self.lenslet_width / self.matrices[i,self.x_map.shape[1]//2,0,0]
        self.sim_offsets = np.empty( (self.sim.im_slit_sz,nm) )
        im_slit_pix_in_microns = (np.arange(self.sim.im_slit_sz) - self.sim.im_slit_sz/2.0) * self.sim.microns_pix
        for i in range(nm):
            self.sim_offsets[:,i] = im_slit_pix_in_microns / self.matrices[i,self.x_map.shape[1]//2,0,0]
        #To aid in 2D extraction, let's explicitly compute the y offsets corresponding to these x offsets...
        #The "matrices" map pixels back to slit co-ordinates. 
        self.slit_tilt = np.zeros( (nm,ny) )
        for i in range(nm):
            for j in range(ny):
                invmat = np.linalg.inv( self.matrices[i,j] )
                #What happens to the +x direction?
                x_dir_map = np.dot(invmat,[1,0])
                self.slit_tilt[i,j] = x_dir_map[1]/x_dir_map[0]
        
    def define_profile(self,fluxes):
        """ Manually define the slit profile as used in lenslet extraction. As this is
        a low-level function, all lenslets must be defined. e.g. by convention, for the
        star lenslets of the high resolution mode, lenslets 0,1 and 21 through 27 would 
        be zero. Also """
        
        if fluxes.shape[0] != self.nl:
            print("Error: {0:s} resolution mode must have {1:d} lenslets".format(self.mode,self.nl))
        else:
            self.square_profile = np.empty( (fluxes.shape[0]*2, fluxes.shape[1]) )
            self.sim_profile = np.empty( (self.sim.im_slit_sz, fluxes.shape[1]) )
            for i in range(fluxes.shape[1]):
                self.square_profile[:,i] = np.array(fluxes[:,i]).repeat(2)
                im_slit=self.sim.make_lenslets(fluxes=fluxes[:,i], mode=self.mode)
                self.sim_profile[:,i] = np.sum(im_slit, axis=0)
        
    def one_d_extract(self, data, badpix=[], lenslet_profile='sim', rnoise=3.0):
        """ Extract flux by integrating down columns (the "y" direction), using an
        optimal extraction method.
        
        Parameters
        ----------
        data: numpy array
            Image data, transposed so that dispersion is in the "y" direction. 
        
        lenslet_profile: 'square' or 'sim'
            Shape of the profile of each fiber as used in the extraction. For a final
        implementation, 'measured' should be a possibility. 'square' assigns each
        pixel uniquely to a single lenslet. For testing only
        
        rnoise: float
            The assumed readout noise.
        
        WARNING: Binning not implemented yet"""
        ny = self.x_map.shape[1]
        nm = self.x_map.shape[0]
        nx = self.sim.szx
        
        #Number of "objects"
        no = self.square_profile.shape[1]
        extracted_flux = np.zeros( (nm,ny,no) )
        extracted_var = np.zeros( (nm,ny,no) )
        
        #Assuming that the data are in photo-electrons, construct a simple model for the
        #pixel inverse variance.
        pixel_inv_var = 1.0/(np.maximum(data,0) + rnoise**2)
        pixel_inv_var[badpix]=0.0
                
        #Loop through all orders then through all y pixels.
        for i in range(nm):
            print("Extracting order: {0:d}".format(i))
            #Based on the profile we're using, create the local offsets and profile vectors
            if lenslet_profile == 'square':
                offsets = self.square_offsets[:,i]
                profile = self.square_profile
            elif lenslet_profile == 'sim':
                offsets = self.sim_offsets[:,i]
                profile = self.sim_profile
            nx_cutout = 2*int( (np.max(offsets) - np.min(offsets))/2 ) + 2
            phi = np.empty( (nx_cutout,no) )
            for j in range(ny):
                #Check for NaNs
                if self.x_map[i,j] != self.x_map[i,j]:
                    extracted_var[i,j,:] = np.nan
                    continue
                #Create our column cutout for the data and the PSF
                x_ix = int(self.x_map[i,j]) - nx_cutout//2 + np.arange(nx_cutout,dtype=int) + nx//2
                for k in range(no):
                    phi[:,k] = np.interp(x_ix - self.x_map[i,j] - nx//2, offsets, profile[:,k])
                    phi[:,k] /= np.sum(phi[:,k])
                #Deal with edge effects...
                ww = np.where( (x_ix >= nx) | (x_ix < 0) )[0]
                x_ix[ww]=0
                phi[ww,:]=0.0
                #Cut out our data and inverse variance.
                col_data = data[j,x_ix]
                col_inv_var = pixel_inv_var[j,x_ix]
                #Fill in the "c" matrix and "b" vector from Sharp and Birchall equation 9
                #Simplify things by writing the sum in the computation of "b" as a matrix
                #multiplication. We can do this because we're content to invert the 
                #(small) matrix "c" here. Equation 17 from Sharp and Birchall 
                #doesn't make a lot of sense... so lets just calculate the variance in the
                #simple explicit way.
                col_inv_var_mat = np.reshape(col_inv_var.repeat(no), (nx_cutout,no) )
                b_mat = phi * col_inv_var_mat
                c_mat = np.dot(phi.T,phi*col_inv_var_mat)
                pixel_weights = np.dot(b_mat,np.linalg.inv(c_mat))
                extracted_flux[i,j,:] = np.dot(col_data,pixel_weights)
                extracted_var[i,j,:] = np.dot(1.0/np.maximum(col_inv_var,1e-12),pixel_weights**2)
        return extracted_flux, extracted_var
        
    def two_d_extract(self,  data, badpix=[], lenslet_profile='sim', rnoise=3.0):
        """ Extract using 2D information. The lenslet model used is a collapsed profile, 
        but where we take into account the slit shear/rotation.  """
        ny = self.x_map.shape[1]
        nm = self.x_map.shape[0]
        nx = self.sim.szx
        
        #Number of "objects"
        no = self.square_profile.shape[1]
        extracted_flux = np.zeros( (nm,ny,no) )
        extracted_var = np.zeros( (nm,ny,no) )
        
        #Assuming that the data are in photo-electrons, construct a simple model for the
        #pixel inverse variance.
        pixel_inv_var = 1.0/(np.maximum(data,0) + rnoise**2)
        pixel_inv_var[badpix]=0.0
                
        #Loop through all orders then through all y pixels.
        for i in range(nm):
            print("Extracting order: {0:d}".format(i))
            #Based on the profile we're using, create the local offsets and profile vectors
            if lenslet_profile == 'sim':
                offsets = self.sim_offsets[:,i]
                profile = self.sim_profile
            else:
                raise userwarning
            nx_cutout = 2*int( (np.max(offsets) - np.min(offsets))/2 ) + 2
            ny_cutout = 2*int(nx_cutout * np.nanmax(np.abs(self.slit_tilt)) / 2) + 3
            for j in range(ny):
                phi = np.zeros( (ny_cutout,nx_cutout,no) )
                #Check for NaNs
                if self.x_map[i,j] != self.x_map[i,j]:
                    extracted_var[i,j,:] = np.nan
                    continue
                #Create our column cutout for the data and the PSF
                x_ix = int(self.x_map[i,j]) - nx_cutout//2 + np.arange(nx_cutout,dtype=int) + nx//2
                y_ix = j + np.arange(ny_cutout, dtype=int) - ny_cutout//2
                for k in range(no):
                    x_prof = np.interp(x_ix - self.x_map[i,j] - nx//2, offsets, profile[:,k])
                    y_pix = (x_ix - self.x_map[i,j] - nx//2) * self.slit_tilt[i,j] + ny_cutout//2
                    frac_y_pix = y_pix - y_pix.astype(int)
                    subx_ix = np.arange(nx_cutout,dtype=int)
                    phi[y_pix.astype(int),subx_ix,k] = (1-frac_y_pix)*x_prof
                    phi[y_pix.astype(int)+1,subx_ix,k] = frac_y_pix*x_prof
                    phi[:,:,k] /= np.sum(phi[:,:,k])
                #Deal with edge effects...
                ww = np.where( (x_ix >= nx) | (x_ix < 0) )[0]
                x_ix[ww]=0
                phi[:,ww,:]=0.0
                ww = np.where( (y_ix >= ny) | (y_ix < 0) )[0]
                y_ix[ww]=0
                phi[ww,:,:]=0.0
                xy = np.meshgrid(y_ix, x_ix)
                #Cut out our data and inverse variance.
                col_data = data[xy].flatten()
                col_inv_var = pixel_inv_var[xy].flatten()
                #Fill in the "c" matrix and "b" vector from Sharp and Birchall equation 9
                #Simplify things by writing the sum in the computation of "b" as a matrix
                #multiplication. We can do this because we're content to invert the 
                #(small) matrix "c" here. Equation 17 from Sharp and Birchall 
                #doesn't make a lot of sense... so lets just calculate the variance in the
                #simple explicit way.
                col_inv_var_mat = np.reshape(col_inv_var.repeat(no), (ny_cutout*nx_cutout,no) )
                phi = phi.reshape( (ny_cutout*nx_cutout,no) )
                b_mat = phi * col_inv_var_mat
                c_mat = np.dot(phi.T,phi*col_inv_var_mat)
                pixel_weights = np.dot(b_mat,np.linalg.inv(c_mat))
                extracted_flux[i,j,:] = np.dot(col_data,pixel_weights)
                extracted_var[i,j,:] = np.dot(1.0/np.maximum(col_inv_var,1e-12),pixel_weights**2)
        return extracted_flux, extracted_var
        
        