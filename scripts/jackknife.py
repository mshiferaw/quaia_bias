import numpy as np
import healpy as hp
from healpy.newvisufunc import projview
from astropy.table import Table

import matplotlib

from Corrfunc import mocks
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from scipy import constants
from scipy import stats
import scipy.interpolate as interp
from random import sample
import matplotlib.gridspec as gridspec
from Corrfunc import theory
from astropy.cosmology import FlatLambdaCDM
from Corrfunc.utils import convert_3d_counts_to_cf
import astropy.cosmology.units as cu
import astropy.units as u
import site
site.addsitedir('/oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts') 
import quaia
from astropy.coordinates import SkyCoord
import seaborn as sns
from astropy.cosmology import Planck18
import pyccl as ccl

def main():

    # Create the bins array
    rbins = np.linspace(rmin, rmax, nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202

    cosmo_planck18 = ccl.Cosmology(Omega_c=0.265, Omega_b=Planck18.Ob0, h=Planck18.h, Neff = Planck18.Neff, T_CMB = Planck18.Tcmb0.value, m_nu = Planck18.m_nu, sigma8=0.8, n_s=0.95)
    cosmo = 'Planck' #'Planck'

    # Use delete-one jackknife resampling
    ## Split the sky into RA stripes
    ra_bins = np.linspace(0, 360, N+1)

    ## Compute $\xi(s)$ in each jackknife bin
    cf_jackknife = {zbin: [] for zbin in zbin_array[:-1]}
    for j in range(len(zbin_array)-1):
        
        tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, _, _= quaia.make_bins(G_hi, j, True, 'z', bins = zbin_array, 
                                                                                           tab_gcat_type = 'data', method = 'minmax', 
                                                                                        fac_rand = fac_rand, mask = mask*100, percentile = True)
    
        for i in range(N): #len(ra_bins)-1):
        
            jackknife_data = (tab_datahi_mask_zbin0['ra'] >= ra_bins[i+1]) | (tab_datahi_mask_zbin0['ra'] < ra_bins[i])
            jackknife_rand = (tab_randhi_mask_zbin0['ra'] >= ra_bins[i+1]) | (tab_randhi_mask_zbin0['ra'] < ra_bins[i])
    
            cf_jackknife[zbin_array[j]].append(quaia.xi_s(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], key_zbin0, 
                                                 nthreads = nthreads, rbins = rbins, Om0 = Planck18.Om0, #Om0[cosmo], 
                                                          h = Planck18.h)) # h = h, Om0 = Omega_m_planck20))

    ## Compute the jackknife covariance matrix
    C = {zbin: [] for zbin in zbin_array[:-1]}
    for zbin in zbin_array[:-1]:
        
        cf_mean = np.nanmean(cf_jackknife[zbin], axis = 0)  # is it ok to ignore nan values?
        delta = cf_jackknife[zbin]-cf_mean
        C[zbin] = (N-1)/N*(delta.T @ delta)

    ## Save the full covariance matrix
    for i in range(len(zbin_array)-1):
    
        np.save('../results/cov_xi_G{}_rmin{}_rmax{}_nbins{}_zmin{}zmax{}_N{}_{}'.format(G_hi, rmin, rmax, nbins, zbin_array[i], zbin_array[i+1], N, cosmo), C[zbin_array[i]])

if __name__=='__main__':
    main()