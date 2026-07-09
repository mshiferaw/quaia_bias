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
import argparse

## Define parameters
nthreads = 18
Delta = 200
delta_c = 1.686
NSIDE = 64
M_max = 16.8 #5
dpi = 0.01
cosmo = 'Planck'
Om0 = {'Planck': Planck18.Om0, 'Quaia': 0.343}
cosmo_planck18 = ccl.Cosmology(Omega_c=0.265, Omega_b=Planck18.Ob0, h=Planck18.h, Neff = Planck18.Neff, T_CMB = Planck18.Tcmb0.value, m_nu = Planck18.m_nu, sigma8=0.8, n_s=0.95)
hmf = ccl.halos.MassFuncTinker10()
fac_rand = 25
nlive = 3000
dlogz = 0.001
bins = np.logspace(9, M_max, 79)

# Define a custom argument type for a list of floats
def list_of_floats(arg):
    return np.array(list(map(float, arg.split(','))))
    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("corrfunc", help="two-point statistic", type = str)
    parser.add_argument("G", help="magnitude limit", type = float)
    parser.add_argument("rmin", help="minimum separation", type = float)
    parser.add_argument("rmax", help="maximum separation", type = float)
    parser.add_argument("nbins", help="number of separation bins", type = int)
    parser.add_argument("scale", help="linear or log scale", type = str)
    parser.add_argument("--L_method", help="max, minmax, split, or maxsplit", type = str)
    parser.add_argument("N_jk", help="number of jackkife regions", type = int)
    parser.add_argument("--mask", help="threshold for pixel-based selection function mask", type = float)
    parser.add_argument("--b", help="galactic-latitude based mask", type = int)
    parser.add_argument("cosmo", help="cosmology", type = str)
    parser.add_argument("--z_bins", help="redshift bins", type = float, nargs = '*')
    parser.add_argument("--pimax", help="pimax", type = int, default = 80)
    parser.add_argument("--n_zbins", help="number of z bins", type = int)
    parser.add_argument("--n_Lbins", help="number of L bins", type = int)
    args = parser.parse_args()

    print('\nG{}: computing {} in {} bins from {} to {} on a {} scale for {} jackknife regions...'.format(args.G, args.corrfunc, args.nbins, args.rmin, args.rmax, args.scale, args.N_jk))
    print('--L: {}\n'.format(args.L_method))

    if args.n_zbins is None:
        z_array = range(len(args.z_bins)-1)
        z_method = 'minmax'
    else:
        z_array = range(args.n_zbins)
        z_method = 'split'
    if args.n_Lbins is None:
        L_array = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])
        # L_bins = L_array[:0:-1]
        if args.L_method=='max':
            L_array = L_array[:0:-1]
        L_bins = L_array
    else:
        L_array = range(args.n_Lbins)
        L_bins = None
    pibins = np.arange(0, args.pimax+1, dpi)

    if args.scale == 'linear':
        rbins_chi2 = np.linspace(args.rmin, args.rmax, args.nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    else:
        rbins_chi2 = np.logspace(np.log10(args.rmin), np.log10(args.rmax), args.nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202

    # Use delete-one jackknife resampling
    ## Split the sky into RA stripes
    ra_bins = np.linspace(0, 360, args.N_jk+1)

    ## Compute $\xi(s)$ in each jackknife bin
    # cf_jackknife = {Lbin: [] for Lbin in L_bins}
    file = '_G{}_rmin{}_rmax{}_nbins{}'.format(args.G, args.rmin, args.rmax, args.nbins)
    
    if args.mask is not None:
        mask = 'mask{}'.format(args.mask/100)
        percentile = True
        mask_type = 'selfunc'
    else:
        mask = 'b{}'.format(args.b)
        percentile = False
        mask_type = 'b'

    print('z:', args.z_bins, z_method, args.n_zbins)
    print('L:', L_bins, args.L_method, args.n_Lbins)

    for j in z_array: 

        if args.L_method is None:

            tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j], True, ['z'], 
                                                                                                     bins = [args.z_bins], tab_gcat_type = 'data', 
                                                                                                      method = ['minmax'], fac_rand = fac_rand, 
                                                                                                  mask = args.mask, percentile = percentile, n_bins = [None], mask_type = mask_type, b = args.b)
            L = ''

            cf_jackknife=[]
        
            for i in range(args.N_jk): #len(ra_bins)-1):
            
                jackknife_data = (tab_datahi_mask_zbin0['ra'] >= ra_bins[i+1]) | (tab_datahi_mask_zbin0['ra'] < ra_bins[i])
                jackknife_rand = (tab_randhi_mask_zbin0['ra'] >= ra_bins[i+1]) | (tab_randhi_mask_zbin0['ra'] < ra_bins[i])
    
                if args.corrfunc == 'xi':
    
                    cf_jackknife.append(quaia.xi_s(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], key_zbin0, 
                                                     nthreads = nthreads, rbins = rbins_chi2, Om0 = Om0[args.cosmo], #Om0[cosmo], 
                                                              h = Planck18.h)) # h = h, Om0 = Omega_m_planck20))
                    
                    pi = '_zmin{}zmax{}{}N{}_{}_MC_{}'.format(args.z_bins[j], args.z_bins[j+1], L, args.N_jk, args.cosmo, mask)
                    
                elif args.corrfunc == 'wp':
                    
                    cf_jackknife.append(quaia.wp_rp(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], key_zbin0, nthreads = nthreads, 
                                                                              rbins = rbins_chi2, nbins = args.nbins, pimax = args.pimax, Om0 = Om0[args.cosmo], h = Planck18.h)[0])
                    pi = '_pimax{}_zmin{}zmax{}{}N{}_{}_MC_{}'.format(args.pimax, args.z_bins[j], args.z_bins[j+1], L, args.N_jk, args.cosmo, mask)
                    
                else:
                    cf_jackknife.append(quaia.w_theta(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], nthreads = nthreads, thetabins = rbins_chi2)) # h = h, Om0 = Omega_m_planck20))

                    pi = '_zmin{}zmax{}{}N{}_{}_MC_{}'.format(args.z_bins[j], args.z_bins[j+1], L, args.N_jk, args.cosmo, mask)

            cf_mean = np.nanmean(cf_jackknife, axis = 0)  # is it ok to ignore nan values?
            delta = cf_jackknife-cf_mean
            C = (args.N_jk-1)/args.N_jk*(delta.T @ delta)

            np.save('../results/cov_{}{}{}'.format(args.corrfunc, file, pi), C)
            
            print('saved to {}{}'.format(file, pi))
            
        else:

            if args.z_bins is None:
                z = '_zsplit{}bin{}'.format(args.n_zbins, j)
            else:
                z = '_zmin{}zmax{}'.format(args.z_bins[j], args.z_bins[j+1])

            # for k in range(len(L_array)-1):
            for k in range(len(L_array)):
                
                try: 
                        
                    # if args.L == 'minmax':
                    #     tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j, k], True, ['z', 'L'], 
                    #                                                                                      bins = [args.z_bins, L_array], tab_gcat_type = 'data', 
                    #                                                                                       method = ['minmax', args.L], fac_rand = fac_rand, 
                    #                                                                                       mask = args.mask, percentile = percentile, n_bins = [None, None], mask_type = mask_type, b = args.b)
                    #     L = '_Lmin{}Lmax{}_'.format(L_array[k], L_array[k+1])
                    # else:
        
                    #     tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j, k], True, ['z', 'L'], 
                    #                                                                                              bins = [args.z_bins, L_bins], tab_gcat_type = 'data', 
                    #                                                                                               method = ['minmax', args.L], fac_rand = fac_rand, 
                    #                                                                                               mask = args.mask, percentile = percentile, n_bins = [None, None], mask_type = mask_type, b = args.b)
                    #     L = '_Lmax{}_'.format(L_bins[k])

                    tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j, k], True, ['z', 'L'], bins = [args.z_bins, L_bins], tab_gcat_type = 'data', method = [z_method, args.L_method], fac_rand = fac_rand, mask = args.mask, percentile = percentile, n_bins = [args.n_zbins, args.n_Lbins], mask_type = mask_type, b = args.b)

                    if args.L_method == 'minmax':
                        L = '_Lmin{}Lmax{}_'.format(L_array[k], L_array[k+1])
                    elif args.L_method == 'max':
                        L = '_Lmax{}_'.format(L_bins[k])
                    elif args.L_method == 'split':
                        L = '_Lsplit{}bin{}_'.format(args.n_Lbins, k)
                    else:
                        L = '_Lmaxsplit{}bin{}_'.format(args.n_Lbins, k)

                    cf_jackknife=[]
        
                    for i in range(args.N_jk): #len(ra_bins)-1):
                    
                        jackknife_data = (tab_datahi_mask_zbin0['ra'] >= ra_bins[i+1]) | (tab_datahi_mask_zbin0['ra'] < ra_bins[i])
                        jackknife_rand = (tab_randhi_mask_zbin0['ra'] >= ra_bins[i+1]) | (tab_randhi_mask_zbin0['ra'] < ra_bins[i])
            
                        if args.corrfunc == 'xi':
            
                            cf_jackknife.append(quaia.xi_s(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], key_zbin0, #[L_bins[k]][j]
                                                             nthreads = nthreads, rbins = rbins_chi2, Om0 = Om0[args.cosmo], #Om0[cosmo], 
                                                                      h = Planck18.h)) # h = h, Om0 = Omega_m_planck20))
                            
                            # pi = '_zmin{}zmax{}{}N{}_{}_MC_{}'.format(args.z_bins[j], args.z_bins[j+1], L, args.N_jk, args.cosmo, mask)
                            pi = '{}{}N{}_{}_MC_{}'.format(z, L, args.N_jk, args.cosmo, mask)
                            
                        elif args.corrfunc == 'wp':
                            
                            cf_jackknife.append(quaia.wp_rp(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], key_zbin0, nthreads = nthreads, 
                                                                                      rbins = rbins_chi2, nbins = args.nbins, pimax = args.pimax, Om0 = Om0[args.cosmo], h = Planck18.h)[0])
                            # pi = '_pimax{}_zmin{}zmax{}{}N{}_{}_MC_{}'.format(args.pimax, args.z_bins[j], args.z_bins[j+1], L, args.N_jk, args.cosmo, mask)
                            pi = '_pimax{}{}{}N{}_{}_MC_{}'.format(args.pimax, z, L, args.N_jk, args.cosmo, mask)
                            
                        else:
                            cf_jackknife.append(quaia.w_theta(tab_datahi_mask_zbin0[jackknife_data], tab_randhi_mask_zbin0[jackknife_rand], nthreads = nthreads, thetabins = rbins_chi2)) # h = h, Om0 = Omega_m_planck20))

                            # pi = '_zmin{}zmax{}{}N{}_{}_MC_{}'.format(args.z_bins[j], args.z_bins[j+1], L, args.N_jk, args.cosmo, mask)
                            pi = '{}{}N{}_{}_MC_{}'.format(z, L, args.N_jk, args.cosmo, mask)
        
            # ## Compute the jackknife covariance matrix
                
                    cf_mean = np.nanmean(cf_jackknife, axis = 0)  # is it ok to ignore nan values?
                    delta = cf_jackknife-cf_mean
                    C = (args.N_jk-1)/args.N_jk*(delta.T @ delta)
        
            ## Save the full covariance matrix
                    np.save('../results/cov_{}{}{}'.format(args.corrfunc, file, pi), C)
                    
                    print('saved to {}{}\n'.format(file, pi))
                    
                except Exception as e:
    
                    # print('G{}_zmin{}zmax{}_Lmin{}Lmax{}:                   no selection function'.format(args.G, args.z_bins[j], args.z_bins[j+1], L_array[k], 
                    #                                                                                           L_array[k+1]))
                    print(e)
                    
    print('done!')
          
if __name__=='__main__':
    main()