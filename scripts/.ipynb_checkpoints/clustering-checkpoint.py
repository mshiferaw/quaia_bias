import numpy as np
import healpy as hp

import matplotlib
from matplotlib import pyplot as plt

from scipy import stats
import matplotlib.gridspec as gridspec
import astropy.cosmology.units as cu
import astropy.units as u
import site
site.addsitedir('/oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts') 
import quaia
import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
import pyccl as ccl
from matplotlib.lines import Line2D
from matplotlib import cm
from scipy.optimize import minimize
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
from scipy import linalg
from scipy.special import gamma
import argparse
from scipy import special
import time

# Model using step-function HOD
## Define parameters
nthreads = 18
Delta = 200
delta_c = 1.686
NSIDE = 64
M_max = 16.8 #5
dpi = 0.01
cosmo_planck18 = ccl.Cosmology(Omega_c=0.265, Omega_b=Planck18.Ob0, h=Planck18.h, Neff = Planck18.Neff, T_CMB = Planck18.Tcmb0.value, m_nu = Planck18.m_nu, sigma8=0.8, n_s=0.95)
Om0 = {'Planck': Planck18.Om0, 'Quaia': 0.343}
hmf = ccl.halos.MassFuncTinker10()
fac_rand = 25
nlive = 3000
dlogz = 0.001
bins = np.logspace(9, M_max, 79)

## Define functions
def a_median(tab_datahi_mask_zbin0, key_zbin0 = ''):
    return 1/(1+np.median(tab_datahi_mask_zbin0['redshift_quaia'+key_zbin0]))

def b_qso(xi_qq, xi_mm, f):
    return np.sqrt(xi_qq/xi_mm-4*f**2/45)-f/3

def nu(M, delta_c = delta_c, a = 1, cosmo = cosmo_planck18):
    sigma = cosmo.sigmaM(M, a = a)
    return delta_c/sigma

def n_g_theory(bins, n_total, a = 1, cosmo = cosmo_planck18, hmf = hmf, recenter = True):

    if recenter == True:
        bins = quaia.recenter(bins)
    nm = hmf(cosmo, bins, a)/cosmo['h']**3
    n_g = nm/(bins*np.log(10))*n_total
    n_g[np.isnan(n_g)]=0
    
    return n_g, nm

def bias(bins, n_m, a = 1, cosmo = cosmo_planck18, mod = False, recenter = True):

    n_g, nm = n_g_theory(bins, n_m, a = a, cosmo = cosmo, recenter = recenter)
    if recenter == True:
        print('recentering')
        bins = quaia.recenter(bins)
    n_g = np.trapz(n_g, bins)
    dn_dm = nm/(bins*np.log(10))

    b1_E_tinker = b1_E_func(nu(bins, a = a, cosmo = cosmo), mod = mod)

    n_m[~np.isfinite(n_m)]=0 
    integrand = dn_dm*n_m
    integrand[np.isnan(integrand)]=0
    
    return bias_function(n_g, integrand, b1_E_tinker, bins, recenter = recenter), n_g

def b1_E_func(nu, Delta = Delta, delta_c = delta_c, mod = False):

    delta_c = 1.686
    
    if mod == False:
        y = np.log10(Delta)
        A =  1.0+0.24*y*np.exp(-(4/y)**4)
        a = 0.44*y-0.88
        B = 0.183
        b =  1.5
        C = 0.019+0.107*y+0.19*np.exp(-(4/y)**4)
        c = 2.4
        
    else:
        A =  1.0
        a = 0.0906
        B = -4.5002
        b =  2.1419
        C = 4.9148
        c = 2.1419
        
    return 1-A*nu**a/(nu**a+delta_c**a)+B*nu**b+C*nu**c # Jose et al
    
def bias_function(n_g, integrand, b, bins, recenter = True):
    
    if recenter == True:
        print('recentering')
        bins = quaia.recenter(bins)
        
    return 1/n_g*np.trapz(integrand*b, bins)

def cf_2h(bins, n_m, a, cosmo, cf_matter, recenter = True):
    
    b_q, n_g = bias(bins, n_m, a, cosmo, recenter = recenter)
    return b_q**2*cf_matter, b_q

def M(bins, n_m, a = 1, cosmo = cosmo_planck18, mod = False, recenter = True, hmf = hmf):

    n_g, nm = n_g_theory(bins, n_m, a = a, cosmo = cosmo, recenter = recenter, hmf = hmf)
    if recenter == True:
        n_g = np.trapz(n_g, quaia.recenter(bins))
        dm = np.diff(bins) 
        dn_dm = nm/(quaia.recenter(bins)*np.logf(10))
    else:
        n_g = np.trapz(n_g, bins)
        dn_dm = nm/(bins*np.log(10))

    n_m[~np.isfinite(n_m)]=0 
    integrand = dn_dm*n_m
    integrand[np.isnan(integrand)]=0
    
    if recenter == True:
        return bias_function(n_g, integrand, quaia.recenter(bins), bins, recenter = recenter), n_g
    else:
        return bias_function(n_g, integrand, bins, bins, recenter = recenter), n_g

def HOD(x, ndim, bins_m = bins):

    if ndim == 2:
        fduty, M_min = x
        bins_m = np.logspace(M_min, M_max, 79)
        n_m=np.ones(np.shape(bins_m))*10**fduty
    else:
        fduty, M_min, sigma = x
        n_m = 10**fduty/2*special.erfc(np.log(10**M_min/bins_m)/(np.sqrt(2)*np.log(10)*sigma))
        
    return n_m, bins_m

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
    parser.add_argument("--L_method", help="threshold (max) or bin (minmax)", type = str)
    parser.add_argument("--mask", help="threshold for pixel-based selection function mask", type = float)
    parser.add_argument("MC", help="LMC/SMC mask", type = str)
    parser.add_argument("--b", help="galactic-latitude based mask", type = int)
    parser.add_argument("cosmo", help="cosmology", type = str)
    parser.add_argument("--z_bins", help="redshift bins", type = float, nargs = '*')
    parser.add_argument("--pimax", help="pimax", type = int, default = 80)
    parser.add_argument("--n_zbins", help="number of z bins", type = int)
    parser.add_argument("--n_Lbins", help="number of L bins", type = int)
    args = parser.parse_args()

    print('\nG{}: computing {} in {} bins from {} to {} on a {} scale...\n'.format(args.G, args.corrfunc, args.nbins, args.rmin, args.rmax, args.scale))
    
    # z_array = range(len(args.z_bins)-1)
    # L_array = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])
    # L_bins = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])[:0:-1]
    if args.n_zbins is None:
        z_array = range(len(args.z_bins)-1)
        z_method = 'minmax'
    else:
        z_array = range(args.n_zbins)
        z_method = 'split'
    if args.n_Lbins is None:
        L_array = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])
        if args.L_method=='max':
            L_array = L_array[:0:-1]
        L_bins = L_array
    else:
        L_array = range(args.n_Lbins)
        L_bins = None
    pibins = np.arange(0, args.pimax+1, dpi)
    
    # ## Calculate number density as $n_q \approx \int_0^{\infty} d M \frac{d n}{d M}\langle N(M)\rangle$
    ## Limit scales to $30\geq s \geq 80h^{-1}$ Mpc
    # Create the bins array
    if args.scale == 'linear':
        rbins_chi2 = np.linspace(args.rmin, args.rmax, args.nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    else:
        rbins_chi2 = np.logspace(np.log10(args.rmin), np.log10(args.rmax), args.nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202

    # cf_chi2, cf_matter_chi2 = {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # yerr_cf_chi2 = {L: [] for L in L_bins}
    # a_zbin = {L: [] for L in L_bins}
    # rpavg_chi2 = {L: [] for L in L_bins}
    H0 = Planck18.H0 # cosmo_planck18['H0'] * u.km/u.s/u.Mpc
    file = '_G{}_rmin{}_rmax{}_nbins{}'.format(args.G, args.rmin, args.rmax, args.nbins)
    
    if args.mask is not None:
        mask = 'mask{}'.format(args.mask/100)
        percentile = True
        mask_type = 'selfunc'
    else:
        mask = 'b{}'.format(args.b)
        percentile = False
        mask_type = 'b'

    for j in z_array:

        if args.L_method is None:

            tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j], True, ['z'], 
                                                                                                     bins = [args.z_bins], tab_gcat_type = 'data', 
                                                                                                      method = ['minmax'], fac_rand = fac_rand, 
                                                                                                  mask = args.mask, percentile = percentile, n_bins = [None], MC = args.MC, mask_type = mask_type, b = args.b)
            cf = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, nthreads = nthreads, thetabins = rbins_chi2)
            pi = '_zmin{}zmax{}_{}{}{}'.format(args.z_bins[j], args.z_bins[j+1], args.cosmo, args.MC, mask)

            np.save('../results/{}{}{}'.format(args.corrfunc, file, pi), cf)
                    
            print('saved to {}{}'.format(file, pi))

            a = a_median(tab_datahi_mask_zbin0)
            np.save('../results/a_G{}_zmin{}zmax{}_{}{}{}'.format(args.G, args.z_bins[j], args.z_bins[j+1], args.cosmo, args.MC, mask), a)

        else:

            if args.z_bins is None:
                z = '_zsplit{}bin{}'.format(args.n_zbins, j)
            else:
                z = '_zmin{}zmax{}'.format(args.z_bins[j], args.z_bins[j+1])


            # for i in range(len(L_array)-1):
            for i in range(len(L_array)):
    
                try:
                    
                    # if args.L == 'minmax':
        
                    #     tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j, i], True, ['z', 'L'], 
                    #                                                                                      bins = [args.z_bins, L_array], tab_gcat_type = 'data', 
                    #                                                                                       method = ['minmax', args.L], fac_rand = fac_rand, 
                    #                                                                                   mask = args.mask, percentile = percentile, n_bins = [None, None], MC = args.MC, mask_type = mask_type, b = args.b)
                    
                    #     L = '_Lmin{}Lmax{}_'.format(L_array[i], L_array[i+1])

                    #     a = a_median(tab_datahi_mask_zbin0)
                        
                    # else:
                    #     tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(args.G, [j, i], True, ['z', 'L'], 
                    #                                                                                      bins = [args.z_bins, L_bins], tab_gcat_type = 'data', 
                    #                                                                                       method = ['minmax', args.L], fac_rand = fac_rand, 
                    #                                                                                       mask = args.mask, percentile = percentile, n_bins = [None, None], MC = args.MC, mask_type = mask_type, b = args.b)
                    #     L = '_Lmax{}_'.format(L_bins[i])

        
                    tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(
                        args.G, [j, i], True, ['z', 'L'], bins = [args.z_bins, L_bins], tab_gcat_type = 'data', method = [z_method, args.L_method], fac_rand = fac_rand, mask = args.mask, percentile = percentile, n_bins = [args.n_zbins, args.n_Lbins], MC = args.MC, mask_type = mask_type, b = args.b)
                    
                    if args.L_method == 'minmax':
                        L = '_Lmin{}Lmax{}_'.format(L_array[i], L_array[i+1])
                    elif args.L_method == 'max':
                        L = '_Lmax{}_'.format(L_bins[i])
                    elif args.L_method == 'split':
                        L = '_Lsplit{}bin{}_'.format(args.n_Lbins, i)
                    else:
                        L = '_Lmaxsplit{}bin{}_'.format(args.n_Lbins, i)

                    a = a_median(tab_datahi_mask_zbin0)
                    # np.save('../results/a_G{}_zmin{}zmax{}{}{}{}{}'.format(args.G, args.z_bins[j], args.z_bins[j+1], L, args.cosmo, args.MC, mask), a)
                    np.save('../results/a_G{}{}{}{}{}{}'.format(args.G, z, L, args.cosmo, args.MC, mask), a)
                    
                    # quasar correlation function
                    if args.corrfunc == 'xi':
                        
                        cf, qq, qr = quaia.xi_s(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, 
                                   rbins = rbins_chi2, Om0 = Om0[args.cosmo], h = Planck18.h, error = True)
                        # pi = '_zmin{}zmax{}{}{}{}{}'.format(args.z_bins[j], args.z_bins[j+1], L, args.cosmo, args.MC, mask)
                        pi = '{}{}{}{}{}'.format(z, L, args.cosmo, args.MC, mask)
        
                    elif args.corrfunc == 'wp':
                        
                        # quasar correlation function
                        cf, rpavg_zbin_i = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, 
                                                                                      rbins = rbins_chi2, nbins = args.nbins, pimax = args.pimax, Om0 = Om0[args.cosmo], h = Planck18.h)
                        
                        # pi = '_pimax{}_zmin{}zmax{}{}{}{}{}'.format(args.pimax, args.z_bins[j], args.z_bins[j+1], L, args.cosmo, args.MC, mask)
                        pi = '_pimax{}{}{}{}{}{}'.format(args.pimax, z, L, args.cosmo, args.MC, mask)
        
                        np.save('../results/rpavg{}{}'.format(file, pi), rpavg_zbin_i)
        
                    else:
        
                        cf = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, nthreads = nthreads, thetabins = rbins_chi2)
                        # pi = '_zmin{}zmax{}{}{}{}{}'.format(args.z_bins[j], args.z_bins[j+1], L, args.cosmo, args.MC, mask)
                        pi = '{}{}{}{}{}'.format(z, L, args.cosmo, args.MC, mask)
        
                    np.save('../results/{}{}{}'.format(args.corrfunc, file, pi), cf)
                    
                    print('saved to {}{}'.format(file, pi))
    
                except Exception as e:
                    
                    print(e)

    print('done!')

if __name__=='__main__':
    main()
    