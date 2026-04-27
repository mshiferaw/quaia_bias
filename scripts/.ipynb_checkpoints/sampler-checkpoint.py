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
from dynesty import NestedSampler
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
from astropy.coordinates import SkyCoord
from scipy import linalg
from scipy.special import gamma
import argparse
from scipy import special
import time

# Model using step-function HOD
## Define parameters
G_hi = 20.5
nthreads = 18
mask = 0.5
Delta = 200
delta_c = 1.686
NSIDE = 64
M_max = 15
pimax = 80 # 40 # 5000.0 # 40.0\
dpi = 0.01
cosmo = 'Planck'
N_jk = 50 #40 #20 #65 #50 #30 #15 #80
cosmo_planck18 = ccl.Cosmology(Omega_c=0.265, Omega_b=Planck18.Ob0, h=Planck18.h, Neff = Planck18.Neff, T_CMB = Planck18.Tcmb0.value, m_nu = Planck18.m_nu, sigma8=0.8, n_s=0.95)
hmf = ccl.halos.MassFuncTinker10()
rmin = 30 #-delta_min*offset # 0.1 # start higher # 1.4 # 30
rmax = 80 #+delta_max*offset # 200
nbins = 10 #20 + int(delta_min+delta_max) #90 #20
hartlap = (N_jk-nbins-2)/(N_jk-1)
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

    # n_m=np.zeros(np.shape(quaia.recenter(bins)))
    # n_m[quaia.recenter(bins)>=M_min]=fduty
    if ndim == 2:
        # bins_m, fduty = x
        fduty, M_min = x
        bins_m = np.logspace(M_min, M_max, 79)
        n_m=np.ones(np.shape(bins_m))*10**fduty
    else:
        # fduty, M_min, M, sigma = x
        fduty, M_min, sigma = x
        n_m = 10**fduty/2*special.erfc(np.log(10**M_min/bins_m)/(np.sqrt(2)*np.log(10)*sigma))
        
    return n_m, bins_m
    
# def HOD_step(bins_m, fduty): #M_max = M_max):

#     # bins_m = np.logspace(M_min, M_max, 79)
#     n_m=np.ones(np.shape(bins_m))*10**fduty
#     return n_m

# def HOD_erfc(fduty, M_min, M, sigma):
    
#     return 10**fduty/2*special.erfc(np.log(10**M_min/M)/(np.sqrt(2)*np.log(10)*sigma))

# def n_DM(cosmo, a, M_min, hmf = hmf, M_max = M_max, N = 79): #16.8
    
#     bins = np.logspace(M_min, M_max, N)
#     nm = hmf(cosmo, bins, a)/cosmo['h']**3 #/u.Mpc**3 # convert to h/Mpc**3
#     dn_dm = nm/(bins*np.log(10))
    
#     return np.trapz(dn_dm, bins)

def loglike(x, a, cf_matter, log_det_C, cf, n_q, nbins, solve, penalty, ndim, cosmo = cosmo_planck18, error = 0.05, M_max = M_max, hartlap = hartlap, N = N_jk):
    
    # if HOD == '_step':
    # if ndim == 2:

    #     fduty, M_min = x
    #     bins_m = np.logspace(M_min, M_max, 79)
    #     # n_m=np.ones(np.shape(bins_m))*10**fduty
    #     # n_m = HOD(bins_m, fduty)
    #     params = bins_m, fduty
        
    # else:

    #     fduty, M_min, sigma = x
    #     params = 
        
    n_m, bins_m = HOD(x, ndim)
    n_q_theory = M(bins_m, n_m, a = a, cosmo = cosmo, recenter = False)[1]
    if penalty == '_soft':
        lpenalty = 0.5*((n_q_theory-n_q)/(error*n_q))**2
    else:
        lpenalty = np.abs(n_q_theory-n_q)/n_q
        if lpenalty > error:
            return -np.inf
        
    cf_2h_zbin, _ = cf_2h(bins_m, n_m, a, cosmo, cf_matter, recenter = False)
    delta = cf-cf_2h_zbin
                
    c_p = gamma(N/2)/((np.pi*(N-1))**(nbins/2)*gamma((N-nbins)/2))
    lnorm = np.log(c_p)-0.5*log_det_C  
    
    return -N*0.5 * np.log(1+(delta.T @ solve(delta))/(N-1)) + lnorm - lpenalty

# def ptform_step(u):
#     """Transforms the uniform random variables `u ~ Unif[0., 1.)`
#     to the parameters of interest."""

#     x = np.array(u)  # copy u 

#     # Uniform from 0 to 1
#     x[0] =  6* u[0] - 6 # scale and shift to [-6., 0.)

#     # Uniform from 9 to 15
#     x[1] = 6 * u[1] + 9  # scale and shift to [9., 15.)

#     return x

def ptform(u, ndim):
    """Transforms the uniform random variables `u ~ Unif[0., 1.)`
    to the parameters of interest."""

    x = np.array(u)  # copy u 

    # Uniform from -6 to 0
    x[0] = 6* u[0] - 6 # scale and shift to [-6., 0.)

    # Uniform from 9 to 15
    x[1] = 6 * u[1] + 9  # scale and shift to [9., 15.)

    if ndim == 3:
        
        # Uniform from 0 to 1.4
        x[2] = 1.4 * u[2]  # scale to [0., 1.4.)

    return x
    
def f0_GLS_solve(x, a, cf_matter, C, cf, ndim, cosmo = cosmo_planck18, M_max = M_max, hartlap = hartlap, bins_m = bins):  #quaia):
        
    # if HOD == '_step':

    #     fduty, M_min = x
    #     bins_m = np.logspace(M_min, M_max, 79)
    #     # n_m=np.ones(np.shape(bins_m))*10**fduty
    #     n_m = HOD(bins_m, fduty)
        
    # else:

    #     fduty, M_min, sigma = x
    #     n_m = HOD_scatter(fduty, M_min, M, sigma)

    n_m, bins_m = HOD(x, ndim)
    cf_2h_zbin, _ = cf_2h(bins_m, n_m, a, cosmo, cf_matter, recenter = False)
    delta = cf-cf_2h_zbin
        
    return delta.T @ (hartlap*np.linalg.solve(C, delta))

def constraint(x, a, n_q, ndim, cosmo = cosmo_planck18, error = 0.05, M_max = M_max, bins_m = bins):
    
    # if HOD == '_step':

    #     fduty, M_min = x
    #     bins_m = np.logspace(M_min, M_max, 79)
    #     # n_m=np.ones(np.shape(bins_m))*10**fduty
    #     n_m = HOD(bins_m, fduty)
        
    # else:

    #     fduty, M_min, sigma = x
    #     n_m = HOD_scatter(fduty, M_min, M, sigma)

    n_m, bins_m = HOD(x, ndim)
    n_q_theory = M(bins_m, n_m, a = a, cosmo = cosmo, recenter = False)[1]

    return error-np.abs(n_q_theory-n_q)/n_q
    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("corrfunc", help="two-point statistic", type = str)
    parser.add_argument("penalty", help="penalty for number density constraint", type = str)
    # parser.add_argument("HOD", help="HOD type", type = str)
    parser.add_argument("ndim", help="number of HOD parameters", type = int)
    args = parser.parse_args()

    z_bins = np.array([0.0, 1.0, 2.0, 3.0, 4.6])
    z_array = range(len(z_bins)-1)
    L_bins = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])[:0:-1]
    pibins = np.arange(0, pimax+1, dpi)
    
    ## Calculate number density as $n_q \approx \int_0^{\infty} d M \frac{d n}{d M}\langle N(M)\rangle$
    tab_datahi_mask_zbins, selfunc_hi_bins = {L: [] for L in L_bins}, {L: [] for L in L_bins}

    for i in range(len(L_bins)):
            
        for j in z_array:
            
            tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [j, i], True, ['z', 'L'], 
                                                                                             bins = [z_bins, L_bins], tab_gcat_type = 'data', 
                                                                                              method = ['minmax', 'max'], fac_rand = fac_rand, 
                                                                                              mask = mask*100, percentile = True, n_bins = [None, None])
            tab_datahi_mask_zbins[L_bins[i]].append(tab_datahi_mask_zbin0['redshift_quaia'])
            selfunc_hi_bins[L_bins[i]].append(selfunc_hi_bin0)

    LMC = SkyCoord(['05 23 34.6 -69 45 22'], unit=(u.hourangle, u.deg)) #https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Large+Magellanic+Cloud
    LMC_radius = 9*u.deg
    SMC = SkyCoord(['00 52 38.0 -72 48 01'], unit=(u.hourangle, u.deg)) # https://simbad.u-strasbg.fr/simbad/sim-id?Ident=small+Magellanic+Cloud
    SMC_radius = 5*u.deg
    NPIX = hp.nside2npix(NSIDE)
    theta, phi = hp.pix2ang(NSIDE, range(NPIX), lonlat=True)
    
    c = SkyCoord(ra=theta*u.degree, dec=phi*u.degree)
    LMC_idx, _, _, _ = SkyCoord.search_around_sky(LMC, c, LMC_radius)
    SMC_idx, _, _, _ = SkyCoord.search_around_sky(SMC, c, SMC_radius)
    MC_idx = np.append(LMC_idx,SMC_idx)
    MC_mask = np.full_like(range(NPIX), True, dtype = bool)

    n_q, n_q_err, V_eff, N_q = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}

    for i in range(len(L_bins)):
            
        for j in z_array:   
    
            N = len(tab_datahi_mask_zbins[L_bins[i]][j])
            cp = np.percentile(selfunc_hi_bins[L_bins[i]][j][MC_mask], mask*100)
            selfunc_mask = selfunc_hi_bins[L_bins[i]][j]>=cp
            f_sky = 1/(4*np.pi)*np.sum(selfunc_hi_bins[L_bins[i]][j][selfunc_mask & MC_mask]*hp.nside2pixarea(NSIDE))
            print(f_sky) # f_sky)
            V = f_sky*(Planck18.comoving_volume(z_bins[j+1])-Planck18.comoving_volume(z_bins[j]))*Planck18.h**3/u.h**3
            N_q[L_bins[i]].append(N)
            n_q[L_bins[i]].append(N/V.value)
            n_q_err[L_bins[i]].append(np.sqrt(N)/V.value) # Poisson error
            V_eff[L_bins[i]].append(V)
            
    MC_mask[MC_idx] = False

    ## Limit scales to $30\geq s \geq 80h^{-1}$ Mpc
    # Create the bins array
    rbins_chi2 = np.linspace(rmin, rmax, nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202

    cf_chi2, cf_matter_chi2 = {L: [] for L in L_bins}, {L: [] for L in L_bins}
    yerr_cf_chi2 = {L: [] for L in L_bins}
    # , b_cf_chi2 = {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # b_cf_q_chi2 = {L: [] for L in L_bins}
    a_zbin = {L: [] for L in L_bins}
    rpavg_chi2 = {L: [] for L in L_bins}
    H0 = Planck18.H0 # cosmo_planck18['H0'] * u.km/u.s/u.Mpc
        
    for i in range(len(L_bins)):
                
        for j in z_array:
            
            tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [j, i], True, ['z', 'L'], 
                                                                                             bins = [z_bins, L_bins], tab_gcat_type = 'data', 
                                                                                              method = ['minmax', 'max'], fac_rand = fac_rand, 
                                                                                              mask = mask*100, percentile = True, n_bins = [None, None])
    
            # scale factor
            a = a_median(tab_datahi_mask_zbin0)
            a_zbin[L_bins[i]].append(a)
            
            # quasar correlation function
            if args.corrfunc == 'xi':
                
                cf, qq, qr = quaia.xi_s(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, 
                           rbins = rbins_chi2, Om0 = Planck18.Om0, h = Planck18.h, error = True)
                cf_chi2[L_bins[i]].append(cf)
    
                # matter correlation function
                rbins_Mpc = (quaia.recenter(rbins_chi2)*u.Mpc/cu.littleh).to(u.Mpc, cu.with_H0(H0))
                cf_matter = ccl.correlations.correlation_3d(cosmo_planck18, r=rbins_Mpc, a=a)
                pi = ''
                
            else:
                # quasar correlation function
                wp_zbin_i, rpavg_zbin_i = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, 
                                                                              rbins = rbins_chi2, nbins = nbins, pimax = pimax, Om0 = Planck18.Om0, h = Planck18.h)
                cf_chi2[L_bins[i]].append(wp_zbin_i)
                rpavg_chi2[L_bins[i]].append(rpavg_zbin_i)
            
                # matter correlation function   
                s = np.sqrt(np.sum(np.meshgrid(np.array(rpavg_zbin_i)**2, quaia.recenter(pibins)**2), axis = 0))
                s_Mpc = (s*u.Mpc/cu.littleh).to(u.Mpc, cu.with_H0(H0))
                
                xirppi_zbin_i = ccl.correlations.correlation_3d(cosmo_planck18, r=np.ravel(s_Mpc), a=a_median(tab_datahi_mask_zbin0))
                cf_matter = 2*np.trapz(np.reshape(xirppi_zbin_i, np.shape(s_Mpc)), quaia.recenter(pibins), axis = 0)
                pi = 'pimax{}_'.format(pimax)
        
            cf_matter_chi2[L_bins[i]].append(cf_matter)
        
            # # effective bias 
            # b = np.sqrt(cf/cf_matter)
            # b_cf_chi2[L_bins[i]].append(b)
        
            # # qso bias
            # b_cf_q_chi2[L_bins[i]].append(b_qso(cf, cf_matter, ccl.growth_rate(cosmo_planck18, a)))

    # Perform nested sampling with dynesty
    # if args.ndim == 2:
    #     # ndim = 2
    #     ptform = lambda u: ptform_step(u)
    # else:
    #     # ndim = 3
    #     ptform = lambda u: ptform_erfc(u)

    ## Across all $z$
    f_duty, M_min, sigma = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}
    labels = ['log10(f_duty)', 'log10(M_min)', 'sigma']
    err = {i: {L: [] for L in L_bins} for i in labels}
    f_duty_samples, M_min_samples, sigma_samples = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}
    maxL_params = {L: [] for L in L_bins}
    file = '_G{}_rmin{}_rmax{}_nbins{}_{}'.format(G_hi, rmin, rmax, nbins, pi)
        
    for i in range(len(L_bins)):

        print('M_i<={}'.format(L_bins[i]))
        
        for j in z_array:

            C = np.load('../results/cov_{}{}zmin{}zmax{}_Lmax{}_N{}_{}_MC_mask.npy'.format(args.corrfunc, file, z_bins[j], z_bins[j+1], L_bins[i], N_jk, cosmo))

            # check if cov is hermitian and positive definite
            if linalg.ishermitian(C) and np.all(np.linalg.eigvals(C) > 0):
                C_factor = linalg.cho_factor(C)
                log_det_C = 2 * np.sum(np.log(np.diag(C_factor[0])))
                solve = lambda rhs: linalg.cho_solve(C_factor, rhs)
            else:
                print('using LU decomposition')
                C_factor = linalg.lu_factor(C)
                log_det_C = np.sum(np.log(np.abs(np.diag(C_factor[0])))) # does abs value make sense
                solve = lambda rhs: linalg.lu_solve(C_factor, rhs)
        
            # initialize our nested sampler
            # if args.ndim == 2:
            sampler = NestedSampler(loglike, ptform, args.ndim, logl_args=[a_zbin[L_bins[i]][j], cf_matter_chi2[L_bins[i]][j], log_det_C, cf_chi2[L_bins[i]][j], n_q[L_bins[i]][j], nbins, solve, args.penalty, args.ndim], nlive = nlive, ptform_args = [args.ndim])
            # else:
            #     sampler = NestedSampler(loglike, ptform, ndim, logl_args=[a_zbin[L_bins[i]][j], cf_matter_chi2[L_bins[i]][j], log_det_C, cf_chi2[L_bins[i]][j], n_q[L_bins[i]][j], nbins, solve, args.penalty, args.HOD], nlive = nlive)
            
            # run the sampler with checkpointing
            start = time.time()
            sampler.run_nested(dlogz = dlogz, print_progress = False)
            results = sampler.results
            print(f"Iterations: {results.niter}")
            print(f"Time: {time.time() - start:.2f}s")
            
            # Print out a summary of the results.
            results.summary()        
            samples, weights = results.samples, results.importance_weights()
        
            # Get median and 1-sigma (16th, 50th, 84th percentiles) for each parameter
            quantiles = [dyfunc.quantile(samples[:, i], [0.16, 0.5, 0.84], weights=weights) 
                         for i in range(samples.shape[1])]
        
            print('\n{:.1f}<z<={:.1f}'.format(z_bins[j], z_bins[j+1]))
            for k, (lo, mid, hi) in enumerate(quantiles):
                print(f"{labels[k]}: {mid:.2f} + {hi-mid:.2f}/-{mid-lo:.2f}")
                err[labels[k]][L_bins[i]].append((lo, hi))
            print()
            
            equal_weight_samples = dyfunc.resample_equal(samples, weights)
            ind = np.random.choice(range(len(equal_weight_samples)), size = 1000, replace = False)
            f_duty_samples[L_bins[i]].append(equal_weight_samples[ind][:,0])
            M_min_samples[L_bins[i]].append(equal_weight_samples[ind][:,1])
            
            f_duty[L_bins[i]].append(quantiles[0][1])
            M_min[L_bins[i]].append(quantiles[1][1]) 

            x = [10**f_duty[L_bins[i]][j], M_min[L_bins[i]][j]]
            if args.ndim == 3:
                sigma_samples[L_bins[i]].append(equal_weight_samples[ind][:,2])
                sigma[L_bins[i]].append(quantiles[2][1]) 
                x.append(sigma[L_bins[i]][j])

            # Get the maximum likelihood sample from dynesty results
            idx_max = np.argmax(results.logl)
            maxL_params[L_bins[i]].append(samples[idx_max])

            # file = '_G{}_rmin{}_rmax{}_nbins{}_{}Lmax{}_N{}_{}_MC_mask'.format(G_hi, rmin, rmax, nbins, pi, L_bins[i], N_jk, cosmo)
                
            # np.save('../results/M_min_{}{}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), M_min[L_bins[i]])
            # np.save('../results/log10(M_min){}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), np.array(err['log10(M_min)'][L_bins[i]]))
    
            # np.save('../results/f_duty{}_{}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), f_duty[L_bins[i]])
            # np.save('../results/log10(f_duty){}_{}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), np.array(err['log10(f_duty)'][L_bins[i]]))
    
            # np.save('../results/zbin_{}'.format(file), 1/np.array(a_zbin[L_bins[i]])-1)
        
    # for i in range(len(L_bins)):
    
    #     print('M_i<={}'.format(L_bins[i]))
    
    #     for j in z_arr
            
            # print('{:.1f}<z<={:.1f}'.format(z_bins[j], z_bins[j+1]))
            # x = 10**f_duty[L_bins[i]][j], M_min[L_bins[i]][j]
            
            # C = np.load('../results/cov_{}_G{}_rmin{}_rmax{}_nbins{}_{}zmin{}zmax{}_Lmax{}_N{}_{}_MC_mask.npy'.format(args.corrfunc, G_hi, rmin, rmax, nbins, pi, z_bins[j], z_bins[j+1], L_bins[i], N_jk, cosmo))
            chi2 = f0_GLS_solve(x, a_zbin[L_bins[i]][j], cf_matter_chi2[L_bins[i]][j], C, cf_chi2[L_bins[i]][j], args.ndim)
            df = nbins-len(x)
        
            print('           chi^2: {:.2f}'.format(chi2))
            print('   reduced chi^2: {:.2f}'.format(chi2/df))
            print('         delta n: {:.2e}'.format(constraint(x, a_zbin[L_bins[i]][j], n_q[L_bins[i]][j], args.ndim)))
            print('condition number: {:.2f}'.format(np.linalg.cond(C)))
            print('             PTE: {:.2f}'.format(stats.chi2.sf(chi2, df)))

        np.save('../results/M_min_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), M_min[L_bins[i]])
        np.save('../results/log10(M_min)_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(err['log10(M_min)'][L_bins[i]]))

        np.save('../results/f_duty_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), f_duty[L_bins[i]])
        np.save('../results/log10(f_duty)_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(err['log10(f_duty)'][L_bins[i]]))

        if args.ndim == 3:
            np.save('../results/sigma_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), sigma[L_bins[i]])
            np.save('../results/sigma_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(err['sigma'][L_bins[i]]))

        np.save('../results/maxL_params_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(maxL_params[L_bins[i]]))
        np.save('../results/zbin{}Lmax{}_{}_MC_mask'.format(file, L_bins[i], cosmo), 1/np.array(a_zbin[L_bins[i]])-1)

if __name__=='__main__':
    main()
    