import numpy as np
import healpy as hp

# import matplotlib
# from matplotlib import pyplot as plt

from scipy import stats
# import matplotlib.gridspec as gridspec
import astropy.cosmology.units as cu
import astropy.units as u
import site
site.addsitedir('/oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts') 
import quaia
# import seaborn as sns
# from matplotlib.colors import LogNorm
# import matplotlib.patches as mpatches
import pyccl as ccl
# from matplotlib.lines import Line2D
# from matplotlib import cm
# from scipy.optimize import minimize
from astropy.cosmology import Planck18
from dynesty import NestedSampler
# from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
# from astropy.coordinates import SkyCoord
from scipy import linalg
from scipy.special import gamma
# import argparse
from scipy import special
# import time
import pickle
from astropy import constants
import mcfit

# Model using step-function HOD
# mask = 0.5
cosmo_planck18 = ccl.Cosmology(
    Omega_c=Planck18.Odm0,
    Omega_b=Planck18.Ob0,
    h=Planck18.h,
    Neff=Planck18.Neff,
    T_CMB=Planck18.Tcmb0.value,
    m_nu=Planck18.m_nu,
    sigma8=0.8102,   # or pull from wherever your actual source specifies Planck18's sigma8
    n_s=0.9665        # same
)
Delta = 200
delta_c = 1.686
mass_def = '200c'
hmf = ccl.halos.MassFuncTinker10(mass_def = mass_def) 
M_max = 16.8 #15

## Define functions
# method = 'pyccl' #'pyccl' #'manual'
b1_E_func_tinker10 = ccl.halos.hbias.tinker10.HaloBiasTinker10(mass_def=mass_def)
n_bins = 201

# Implement halo model
bins = np.logspace(9, M_max, n_bins)
c_m = ccl.halos.concentration.diemer15.ConcentrationDiemer15()

## Limit angular scales to $0.001^\circ\geq \theta \geq 10^\circ$
rmin = 0.001
rmax = 10.0 #80 #+delta_max*offset # 200
nbins = 30 #10 #20 + int(delta_min+delta_max) #90 #20
thetabins = np.logspace(np.log10(rmin), np.log10(rmax), nbins+1)
log_thetamin = -2 #-1
log_thetamax = 0
theta_mask=(quaia.recenter(np.log10(thetabins))>=log_thetamin) & (quaia.recenter(np.log10(thetabins))<=log_thetamax)
# pimax = 80 # 40 # 5000.0 # 40.0\
# dpi = 0.01

### Translate mass definition
dz = np.linspace(0, 4.6)
# Implement Limber approximation: $w(\theta)=\int d z \frac{H(z)}{c}\left(\frac{d N}{d z}\right)^2 \int \frac{d k k}{2 \pi} P_{q q} J_0[k \theta \chi(z)]=\int d z \frac{H(z)}{c}\left(\frac{d N}{d z}\right)^2 w_p(\theta \chi(z))$
# rpbins = [(thetabins*u.deg).to(u.rad)*quaia.comoving_dist(z, units = 'Mpc') for z in dz]
# P_mm = [ccl.power.linear_matter_power(cosmo_planck18, k_log, 1/(1+z)) for z in dz]
k_log = np.logspace(-2.5, 3) #-4, 2)
da = 1/(1+dz)
P_mm = ccl.power.linear_matter_power(cosmo_planck18, k_log, da)
rpbins = (thetabins*u.deg).to(u.rad)[:, None]*quaia.comoving_dist(dz, units = 'Mpc')

## Pick cosmology
cosmo = 'Planck'
N_jk = 50 #40 #20 #65 #50 #30 #15 #80
NFW_200c = ccl.halos.profiles.nfw.HaloProfileNFW(mass_def = '200c', concentration = c_m)
# cosmo_planck18 = ccl.Cosmology(Omega_c=0.265, Omega_b=Planck18.Ob0, h=Planck18.h, Neff = Planck18.Neff, T_CMB = Planck18.Tcmb0.value, m_nu = Planck18.m_nu, sigma8=0.8, n_s=0.95)
hartlap = (N_jk-nbins-2)/(N_jk-1)

## Define functions
# def a_median(tab_datahi_mask_zbin0, key_zbin0 = ''):
#     return 1/(1+np.median(tab_datahi_mask_zbin0['redshift_quaia'+key_zbin0]))

# def b_qso(xi_qq, xi_mm, f):
#     return np.sqrt(xi_qq/xi_mm-4*f**2/45)-f/3

# def nu(M, delta_c = delta_c, a = 1, cosmo = cosmo_planck18):
#     sigma = cosmo.sigmaM(M, a = a)
#     return delta_c/sigma

def n_g_theory(bins, n_total, a = 1, cosmo = cosmo_planck18, hmf = hmf, recenter = True, units = '(h/Mpc)^3'):

    if recenter == True:
        bins = quaia.recenter(bins)
    if units == '(h/Mpc)^3':
        nm = hmf(cosmo, bins, a)/cosmo['h']**3*(cu.littleh/u.Mpc)**3
    else:
        nm = hmf(cosmo, bins, a)*(1/u.Mpc)**3
    n_g = nm/(bins*np.log(10))*n_total
    n_g[np.isnan(n_g)]=0
    
    return n_g, nm # (units of Mpc/h)^3

def bias(bins, n_m, a = 1, cosmo = cosmo_planck18, mod = False, recenter = True, hmf = hmf, #method = method,
         b1_E_func_tinker10 = b1_E_func_tinker10, b1_E_tinker = None, n_g_manual = None, calc_bias = True, units = '(h/Mpc)^3'):

    n_g, nm = n_g_theory(bins, n_m, a = a, cosmo = cosmo, recenter = recenter, hmf = hmf, units = units)
    if recenter == True:
        print('recentering')
        bins = quaia.recenter(bins)
        
    if n_g_manual is None:
        n_g = np.trapz(n_g, bins)
        n_g_manual = n_g
        
    if calc_bias == False:
        return n_g

    else:
        
        dn_dm = nm/(bins*np.log(10))
    
        if b1_E_tinker is None:
            
            # b1_E_tinker = b1_E_func(bins, a = a, cosmo = cosmo, mod = mod, b1_E_func_tinker10 = b1_E_func_tinker10) #method = method,
            b1_E_tinker = b1_E_func_tinker10(cosmo, bins, a)
            
        n_m[~np.isfinite(n_m)]=0 
        integrand = dn_dm*n_m
        integrand[np.isnan(integrand)]=0

        return bias_function(n_g_manual, integrand, b1_E_tinker, bins, recenter = recenter), n_g

# def b1_E_func(bins, a = 1, cosmo = cosmo_planck18, Delta = Delta, delta_c = delta_c, mod = False, b1_E_func_tinker10 = b1_E_func_tinker10): #nu #method = method,

#     if method == 'pyccl':

#         return b1_E_func_tinker10(cosmo, bins, a)

#     else:

#         Nu = nu(bins, a = a, cosmo = cosmo)
        
#         if mod == False:
#             y = np.log10(Delta)
#             A =  1.0+0.24*y*np.exp(-(4/y)**4)
#             a = 0.44*y-0.88
#             B = 0.183
#             b =  1.5
#             C = 0.019+0.107*y+0.19*np.exp(-(4/y)**4)
#             c = 2.4
            
#         else:
#             A =  1.0
#             a = 0.0906
#             B = -4.5002
#             b =  2.1419
#             C = 4.9148
#             c = 2.1419
        
#         return 1-A*Nu**a/(Nu**a+delta_c**a)+B*Nu**b+C*Nu**c # Jose et al
    
def bias_function(n_g, integrand, b, bins, recenter = True):
    
    if recenter == True:
        print('recentering')
        bins = quaia.recenter(bins)
        
    return 1/n_g*np.trapz(integrand*b, bins)

# def cf_2h(bins, n_m, a, cosmo, cf_matter, recenter = True, method = method, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10):
    
#     b_q, n_g = bias(bins, n_m, a, cosmo, recenter = recenter, method = method, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10)
    
#     return b_q**2*cf_matter, b_q

# def M(bins, n_m, a = 1, cosmo = cosmo_planck18, mod = False, recenter = True, hmf = hmf, units = '(h/Mpc)^3'):

#     n_g, nm = n_g_theory(bins, n_m, a = a, cosmo = cosmo, recenter = recenter, hmf = hmf, units = units)
#     if recenter == True:
#         n_g = np.trapz(n_g, quaia.recenter(bins))
#         dm = np.diff(bins) 
#         dn_dm = nm/(quaia.recenter(bins)*np.log(10))
#     else:
#         n_g = np.trapz(n_g, bins)
#         dn_dm = nm/(bins*np.log(10))

#     n_m[~np.isfinite(n_m)]=0 
#     integrand = dn_dm*n_m
#     integrand[np.isnan(integrand)]=0
    
#     if recenter == True:
#         return bias_function(n_g, integrand, quaia.recenter(bins), bins, recenter = recenter), n_g
#     else:
#         return bias_function(n_g, integrand, bins, bins, recenter = recenter), n_g

# def HOD(x, ndim, bins_m = bins):

#     # n_m=np.zeros(np.shape(quaia.recenter(bins)))
#     # n_m[quaia.recenter(bins)>=M_min]=fduty
#     if ndim == 2:
#         # bins_m, fduty = x
#         fduty, M_min = x
#         bins_m = np.logspace(M_min, M_max, 79)
#         n_m=np.ones(np.shape(bins_m))*10**fduty
#     else:
#         # fduty, M_min, M, sigma = x
#         fduty, M_min, sigma = x
#         n_m = 10**fduty/2*special.erfc(np.log(10**M_min/bins_m)/(np.sqrt(2)*np.log(10)*sigma))
        
#     return n_m, bins_m
    
# def HOD_step(bins_m, fduty): #M_max = M_max):

#     # bins_m = np.logspace(M_min, M_max, 79)
#     n_m=np.ones(np.shape(bins_m))*10**fduty
#     return n_m

# def HOD_step(bins, fduty): #, M_min, recenter = False):

#     return np.ones(np.shape(bins))*fduty

def HOD_scatter(bins, M_min, fduty, sigma):
    
    return fduty/2*special.erfc(np.log(10**M_min/bins)/(np.sqrt(2)*np.log(10)*sigma))

def N_c(bins_m, M_min, f_duty, M_max = M_max, sigma = 0):#, step = False): #, recenter = False): #bins_m

    n_m = HOD_scatter(bins_m, M_min, f_duty, sigma)
    
    # if np.all(bins_m[0]==10**M_min) and np.all(sigma==0): # for step function HOD

    #     return HOD_step(bins_m, f_duty) # M_min, f_duty, recenter = recenter)

    # if np.all(sigma==0):

    #     n_m[bins_m==10**M_min]=f_duty

    # else:
    # try:
        
    #     n_m[(bins_m==10**M_min) & (sigma[None, :]==0)]=f_duty/2 # for multiple values

    # except:

    # if sigma !=0:

    if sigma == 0:
        
        n_m[bins_m==10**M_min]=f_duty/2 # assuming im just feeding in the scatter HOD and correcting whenever it samples sigma = 0 
        
    return n_m
    
# def HOD_erfc(fduty, M_min, M, sigma):
    
#     return 10**fduty/2*special.erfc(np.log(10**M_min/M)/(np.sqrt(2)*np.log(10)*sigma))

# def n_DM(cosmo, a, M_min, hmf = hmf, M_max = M_max, N = 79): #16.8
    
#     bins = np.logspace(M_min, M_max, N)
#     nm = hmf(cosmo, bins, a)/cosmo['h']**3 #/u.Mpc**3 # convert to h/Mpc**3
#     dn_dm = nm/(bins*np.log(10))
    
#     return np.trapz(dn_dm, bins)

# def n_DM(cosmo, a, M_min, hmf = hmf[mass_def], M_max = M_max, N = 79): #16.8
    
#     bins = np.logspace(M_min, M_max, N)
#     nm = hmf(cosmo, bins, a)/cosmo['h']**3 #/u.Mpc**3 # convert to h/Mpc**3
#     dn_dm = nm/(bins*np.log(10))
    
#     return np.trapz(dn_dm, bins)

def N_s(bins_m, fduty, M0, M1, alpha, exp = False, step = True): #bins_m

    if exp == False:

        n_m = fduty*((bins_m-M0)/M1)**alpha
        
        if step == True:
            n_m[bins_m<M0]=0
            
    else:
        # n_m = (bins_m/M1)**alpha*np.exp(-M0/bins_m)
        n_m = N_s_exp(bins_m, M1, alpha, M0, fduty)
        
    return n_m

def N_s_exp(M, M1, alpha, Mcut, fduty):

    return fduty*(M/M1)**alpha*np.exp(-Mcut/M)

# def NFW(bins_m, a, k = k_log, Delta = Delta, cosmo = cosmo_planck18, c_m = c_m, mass_def = '200c'):

#     c = c_m(cosmo, bins_m, a) # assumes 200c mass definition
#     f = 1/(np.log(1+c)-c/(1+c))

#     rho_bar = ccl.background.rho_x(cosmo, a, 'critical', is_comoving = True)  #'matter'  # Rvir = (3*bins_m*u.Msun/(4*np.pi*Delta*rho_bar))**(1/3) # Mpc
#     Rvir = (3*bins_m/(4*np.pi*Delta*rho_bar))**(1/3) # Mpc
#     kappa = k[:, None]*Rvir/c
#     Si, Ci = sici(kappa*(1+c))
#     Si_kappa, Ci_kappa = sici(kappa)
    
#     u = f*(np.sin(kappa)*(Si-Si_kappa)+np.cos(kappa)*(Ci-Ci_kappa)-np.sin(kappa*c)/(kappa*(1+c)))
#     u[k==0]=1
    
#     return u

# def HOD(bins, M_min, fduty, sigma, a, M1, alpha, method = method, NFW_200c = NFW_200c, k = k_log, cosmo = cosmo_planck18): #, recenter = False):
    
#     if method == 'manual':
        
#         nfw = NFW(bins, a, k = k)

#     else:

#         nfw = NFW_200c.fourier(cosmo, k, bins, a).T/bins

#     return N_c(bins, M_min, 10**fduty, sigma = sigma)+nfw*N_s(bins, 10**fduty, 10**M_min, 10**M1, alpha) #, recenter = recenter

def HOD_2h(bins, M_min, fduty, sigma, a, M1, alpha, NFW_200c = NFW_200c, k = k_log, cosmo = cosmo_planck18): #, recenter = False): # method = method, 
    
    # if method == 'manual':
        
        # nfw = NFW(bins, a, k = k)

    # else:

    nfw = NFW_200c.fourier(cosmo, k, bins, a).T/bins

    n_c = N_c(bins, M_min, 10**fduty, sigma = sigma)
    n_s = N_s(bins, 10**fduty, 10**M_min, 10**M1, alpha)
    
    return n_c+nfw*n_s, n_c, n_s #, recenter = recenter

def P_1h(bins, a, n_g, params_n_m = None, params_HOD = None, k = k_log, cosmo = cosmo_planck18, recenter = False, hmf = hmf, #n_c, n_s #method = method, 
         b1_E_func_tinker10 = b1_E_func_tinker10, NFW_200c = NFW_200c, damping = False, k_star = 1e-2/Planck18.h, units = '1/Mpc^3'): #[mass_def]

    if params_n_m is not None:
        
        n_c, n_s = params_n_m

    else:
        
        fduty, M_min, sigma, M1, alpha = params_HOD
        n_c = N_c(bins, M_min, 10**fduty, sigma = sigma)
        n_s = N_s(bins, 10**fduty, 10**M_min, 10**M1, alpha)

    # if method == 'manual':
    #     nfw = NFW(bins, a, k = k)

    # else:
        
    nfw = NFW_200c.fourier(cosmo, k, bins, a).T/bins
    
    n_m = 2*nfw*n_c*n_s+nfw**2*n_s**2 #(n_s-1) # assume poisson distribution for satellites # should i implement central satellite condition? i think no
    P_qq_1h, _ = bias(bins, n_m, a, cosmo, recenter = recenter, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, #method = method,
                           b1_E_tinker = np.ones(np.shape(bins)), n_g_manual = n_g, units = units)

    if damping == True:
        P_qq_1h *= (k/k_star)**4/(1+(k/k_star)**4)

    return P_qq_1h/n_g #**2
    
# def P_1h(bins, a, n_c, n_s, k = k_log, cosmo = cosmo_planck18, recenter = False, method = method, hmf = hmf[mass_def], 
#          b1_E_func_tinker10 = b1_E_func_tinker10[mass_def], NFW_200c = NFW_200c):

#     if method == 'manual':
#         nfw = NFW(bins, a, k = k)*n_s

#     else:
        
#         nfw = NFW_200c.fourier(cosmo, k, bins, a).T/bins*n_s
                
#     n_m = 2*nfw*n_c*n_s+nfw**2*n_s*(n_s-1)
#     P_qq_1h, n_g = bias(bins, n_m, a, cosmo, recenter = recenter, method = method, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, 
#                            b1_E_tinker = np.ones(np.shape(bins)))

#     return P_qq_1h/n_g

# def wtheta_2h_optimized(bins, H0_c_dN_dz, params_n_m = None, params_HOD = None, # P_mm = P_mm, 
#                         # wp = False, 
#                         thetabins = thetabins, dz = dz, k = k_log, recenter = False, #n_m
#                         theta_mask = False, cosmo = cosmo_planck18, # rpbins = rpbins, np.full_like(quaia.recenter(thetabins), True, dtype = bool)
#                         hankel = mcfit.Hankel(k_log, lowring = True), recenter_thetabins = False, calc_bias = False, wtheta_mm = False, method = method, 
#                         hmf = hmf[mass_def], b1_E_func_tinker10 = b1_E_func_tinker10[mass_def], NFW_200c = NFW_200c, one_halo = False):

#     # compute for each z... for statement?
#     b_q, k_integral = [], []
#     for z in dz:
#         a = 1/(1+z)
#         P_mm = ccl.power.linear_matter_power(cosmo, k, a)
#         rpbins = (thetabins*u.deg).to(u.rad)*quaia.comoving_dist(z, units = 'Mpc')
        
#         if params_n_m is not None:
            
#             n_c, n_s = params_n_m

#             # print(np.shape(bins))
#             if method == 'manual':
#                 n_m = n_c + NFW(bins, a, k = k)*n_s

#             else:
#                 n_m = n_c + NFW_200c.fourier(cosmo, k, bins, a).T/bins*n_s
#         else:
            
#             fduty, M_min, sigma, M1, alpha = params_HOD
#             n_m = HOD(bins, M_min, fduty, sigma, a, M1, alpha, method = method, cosmo = cosmo, k = k) #, recenter = recenter)
#         if wtheta_mm == True:
#             b_z = 1
#         else:
#             b_z, n_z = bias(bins, n_m, a, cosmo, recenter = recenter, method = method, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10)
#         # if wp: # need to implement real space version, do centrals only for now
#         #     wp_z, _ = wp_2h(bins, n_c, a, cosmo, False, rpbins, recenter = recenter_thetabins, method = method, hmf = hmf, 
#         #                     b1_E_func_tinker10 = b1_E_func_tinker10) #n_m #recenter
#         else:
            
#             if one_halo == False:
                
#                 P_qq = b_z**2*P_mm
                
#             else:
                
#                 P_qq_2h = b_z**2*P_mm
#                 P_qq_1h = P_1h(bins, a, n_c, n_s, k = k, cosmo = cosmo, recenter = recenter, method = method, hmf = hmf, 
#                                     b1_E_func_tinker10 = b1_E_func_tinker10, NFW_200c = NFW_200c)
#                 P_qq = P_qq_2h + P_qq_1h
                
#             rpbins_k, wp_k = hankel(P_qq/(2*np.pi))
#             if recenter_thetabins == True:
#                 rp = quaia.recenter(rpbins.value)

#             else:
#                 rp = rpbins.value
#             if isinstance(theta_mask, bool):
#                 theta_mask = np.full_like(rp, True, dtype = bool)
#             wp_z = np.interp(rp[theta_mask], rpbins_k, wp_k)
            
#         k_integral.append(wp_z) # as a function of rp
#         b_q.append(b_z)
        
#     z_integrand = [H0_c_dN_dz[i]*k_integral[i] for i in range(len(dz))]
#     if calc_bias == False:
    
#         return np.trapz(z_integrand, dz, axis = 0)#, np.trapz(b_q*dN_dz, z_array) #np.trapz(b_q*np.array(n_g)*dN_dz, z_array)/np.trapz(n_g*dN_dz, z_array)

#     else:
#         return np.trapz(z_integrand, dz, axis = 0), np.trapz(b_q*dN_dz, dz)

def wtheta_2h_dynesty(bins, H0_c_dN_dz, params_HOD, thetabins = thetabins, k = k_log, recenter = False, theta_mask = theta_mask, cosmo = cosmo_planck18, hankel = mcfit.Hankel(k_log, lowring = True), recenter_thetabins = True, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, #method = method
                      NFW_200c = NFW_200c, P_mm = P_mm, da = da, rpbins = rpbins, dz = dz): #, k_integral = np.empty(np.shape(rpbins.T)), z_integrand = np.empty(np.shape(rpbins.T))):

    # compute for each z... for statement?
    b_q, P_qq = [], [] #k_integral #, []
    fduty, M_min, sigma, M1, alpha = params_HOD
    for i, a in enumerate(da):
        
        n_m, n_c, n_s = HOD_2h(bins, M_min, fduty, sigma, a, M1, alpha, cosmo = cosmo, k = k) #, recenter = recenter)
    
        n_g = bias(bins, n_c+n_s, a, cosmo, recenter = recenter, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, 
                               b1_E_tinker = np.ones(np.shape(bins)), calc_bias = False, units = '1/Mpc^3')
        
        b_z, _ = bias(bins, n_m, a, cosmo, recenter = recenter, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, n_g_manual = n_g, units = '1/Mpc^3')
            
        P_qq_2h = b_z**2*P_mm[i]*u.Mpc**3
        P_qq_1h = P_1h(bins, a, n_g, params_HOD = params_HOD, k = k, cosmo = cosmo, recenter = recenter, hmf = hmf, 
                            b1_E_func_tinker10 = b1_E_func_tinker10, NFW_200c = NFW_200c) #n_c, n_s
        P_qq.append(P_qq_2h + P_qq_1h)
        b_q.append(b_z)
    
    rpbins_k, wp_k = hankel(np.array(P_qq)/(2*np.pi))
    
    if recenter_thetabins == True:
        rp = quaia.recenter(rpbins.value)

    else:
        rp = rpbins.value

    if isinstance(theta_mask, bool):
        theta_mask = np.full_like(thetabins, True, dtype = bool)

    k_integral = [np.interp(rp[theta_mask][:,i], rpbins_k, wp_k[i]) for i in range(len(dz))]
    z_integrand = H0_c_dN_dz[:, None] * k_integral

    return np.trapz(z_integrand, dz, axis = 0)

def log_det(C, ix_grid):  
    
    # check if cov is hermitian and positive definite
    if linalg.ishermitian(C[ix_grid]) and np.all(np.linalg.eigvals(C[ix_grid]) > 0):
        C_factor = linalg.cho_factor(C[ix_grid])
        log_det_C = 2 * np.sum(np.log(np.diag(C_factor[0])))
        solve = lambda rhs: linalg.cho_solve(C_factor, rhs)
    else:
        print('using LU decomposition')
        C_factor = linalg.lu_factor(C[ix_grid])
        log_det_C = np.sum(np.log(np.abs(np.diag(C_factor[0])))) # does abs value make sense
        solve = lambda rhs: linalg.lu_solve(C_factor, rhs)

    return log_det_C, solve
    
# def loglike(x, a, cf_matter, log_det_C, cf, n_q, nbins, solve, penalty, ndim, cosmo = cosmo_planck18, error = 0.05, M_max = M_max, hartlap = hartlap, N = N_jk):
    
#     # if HOD == '_step':
#     # if ndim == 2:

#     #     fduty, M_min = x
#     #     bins_m = np.logspace(M_min, M_max, 79)
#     #     # n_m=np.ones(np.shape(bins_m))*10**fduty
#     #     # n_m = HOD(bins_m, fduty)
#     #     params = bins_m, fduty
        
#     # else:

#     #     fduty, M_min, sigma = x
#     #     params = 
        
#     n_m, bins_m = HOD(x, ndim)
#     n_q_theory = M(bins_m, n_m, a = a, cosmo = cosmo, recenter = False)[1]
#     if penalty == '_soft':
#         lpenalty = 0.5*((n_q_theory-n_q)/(error*n_q))**2
#     else:
#         lpenalty = np.abs(n_q_theory-n_q)/n_q
#         if lpenalty > error:
#             return -np.inf
        
#     cf_2h_zbin, _ = cf_2h(bins_m, n_m, a, cosmo, cf_matter, recenter = False)
#     delta = cf-cf_2h_zbin
                
#     c_p = gamma(N/2)/((np.pi*(N-1))**(nbins/2)*gamma((N-nbins)/2))
#     lnorm = np.log(c_p)-0.5*log_det_C  
    
#     return -N*0.5 * np.log(1+(delta.T @ solve(delta))/(N-1)) + lnorm - lpenalty

def loglike(x, log_det_C, wp, nbins, solve, dN_dz, theta_mask = theta_mask, M_max = M_max, N = N_jk, fduty = 0, recenter_thetabins = True, n_bins = n_bins, bins_m = bins):
    
    M_min, Alpha, sigma = x
        
    wp_2h_zbin = wtheta_2h_dynesty(bins_m, dN_dz, params_HOD = [0, M_min, sigma, np.log10(12*10**M_min), Alpha])

    delta = wp-wp_2h_zbin
    
    c_p = gamma(N/2)/((np.pi*(N-1))**(nbins/2)*gamma((N-nbins)/2))
    lnorm = np.log(c_p)-0.5*log_det_C  

    return -N*0.5 * np.log(1+(delta.T @ solve(delta))/(N-1)) + lnorm #- penalty
    
# def ptform_step(u):
#     """Transforms the uniform random variables `u ~ Unif[0., 1.)`
#     to the parameters of interest."""

#     x = np.array(u)  # copy u 

#     # Uniform from 0 to 1
#     x[0] =  6* u[0] - 6 # scale and shift to [-6., 0.)

#     # Uniform from 9 to 15
#     x[1] = 6 * u[1] + 9  # scale and shift to [9., 15.)

#     return x

# def ptform(u, ndim):
#     """Transforms the uniform random variables `u ~ Unif[0., 1.)`
#     to the parameters of interest."""

#     x = np.array(u)  # copy u 

#     # Uniform from -6 to 0
#     x[0] = 6* u[0] - 6 # scale and shift to [-6., 0.)

#     # Uniform from 9 to 15
#     x[1] = 6 * u[1] + 9  # scale and shift to [9., 15.)

#     if ndim == 3:
        
#         # Uniform from 0 to 1.4
#         x[2] = 1.4 * u[2]  # scale to [0., 1.4.)

#     return x

def ptform(u):
    """Transforms the uniform random variables `u ~ Unif[0., 1.)`
    to the parameters of interest."""

    x = np.array(u)  # copy u 

    # Uniform from 9 to 16
    x[0] = 7 * u[0] + 9  # scale and shift to [9., 16.)

    # Uniform from 0 to 4
    x[1] = 4 * u[1]  # scale and shift to [0, 4.)

    # Uniform from 1 to 60
    x[2] = 1.4 * u[2]   # scale and shift to [0, 1.4.)

    return x
    
# def f0_GLS_solve(x, a, cf_matter, C, cf, ndim, cosmo = cosmo_planck18, M_max = M_max, hartlap = hartlap, bins_m = bins):  #quaia):
        
#     # if HOD == '_step':

#     #     fduty, M_min = x
#     #     bins_m = np.logspace(M_min, M_max, 79)
#     #     # n_m=np.ones(np.shape(bins_m))*10**fduty
#     #     n_m = HOD(bins_m, fduty)
        
#     # else:

#     #     fduty, M_min, sigma = x
#     #     n_m = HOD_scatter(fduty, M_min, M, sigma)

#     n_m, bins_m = HOD(x, ndim)
#     cf_2h_zbin, _ = cf_2h(bins_m, n_m, a, cosmo, cf_matter, recenter = False)
#     delta = cf-cf_2h_zbin
        
#     return delta.T @ (hartlap*np.linalg.solve(C, delta))

# def constraint(x, a, n_q, ndim, cosmo = cosmo_planck18, error = 0.05, M_max = M_max, bins_m = bins):
    
#     # if HOD == '_step':

#     #     fduty, M_min = x
#     #     bins_m = np.logspace(M_min, M_max, 79)
#     #     # n_m=np.ones(np.shape(bins_m))*10**fduty
#     #     n_m = HOD(bins_m, fduty)
        
#     # else:

#     #     fduty, M_min, sigma = x
#     #     n_m = HOD_scatter(fduty, M_min, M, sigma)

#     n_m, bins_m = HOD(x, ndim)
#     n_q_theory = M(bins_m, n_m, a = a, cosmo = cosmo, recenter = False)[1]

#     return error-np.abs(n_q_theory-n_q)/n_q
    
def main():

    # parser = argparse.ArgumentParser()
    # parser.add_argument("corrfunc", help="two-point statistic", type = str)
    # parser.add_argument("penalty", help="penalty for number density constraint", type = str)
    # # parser.add_argument("HOD", help="HOD type", type = str)
    # parser.add_argument("ndim", help="number of HOD parameters", type = int)
    # args = parser.parse_args()

    ## Define parameters
    G_hi = 20.5
    nthreads = 18
    z_bins = np.array([0.0, 1.0, 2.0, 3.0, 4.6])
    z_array = range(len(z_bins)-1)
    NSIDE = 64

    # Visualize the luminosity threshold cuts and redshift bins
    L_bins = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])[:0:-1]
    b = 30
    fac_rand = 25
    # pibins = np.arange(0, pimax+1, dpi)
    
    ## Calculate number density as $n_q \approx \int_0^{\infty} d M \frac{d n}{d M}\langle N(M)\rangle$
    # tab_datahi_mask_zbins, selfunc_hi_bins = {L: [] for L in L_bins}, {L: [] for L in L_bins}

    # for i in range(len(L_bins)):
            
    #     for j in z_array:
            
    #         tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [j, i], True, ['z', 'L'], 
    #                                                                                          bins = [z_bins, L_bins], tab_gcat_type = 'data', 
    #                                                                                           method = ['minmax', 'max'], fac_rand = fac_rand, 
    #                                                                                           mask = mask*100, percentile = True, n_bins = [None, None])
    #         tab_datahi_mask_zbins[L_bins[i]].append(tab_datahi_mask_zbin0['redshift_quaia'])
    #         selfunc_hi_bins[L_bins[i]].append(selfunc_hi_bin0)

    # LMC = SkyCoord(['05 23 34.6 -69 45 22'], unit=(u.hourangle, u.deg)) #https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Large+Magellanic+Cloud
    # LMC_radius = 9*u.deg
    # SMC = SkyCoord(['00 52 38.0 -72 48 01'], unit=(u.hourangle, u.deg)) # https://simbad.u-strasbg.fr/simbad/sim-id?Ident=small+Magellanic+Cloud
    # SMC_radius = 5*u.deg
    # NPIX = hp.nside2npix(NSIDE)
    # theta, phi = hp.pix2ang(NSIDE, range(NPIX), lonlat=True)
    
    # c = SkyCoord(ra=theta*u.degree, dec=phi*u.degree)
    # LMC_idx, _, _, _ = SkyCoord.search_around_sky(LMC, c, LMC_radius)
    # SMC_idx, _, _, _ = SkyCoord.search_around_sky(SMC, c, SMC_radius)
    # MC_idx = np.append(LMC_idx,SMC_idx)
    # MC_mask = np.full_like(range(NPIX), True, dtype = bool)
    # MC_mask[MC_idx] = False

    # selfunc_mask = np.abs(c.galactic.b.value)>=b
    
    # n_q, n_q_err, V_eff, N_q = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}

    # for i in range(len(L_bins)):
            
    #     for j in z_array:   

    #         tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [j, i], True, ['z', 'L'], 
    #                                                                                          bins = [z_bins, L_bins], tab_gcat_type = 'data', 
    #                                                                                           method = ['minmax', 'max'], fac_rand = fac_rand, b = b, 
    #                                                                                                                   mask_type = 'b',
    #                                                                                                                   n_bins = [None, None])
    
    #         N = len(tab_datahi_mask_zbin0)
    #         # cp = np.percentile(selfunc_hi_bins[L_bins[i]][j][MC_mask], mask*100)
    #         # selfunc_mask = selfunc_hi_bins[L_bins[i]][j]>=cp
    #         f_sky = 1/(4*np.pi)*np.sum(selfunc_hi_bin0[selfunc_mask & MC_mask]*hp.nside2pixarea(NSIDE))
    #         # print(f_sky) # f_sky)
    #         V = f_sky*(Planck18.comoving_volume(z_bins[j+1])-Planck18.comoving_volume(z_bins[j]))*Planck18.h**3/u.h**3
    #         N_q[L_bins[i]].append(N)
    #         n_q[L_bins[i]].append(N/V.value)
    #         n_q_err[L_bins[i]].append(np.sqrt(N)/V.value) # Poisson error
    #         V_eff[L_bins[i]].append(V)
            
    # MC_mask[MC_idx] = False

    ## Limit scales to $30\geq s \geq 80h^{-1}$ Mpc
    # Create the bins array
    # rbins_chi2 = np.linspace(rmin, rmax, nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202

    # Define the dimensionality of our problem.
    ndim = 3
    diag = False #True # False
    ix_grid = np.ix_(theta_mask, theta_mask)
    
    nlive = 500
    dlogz = 0.01
    nbins_mask = np.sum(theta_mask)
    quantile = [0.16, 0.5, 0.84]
    
    # rstate=np.random.default_rng(0)
    thetamin, thetamax = 10**log_thetamin, 10**log_thetamax

    # cf_chi2, cf_matter_chi2 = {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # yerr_cf_chi2 = {L: [] for L in L_bins}
    # , b_cf_chi2 = {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # b_cf_q_chi2 = {L: [] for L in L_bins}
    # a_zbin = {L: [] for L in L_bins}
    # rpavg_chi2 = {L: [] for L in L_bins}
    # H0 = Planck18.H0 # cosmo_planck18['H0'] * u.km/u.s/u.Mpc
    # wp_chi2 = {L: [] for L in L_bins}

    ## Across all $z$
    # n_q, n_q_err, V_eff, N_q = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # dN_dz_norm, H0_c_dN_dz_norm = {L: [] for L in L_bins}, {L: [] for L in L_bins}
    err_labels = ['log10(f_duty)', 'log10(M_min)', 'alpha', 'Sigma']
    err = {i: {L: [] for L in L_bins} for i in err_labels}
    # best_fit_params = {L: [] for L in L_bins}
    # equal_weight_samples = {L: [] for L in L_bins}
    # results_Mmin = {L: [] for L in L_bins}
    sample = {err_label: {L: [] for L in L_bins} for err_label in err_labels}
    median = {err_label: {L: [] for L in L_bins} for err_label in err_labels}

    for i in range(len(L_bins)):

        print('M_i<={}'.format(L_bins[i]))

        for j in z_array:   

            tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [j, i], True, ['z', 'L'], 
                                                                                             bins = [z_bins, L_bins], tab_gcat_type = 'data', 
                                                                                              method = ['minmax', 'max'], fac_rand = fac_rand, b = b, 
                                                                                                                      mask_type = 'b',
                                                                                                                      n_bins = [None, None])
    
            # N = len(tab_datahi_mask_zbin0)
            # cp = np.percentile(selfunc_hi_bins[L_bins[i]][j][MC_mask], mask*100)
            # selfunc_mask = selfunc_hi_bins[L_bins[i]][j]>=cp
            # f_sky = 1/(4*np.pi)*np.sum(selfunc_hi_bin0[selfunc_mask & MC_mask]*hp.nside2pixarea(NSIDE))
            # print(f_sky) # f_sky)
            # V = f_sky*(Planck18.comoving_volume(z_bins[j+1])-Planck18.comoving_volume(z_bins[j]))*Planck18.h**3/u.h**3
            N_q = len(tab_datahi_mask_zbin0)
            # n_q[L_bins[i]].append(N_q/V.value)
            # n_q_err[L_bins[i]].append(np.sqrt(N_q)/V.value) # Poisson error
            # V_eff[L_bins[i]].append(V)
            
    # for i in range(len(L_bins)):
                
    #     for j in z_array:

            wp_zbin_i = np.load('../results/wtheta_G{}_rmin{}_rmax{}_nbins{}_zmin{}zmax{}_Lmax{}_{}_MC_b{}.npy'.format(G_hi, rmin, rmax, nbins, z_bins[j], 
                                                                                                                       z_bins[j+1], L_bins[i], cosmo, b))
            # a = float(np.load('../results/a_G{}_zmin{}zmax{}_Lmax{}_{}_MC_b{}.npy'.format(G_hi, z_bins[j], z_bins[j+1], L_bins[i], cosmo, b)))
          
            # wp_chi2[L_bins[i]].append(wp_zbin_i)
    
            # # scale factor
            # a_zbin[L_bins[i]].append(a)
            
            # tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [j, i], True, ['z', 'L'], 
            #                                                                                  bins = [z_bins, L_bins], tab_gcat_type = 'data', 
            #                                                                                   method = ['minmax', 'max'], fac_rand = fac_rand, 
            #                                                                                   mask = mask*100, percentile = True, n_bins = [None, None])
    
            # # scale factor
            # a = a_median(tab_datahi_mask_zbin0)
            # a_zbin[L_bins[i]].append(a)
            
            # # quasar correlation function
            # if args.corrfunc == 'xi':
                
            #     cf, qq, qr = quaia.xi_s(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, 
            #                rbins = rbins_chi2, Om0 = Planck18.Om0, h = Planck18.h, error = True)
            #     cf_chi2[L_bins[i]].append(cf)
    
            #     # matter correlation function
            #     rbins_Mpc = (quaia.recenter(rbins_chi2)*u.Mpc/cu.littleh).to(u.Mpc, cu.with_H0(H0))
            #     cf_matter = ccl.correlations.correlation_3d(cosmo_planck18, r=rbins_Mpc, a=a)
            #     pi = ''
                
            # else:
            #     # quasar correlation function
            #     wp_zbin_i, rpavg_zbin_i = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, 
            #                                                                   rbins = rbins_chi2, nbins = nbins, pimax = pimax, Om0 = Planck18.Om0, h = Planck18.h)
            #     cf_chi2[L_bins[i]].append(wp_zbin_i)
            #     rpavg_chi2[L_bins[i]].append(rpavg_zbin_i)
            
            #     # matter correlation function   
            #     s = np.sqrt(np.sum(np.meshgrid(np.array(rpavg_zbin_i)**2, quaia.recenter(pibins)**2), axis = 0))
            #     s_Mpc = (s*u.Mpc/cu.littleh).to(u.Mpc, cu.with_H0(H0))
                
            #     xirppi_zbin_i = ccl.correlations.correlation_3d(cosmo_planck18, r=np.ravel(s_Mpc), a=a_median(tab_datahi_mask_zbin0))
            #     cf_matter = 2*np.trapz(np.reshape(xirppi_zbin_i, np.shape(s_Mpc)), quaia.recenter(pibins), axis = 0)
            #     pi = 'pimax{}_'.format(pimax)
        
            # cf_matter_chi2[L_bins[i]].append(cf_matter)
        
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

    # for i in range(len(L_bins)):
    #     for j in z_array:
        
            # tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(
            #     G_hi, [j, i], True, ['z', 'L'], bins = [z_bins, L_bins], tab_gcat_type = 'data', method = ['minmax', 'max'], fac_rand = fac_rand, 
            #     n_bins = [None, None], b = b, mask_type = 'b', verbose = False)
            
            dN_dz = np.sum(stats.norm.pdf(dz[None, :], loc=tab_datahi_mask_zbin0['redshift_quaia'][:, None], scale=tab_datahi_mask_zbin0['redshift_quaia_err'][:, None]), axis = 0)
            dN_dz /= np.sum(stats.norm.cdf(dz[-1], loc=tab_datahi_mask_zbin0['redshift_quaia'], scale=tab_datahi_mask_zbin0['redshift_quaia_err']))
            
            # dN_dz_norm[L_bins[i]].append(dN_dz)
            H0_c_dN_dz_norm = ccl.background.h_over_h0(cosmo_planck18, da)*Planck18.H0.to(1/u.s)/constants.c.to(u.Mpc/u.s)*dN_dz**2

    # for i in range(len(L_bins)):
    
        # for j in z_array:

            file = '_G{}_rmin{}_rmax{}_nbins{}_zmin{}zmax{}_Lmax{}_N{}_{}_MC_b{}'.format(G_hi, thetamin, thetamax, nbins, z_bins[j], z_bins[j+1], L_bins[i], N_jk, cosmo, b)
            
            try:
    
                for k, err_label in enumerate(err_labels[1:]):
                    
                    sample[err_label][L_bins[i]].append(np.load('../results/erfc_samples_{}{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                        file, ndim, dlogz, nlive)))
                    median[err_label][L_bins[i]].append(np.load('../results/erfc_median_{}{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                            file, ndim, dlogz, nlive)))
                    err[err_label][L_bins[i]].append(np.load('../results/erfc_err_{}{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                            file, ndim, dlogz, nlive)))
                    
                best_fit_param = np.load('../results/erfc_best_fit_params_wtheta{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(
                            file, ndim, dlogz, nlive))

                # read the results
                with open('../results/erfc_results{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(
                            file, ndim, dlogz, nlive), 'rb') as fp:
                    results = pickle.load(fp)
        
            except:
                
                if (np.round(N_q/1000)>1):
        
                    C = np.load('../results/cov_wtheta{}.npy'.format(
                        file))
                    if diag == True:
                        C*=np.identity(len(C))
            
                    log_det_C, solve = log_det(C, ix_grid)
                    
                    # initialize our nested sampler
                    # H0_c_dN_dz_array = np.array([H0_c_dN_dz_i.value for H0_c_dN_dz_i in H0_c_dN_dz_norm)
                    sampler = NestedSampler(loglike, ptform, ndim, logl_args=[log_det_C, wp_zbin_i[theta_mask], nbins_mask, solve, 
                                                                              H0_c_dN_dz_norm.value], nlive = nlive, rstate=np.random.default_rng(4*i+j))#np.random.default_rng(0)) #H0_c_dN_dz_norm[L_bins[i]][j]
                    
                    # run the sampler with checkpointing
                    sampler.run_nested(dlogz = dlogz)
                    results = sampler.results
                    
                    # Print out a summary of the results.
                    results.summary()
                    
                    # fig, axes = dyplot.runplot(results)  # summary (run) plot
                    # fig.subplots_adjust(hspace=0.25)
                    # fig.suptitle('$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(L_bins[i], z_bins[j], z_bins[j+1]), fontsize = 20, y = 0.9)
        
                    # fig, axes = dyplot.traceplot(results, 
                    #                              show_titles=True,
                    #                              trace_cmap='viridis', connect=True,
                    #                              connect_highlight=range(5), labels=labels, quantiles = quantile, title_quantiles = quantile)
                    # fig.subplots_adjust(hspace=0.7)
                    # fig.suptitle('$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(L_bins[i], z_bins[j], z_bins[j+1]), fontsize = 20, y = 1)
                
                    # # initialize figure
                    # fig, axes = plt.subplots(2, 2, figsize=(10, 10))
                    
                    # # plot initial run (res1; left)
                    # fg, ax = dyplot.cornerpoints(results, cmap='plasma',
                    #                              kde=False, labels=labels, fig=(fig, axes))
                    # fig.suptitle('$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(L_bins[i], z_bins[j], z_bins[j+1]), fontsize = 20, y = 0.975, x = 0.575)
        
                    # # # initialize figure
                    # fig, axes = plt.subplots(3, 3, figsize=(15, 15))
                    
                    # # plot initial run (res1; left)
                    # fg, ax = dyplot.cornerplot(results, color='blue', truths=np.zeros(ndim), hist2d_kwargs = {'contourf_kwargs': contourf_kwargs}, quantiles_2d = quantile_2d, title_quantiles = quantile,
                    #                            truth_color='black', show_titles=True,
                    #                            max_n_ticks=3, quantiles=quantile, labels=labels, fig=(fig, axes[:, :3]))
                    # fig.suptitle('$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(L_bins[i], z_bins[j], z_bins[j+1]), fontsize = 20, y = 1, x = 0.375)
        
                    samples, weights = results.samples, results.importance_weights()
                
                    # Get median and 1-sigma (16th, 50th, 84th percentiles) for each parameter
                    quantiles = [dyfunc.quantile(samples[:, i], quantile, weights=weights) 
                                 for i in range(samples.shape[1])]
                
                    print('{:.1f}<z<={:.1f}'.format(z_bins[j], z_bins[j+1]))
                    for k, (lo, mid, hi) in enumerate(quantiles):
                        print(f"{err_labels[k+1]}: {mid:.2f} + {hi-mid:.2f}/-{mid-lo:.2f}")
                        err[err_labels[k+1]][L_bins[i]].append((lo, hi))
                    print()
                    
                    equal_weight_samples = dyfunc.resample_equal(samples, weights)
                    
                    for k, err_label in enumerate(err_labels[1:]):
                        sample[err_label][L_bins[i]].append(equal_weight_samples[:,k]) #[ind]
                        median[err_label][L_bins[i]].append(quantiles[k][1])
                    
                    # Get the maximum likelihood sample from dynesty results
                    idx_max = np.argmax(results.logl)
                    best_fit_param = samples[idx_max]
                
                else:
        
                    equal_weight_samples = [np.nan]*ndim
        
                    for k, err_label in enumerate(err_labels[1:]):
                        sample[err_label][L_bins[i]].append([np.nan]) #[np.nan]*1000
                        median[err_label][L_bins[i]].append(np.nan)
                        err[err_label][L_bins[i]].append((np.nan, np.nan))
        
                    best_fit_param = [np.nan]*ndim
        
                    results = np.nan
        
                for k, err_label in enumerate(err_labels[1:]):
                    np.save('../results/erfc_samples_{}{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                            file, ndim, dlogz, nlive), sample[err_label][L_bins[i]][j])
                    np.save('../results/erfc_median_{}{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                            file, ndim, dlogz, nlive), median[err_label][L_bins[i]][j])
                    np.save('../results/erfc_err_{}{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                            file, ndim, dlogz, nlive), err[err_label][L_bins[i]][j])
        
                # # Get the maximum likelihood sample from dynesty results
                # best_fit_params[L_bins[i]].append(best_fit_param)
        
                np.save('../results/erfc_best_fit_params{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(
                            file, ndim, dlogz, nlive), best_fit_param)
                
                # results_Mmin[L_bins[i]].append(results)
                # save the results
                with open('../results/erfc_results{}_1h_ndim{}_dlogz{}_nlive{}.npy'.format(
                            file, ndim, dlogz, nlive), 'wb') as fp:
                    pickle.dump(results, fp)
                
    # f_duty, M_min, sigma = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # labels = ['log10(f_duty)', 'log10(M_min)', 'sigma']
    # err = {i: {L: [] for L in L_bins} for i in labels}
    # f_duty_samples, M_min_samples, sigma_samples = {L: [] for L in L_bins}, {L: [] for L in L_bins}, {L: [] for L in L_bins}
    # maxL_params = {L: [] for L in L_bins}
        
    # for i in range(len(L_bins)):

    #     print('M_i<={}'.format(L_bins[i]))
        
    #     for j in z_array:

    #         C = np.load('../results/cov_{}{}zmin{}zmax{}_Lmax{}_N{}_{}_MC_mask.npy'.format(args.corrfunc, file, z_bins[j], z_bins[j+1], L_bins[i], N_jk, cosmo))

    #         # check if cov is hermitian and positive definite
    #         if linalg.ishermitian(C) and np.all(np.linalg.eigvals(C) > 0):
    #             C_factor = linalg.cho_factor(C)
    #             log_det_C = 2 * np.sum(np.log(np.diag(C_factor[0])))
    #             solve = lambda rhs: linalg.cho_solve(C_factor, rhs)
    #         else:
    #             print('using LU decomposition')
    #             C_factor = linalg.lu_factor(C)
    #             log_det_C = np.sum(np.log(np.abs(np.diag(C_factor[0])))) # does abs value make sense
    #             solve = lambda rhs: linalg.lu_solve(C_factor, rhs)
        
    #         # initialize our nested sampler
    #         # if args.ndim == 2:
    #         sampler = NestedSampler(loglike, ptform, args.ndim, logl_args=[a_zbin[L_bins[i]][j], cf_matter_chi2[L_bins[i]][j], log_det_C, cf_chi2[L_bins[i]][j], n_q[L_bins[i]][j], nbins, solve, args.penalty, args.ndim], nlive = nlive, ptform_args = [args.ndim])
    #         # else:
    #         #     sampler = NestedSampler(loglike, ptform, ndim, logl_args=[a_zbin[L_bins[i]][j], cf_matter_chi2[L_bins[i]][j], log_det_C, cf_chi2[L_bins[i]][j], n_q[L_bins[i]][j], nbins, solve, args.penalty, args.HOD], nlive = nlive)
            
    #         # run the sampler with checkpointing
    #         start = time.time()
    #         sampler.run_nested(dlogz = dlogz, print_progress = False)
    #         results = sampler.results
    #         print(f"Iterations: {results.niter}")
    #         print(f"Time: {time.time() - start:.2f}s")
            
    #         # Print out a summary of the results.
    #         results.summary()        
    #         samples, weights = results.samples, results.importance_weights()
        
    #         # Get median and 1-sigma (16th, 50th, 84th percentiles) for each parameter
    #         quantiles = [dyfunc.quantile(samples[:, i], [0.16, 0.5, 0.84], weights=weights) 
    #                      for i in range(samples.shape[1])]
        
    #         print('\n{:.1f}<z<={:.1f}'.format(z_bins[j], z_bins[j+1]))
    #         for k, (lo, mid, hi) in enumerate(quantiles):
    #             print(f"{labels[k]}: {mid:.2f} + {hi-mid:.2f}/-{mid-lo:.2f}")
    #             err[labels[k]][L_bins[i]].append((lo, hi))
    #         print()
            
    #         equal_weight_samples = dyfunc.resample_equal(samples, weights)
    #         ind = np.random.choice(range(len(equal_weight_samples)), size = 1000, replace = False)
    #         f_duty_samples[L_bins[i]].append(equal_weight_samples[ind][:,0])
    #         M_min_samples[L_bins[i]].append(equal_weight_samples[ind][:,1])
            
    #         f_duty[L_bins[i]].append(quantiles[0][1])
    #         M_min[L_bins[i]].append(quantiles[1][1]) 

    #         x = [10**f_duty[L_bins[i]][j], M_min[L_bins[i]][j]]
    #         if args.ndim == 3:
    #             sigma_samples[L_bins[i]].append(equal_weight_samples[ind][:,2])
    #             sigma[L_bins[i]].append(quantiles[2][1]) 
    #             x.append(sigma[L_bins[i]][j])

    #         # Get the maximum likelihood sample from dynesty results
    #         idx_max = np.argmax(results.logl)
    #         maxL_params[L_bins[i]].append(samples[idx_max])

    #         # file = '_G{}_rmin{}_rmax{}_nbins{}_{}Lmax{}_N{}_{}_MC_mask'.format(G_hi, rmin, rmax, nbins, pi, L_bins[i], N_jk, cosmo)
                
    #         # np.save('../results/M_min_{}{}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), M_min[L_bins[i]])
    #         # np.save('../results/log10(M_min){}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), np.array(err['log10(M_min)'][L_bins[i]]))
    
    #         # np.save('../results/f_duty{}_{}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), f_duty[L_bins[i]])
    #         # np.save('../results/log10(f_duty){}_{}_sellentin_nlive{}_dlogz{}{}'.format(args.corrfunc, file, nlive, dlogz, args.penalty), np.array(err['log10(f_duty)'][L_bins[i]]))
    
    #         # np.save('../results/zbin_{}'.format(file), 1/np.array(a_zbin[L_bins[i]])-1)
        
    # # for i in range(len(L_bins)):
    
    # #     print('M_i<={}'.format(L_bins[i]))
    
    # #     for j in z_arr
            
    #         # print('{:.1f}<z<={:.1f}'.format(z_bins[j], z_bins[j+1]))
    #         # x = 10**f_duty[L_bins[i]][j], M_min[L_bins[i]][j]
            
    #         # C = np.load('../results/cov_{}_G{}_rmin{}_rmax{}_nbins{}_{}zmin{}zmax{}_Lmax{}_N{}_{}_MC_mask.npy'.format(args.corrfunc, G_hi, rmin, rmax, nbins, pi, z_bins[j], z_bins[j+1], L_bins[i], N_jk, cosmo))
    #         chi2 = f0_GLS_solve(x, a_zbin[L_bins[i]][j], cf_matter_chi2[L_bins[i]][j], C, cf_chi2[L_bins[i]][j], args.ndim)
    #         df = nbins-len(x)
        
    #         print('           chi^2: {:.2f}'.format(chi2))
    #         print('   reduced chi^2: {:.2f}'.format(chi2/df))
    #         print('         delta n: {:.2e}'.format(constraint(x, a_zbin[L_bins[i]][j], n_q[L_bins[i]][j], args.ndim)))
    #         print('condition number: {:.2f}'.format(np.linalg.cond(C)))
    #         print('             PTE: {:.2f}'.format(stats.chi2.sf(chi2, df)))

    #     np.save('../results/M_min_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), M_min[L_bins[i]])
    #     np.save('../results/log10(M_min)_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(err['log10(M_min)'][L_bins[i]]))

    #     np.save('../results/f_duty_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), f_duty[L_bins[i]])
    #     np.save('../results/log10(f_duty)_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(err['log10(f_duty)'][L_bins[i]]))

    #     if args.ndim == 3:
    #         np.save('../results/sigma_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), sigma[L_bins[i]])
    #         np.save('../results/sigma_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(err['sigma'][L_bins[i]]))

    #     np.save('../results/maxL_params_{}{}Lmax{}_N{}_{}_MC_mask_sellentin_nlive{}_dlogz{}{}_ndim{}'.format(args.corrfunc, file, L_bins[i], N_jk, cosmo, nlive, dlogz, args.penalty, args.ndim), np.array(maxL_params[L_bins[i]]))
    #     np.save('../results/zbin{}Lmax{}_{}_MC_mask'.format(file, L_bins[i], cosmo), 1/np.array(a_zbin[L_bins[i]])-1)

if __name__=='__main__':
    main()
    