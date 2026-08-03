import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import astropy.cosmology.units as cu
import astropy.units as u
import site
site.addsitedir('/oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts') 
import quaia
import seaborn as sns
import pyccl as ccl
from astropy.cosmology import Planck18
from dynesty import NestedSampler
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
from scipy import linalg
from scipy.special import gamma
import argparse
from scipy import special
import pickle
from astropy import constants
import mcfit
import matplotlib.colors as mcolors
from dynesty.pool import Pool
from matplotlib.backends.backend_pdf import PdfPages
from dynesty.utils import kld_error

# Model using step-function HOD
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

### Translate mass definition
dz = np.linspace(0, 4.6)

# Implement Limber approximation: $w(\theta)=\int d z \frac{H(z)}{c}\left(\frac{d N}{d z}\right)^2 \int \frac{d k k}{2 \pi} P_{q q} J_0[k \theta \chi(z)]=\int d z \frac{H(z)}{c}\left(\frac{d N}{d z}\right)^2 w_p(\theta \chi(z))$
k_log = np.logspace(-2.5, 3) #-4, 2)
da = 1/(1+dz)
P_mm = ccl.power.linear_matter_power(cosmo_planck18, k_log, da)
rpbins = (thetabins*u.deg).to(u.rad)[:, None]*quaia.comoving_dist(dz, units = 'Mpc')

## Pick cosmology
cosmo = 'Planck'
N_jk = 50 #40 #20 #65 #50 #30 #15 #80
NFW_200c = ccl.halos.profiles.nfw.HaloProfileNFW(mass_def = '200c', concentration = c_m)
d_hmf = [hmf(cosmo_planck18, bins, a) for a in da]
d_b1 = [b1_E_func_tinker10(cosmo_planck18, bins, a) for a in da]
d_NFW = [NFW_200c.fourier(cosmo_planck18, k_log, bins, a).T/bins for a in da]
hartlap = (N_jk-nbins-2)/(N_jk-1)

# create classes
class CholeskySolver:
    def __init__(self, C_factor):
        self.C_factor = C_factor
    def __call__(self, rhs):
        return linalg.cho_solve(self.C_factor, rhs)

class LUSolver:
    def __init__(self, C_factor):
        self.C_factor = C_factor
    def __call__(self, rhs):
        return linalg.lu_solve(self.C_factor, rhs)
        
## Define functions
def n_g_theory(bins, n_total, a = 1, cosmo = cosmo_planck18, hmf = hmf, recenter = True, units = '(h/Mpc)^3', d_hmf = None):

    if recenter == True:
        bins = quaia.recenter(bins)
    nm_value = d_hmf if d_hmf is not None else hmf(cosmo, bins, a)
    if units == '(h/Mpc)^3':
        nm=nm_value/cosmo['h']**3*(cu.littleh/u.Mpc)**3
    else:
        nm=nm_value*(1/u.Mpc)**3
    n_g = nm/(bins*np.log(10))*n_total
    n_g[np.isnan(n_g)]=0
    
    return n_g, nm # (units of Mpc/h)^3

def bias(bins, n_m, a = 1, cosmo = cosmo_planck18, mod = False, recenter = True, hmf = hmf, #method = method,
         b1_E_func_tinker10 = b1_E_func_tinker10, b1_E_tinker = None, n_g_manual = None, calc_bias = True, units = '(h/Mpc)^3', d_hmf = None, d_b1 = None):

    n_g, nm = n_g_theory(bins, n_m, a = a, cosmo = cosmo, recenter = recenter, hmf = hmf, units = units, d_hmf = d_hmf)
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
            
            # b1_E_tinker = b1_E_func_tinker10(cosmo, bins, a)
            b1_E_tinker = d_b1 if d_b1 is not None else b1_E_func_tinker10(cosmo, bins, a)
            
        n_m[~np.isfinite(n_m)]=0 
        integrand = dn_dm*n_m
        integrand[np.isnan(integrand)]=0

        return bias_function(n_g_manual, integrand, b1_E_tinker, bins, recenter = recenter), n_g

def bias_function(n_g, integrand, b, bins, recenter = True):
    
    if recenter == True:
        print('recentering')
        bins = quaia.recenter(bins)
        
    return 1/n_g*np.trapz(integrand*b, bins)

def HOD_scatter(bins, M_min, fduty, sigma):
    
    return fduty/2*special.erfc(np.log(10**M_min/bins)/(np.sqrt(2)*np.log(10)*sigma))

def N_c(bins_m, M_min, f_duty, ndim, M_max = M_max, sigma = 0):#, step = False): #, recenter = False): #bins_m

    n_m = HOD_scatter(bins_m, M_min, f_duty, sigma)

    if ndim == 2:

        n_m[bins_m==10**M_min]=f_duty # step

    else:

        n_m[bins_m==10**M_min]=f_duty/2 # scatter
        
    return n_m

def N_s(bins_m, fduty, M0, M1, alpha, exp = False, step = True): #bins_m

    if exp == False:
        n_m = fduty*((bins_m-M0)/M1)**alpha
        
        if step == True:
            n_m[bins_m<M0]=0
            
    else:
        n_m = N_s_exp(bins_m, M1, alpha, M0, fduty)
        
    return n_m

def N_s_exp(M, M1, alpha, Mcut, fduty):

    return fduty*(M/M1)**alpha*np.exp(-Mcut/M)

def HOD_2h(bins, M_min, fduty, sigma, a, M1, alpha, ndim, NFW_200c = NFW_200c, k = k_log, cosmo = cosmo_planck18, d_NFW = None): #, recenter = False): # method = method, 

    nfw = d_NFW if d_NFW is not None else NFW_200c.fourier(cosmo, k, bins, a).T/bins

    n_c = N_c(bins, M_min, 10**fduty, ndim, sigma = sigma)
    n_s = N_s(bins, 10**fduty, 10**M_min, 10**M1, alpha)
    
    return n_c+nfw*n_s, n_c, n_s #, recenter = recenter

def P_1h(bins, a, n_g, ndim, params_n_m = None, params_HOD = None, k = k_log, cosmo = cosmo_planck18, recenter = False, hmf = hmf, #n_c, n_s #method = method, 
         b1_E_func_tinker10 = b1_E_func_tinker10, NFW_200c = NFW_200c, damping = False, k_star = 1e-2/Planck18.h, units = '1/Mpc^3', d_hmf = None, d_b1 = None, d_NFW = None): #[mass_def]

    if params_n_m is not None:
        
        n_c, n_s = params_n_m

    else:
        
        fduty, M_min, sigma, M1, alpha = params_HOD
        n_c = N_c(bins, M_min, 10**fduty, ndim, sigma = sigma)
        n_s = N_s(bins, 10**fduty, 10**M_min, 10**M1, alpha)

    nfw = d_NFW if d_NFW is not None else NFW_200c.fourier(cosmo, k, bins, a).T/bins
    
    n_m = 2*nfw*n_c*n_s+nfw**2*n_s**2 #(n_s-1) # assume poisson distribution for satellites # should i implement central satellite condition? i think no
    P_qq_1h, _ = bias(bins, n_m, a, cosmo, recenter = recenter, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, #method = method,
                           b1_E_tinker = np.ones(np.shape(bins)), n_g_manual = n_g, units = units, d_hmf = d_hmf, d_b1 = d_b1)

    if damping == True:
        P_qq_1h *= (k/k_star)**4/(1+(k/k_star)**4)

    return P_qq_1h/n_g #**2
 
def wtheta_2h_dynesty(bins, H0_c_dN_dz, params_HOD, ndim, thetabins = thetabins, k = k_log, recenter = False, theta_mask = theta_mask, cosmo = cosmo_planck18, hankel = mcfit.Hankel(k_log, lowring = True), recenter_thetabins = True, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, #method = method
                      NFW_200c = NFW_200c, P_mm = P_mm, da = da, rpbins = rpbins, dz = dz, d_hmf = d_hmf, d_b1 = d_b1, d_NFW = d_NFW): #, k_integral = np.empty(np.shape(rpbins.T)), z_integrand = np.empty(np.shape(rpbins.T))):

    # compute for each z... for statement?
    b_q, P_qq = [], [] #k_integral #, []
    fduty, M_min, sigma, M1, alpha = params_HOD
    for i, a in enumerate(da):
        
        n_m, n_c, n_s = HOD_2h(bins, M_min, fduty, sigma, a, M1, alpha, ndim, cosmo = cosmo, k = k, NFW_200c = NFW_200c, d_NFW = d_NFW[i]) #, recenter = recenter)
    
        n_g = bias(bins, n_c+n_s, a, cosmo, recenter = recenter, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, 
                               b1_E_tinker = np.ones(np.shape(bins)), calc_bias = False, units = '1/Mpc^3', d_hmf = d_hmf[i], d_b1 = d_b1[i])
        
        b_z, _ = bias(bins, n_m, a, cosmo, recenter = recenter, hmf = hmf, b1_E_func_tinker10 = b1_E_func_tinker10, n_g_manual = n_g, units = '1/Mpc^3', d_hmf = d_hmf[i])
            
        P_qq_2h = b_z**2*P_mm[i]*u.Mpc**3
        P_qq_1h = P_1h(bins, a, n_g, ndim, params_HOD = params_HOD, k = k, cosmo = cosmo, recenter = recenter, hmf = hmf, 
                            b1_E_func_tinker10 = b1_E_func_tinker10, NFW_200c = NFW_200c, d_hmf = d_hmf[i], d_b1 = d_b1[i], d_NFW = d_NFW[i]) #n_c, n_s
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
        solve = CholeskySolver(C_factor) #lambda rhs: linalg.cho_solve(C_factor, rhs)
    else:
        print('using LU decomposition')
        C_factor = linalg.lu_factor(C[ix_grid])
        log_det_C = np.sum(np.log(np.abs(np.diag(C_factor[0])))) # does abs value make sense
        solve = LUSolver(C_factor)

    return log_det_C, solve
   
def loglike(x, log_det_C, wp, nbins, solve, dN_dz, ndim, theta_mask = theta_mask, M_max = M_max, N = N_jk, fduty = 0, recenter_thetabins = True, n_bins = n_bins, bins_m = bins):
        
    M_min, Alpha, sigma = [*x, 0][:3] # assign sigma = 0 if ndim==2
        
    wp_2h_zbin = wtheta_2h_dynesty(bins_m, dN_dz, [0, M_min, sigma, np.log10(12*10**M_min), Alpha], ndim)

    delta = wp-wp_2h_zbin
    
    c_p = gamma(N/2)/((np.pi*(N-1))**(nbins/2)*gamma((N-nbins)/2))
    lnorm = np.log(c_p)-0.5*log_det_C  

    return -N*0.5 * np.log(1+(delta.T @ solve(delta))/(N-1)) + lnorm #- penalty

def ptform(u, ndim):
    """Transforms the uniform random variables `u ~ Unif[0., 1.)`
    to the parameters of interest."""

    x = np.array(u)  # copy u 

    # Uniform from 9 to 16
    x[0] = 7 * u[0] + 9  # scale and shift to [9., 16.)

    # Uniform from 0 to 4
    x[1] = 4 * u[1]  # scale and shift to [0, 4.)

    if ndim == 3:
        
        # Uniform from 1 to 60
        x[2] = 1.4 * u[2]   # scale and shift to [0, 1.4.)

    return x

def make_contour_cmap(color, n_levels, alpha_min=0, alpha_max=1.0):
    
    cmap = []
    alphas = np.linspace(0,1,n_levels+1)
    for alpha in alphas:
        rgba = np.array(mcolors.to_rgba(color))
        rgba[-1] = alpha
        cmap.append(rgba)

    return cmap
    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("ndim", help="number of HOD parameters", type = int)
    parser.add_argument("zmin", help="minimum redshift", type = float)
    parser.add_argument("zmax", help="maximum redshift", type = float)
    parser.add_argument("Lmax", help="luminosity threshold", type = float)
    parser.add_argument('--plot', help="summary plots", action=argparse.BooleanOptionalAction, default = True)
    parser.add_argument('--nlive', help="number of live points", type = int, default = 500)
    args = parser.parse_args()

    matplotlib.rcParams['ytick.labelsize'] = 18
    matplotlib.rcParams['xtick.labelsize'] = 18
    matplotlib.rcParams['axes.labelsize'] = 22
    matplotlib.rcParams['legend.fontsize'] = 18
    matplotlib.rcParams['axes.titlesize'] = 20

    ## Define parameters
    G_hi = 20.5

    # Visualize the luminosity threshold cuts and redshift bins
    L_bins = np.array([-32.0, -27.0, -26.0, -25.0, -20.0])[:0:-1]
    b = 30
    fac_rand = 25

    ## Define the total occupation
    cmap_theta = sns.color_palette("flare_r", n_colors = len(L_bins))

    # Perform nested sampling with dynesty
    # Define the dimensionality of our problem.
    diag = False #True # False
    ix_grid = np.ix_(theta_mask, theta_mask)
    # nlive = 500
    dlogz = 0.01
    nbins_mask = np.sum(theta_mask)
    quantile = [0.16, 0.5, 0.84]
    labels=[r'$\log_{10} M_\mathrm{min}$', r'$\alpha$', r'$\Sigma$'][:args.ndim]
    quantile_2d = 1.0 - np.exp(-0.5 * np.array([1, 2])**2)#, 3])**2)
    contourf_kwargs = {'colors': make_contour_cmap(cmap_theta[list(L_bins).index(args.Lmax)], len(quantile))} #blue
    contour_kwargs = {'colors': [cmap_theta[list(L_bins).index(args.Lmax)]]}
    
    ## Across all $z$
    err_labels = ['log10(f_duty)', 'log10(M_min)', 'alpha', 'Sigma'][:args.ndim+1]
    thetamin, thetamax = 10**log_thetamin, 10**log_thetamax
          
    err = {}#i: {L: [] for L in L_bins} for i in err_labels}
    sample = {}#err_label: {L: [] for L in L_bins} for err_label in err_labels}
    median = {}#err_label: {L: [] for L in L_bins} for err_label in err_labels} 

    tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, selfunc_hi_bin0, _, _ = quaia.make_bins(G_hi, [0, 0], True, ['z', 'L'], 
                                                                                 bins = [[args.zmin, args.zmax], [args.Lmax]], tab_gcat_type = 'data', 
                                                                                  method = ['minmax', 'max'], fac_rand = fac_rand, b = b, 
                                                                                                          mask_type = 'b',
                                                                                                          n_bins = [None, None])

    N_q = len(tab_datahi_mask_zbin0)

    wp_zbin_i = np.load('../results/wtheta_G{}_rmin{}_rmax{}_nbins{}_zmin{}zmax{}_Lmax{}_{}_MC_b{}.npy'.format(G_hi, rmin, rmax, nbins, args.zmin, args.zmax, args.Lmax, #z_bins[args.j], 
                                                                                                               cosmo, b))
    
    dN_dz = np.sum(stats.norm.pdf(dz[None, :], loc=tab_datahi_mask_zbin0['redshift_quaia'][:, None], scale=tab_datahi_mask_zbin0['redshift_quaia_err'][:, None]), axis = 0)
    dN_dz /= np.sum(stats.norm.cdf(dz[-1], loc=tab_datahi_mask_zbin0['redshift_quaia'], scale=tab_datahi_mask_zbin0['redshift_quaia_err']))
    H0_c_dN_dz_norm = ccl.background.h_over_h0(cosmo_planck18, da)*Planck18.H0.to(1/u.s)/constants.c.to(u.Mpc/u.s)*dN_dz**2

    file = '_G{}_rmin{}_rmax{}_nbins{}_zmin{}zmax{}_Lmax{}_N{}_{}_MC_b{}'.format(G_hi, thetamin, thetamax, nbins_mask, args.zmin, args.zmax, args.Lmax, N_jk, cosmo, b) #z_bins[args.j], z_bins[args.j+1], 
    
    try:

        for k, err_label in enumerate(err_labels[1:]):
            
            sample[err_label] = np.load('../results/samples_{}{}_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                file, args.ndim, dlogz, args.nlive))
            median[err_label] = np.load('../results/median_{}{}_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                    file, args.ndim, dlogz, args.nlive))
            err[err_label] = np.load('../results/err_{}{}_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                    file, args.ndim, dlogz, args.nlive))
            
        best_fit_param = np.load('../results/best_fit_params_wtheta{}_ndim{}_dlogz{}_nlive{}.npy'.format(
                    file, args.ndim, dlogz, args.nlive))

        # read the results
        with open('../results/results{}_ndim{}_dlogz{}_nlive{}.npy'.format(
                    file, args.ndim, dlogz, args.nlive), 'rb') as fp:
            results = pickle.load(fp)

        print('sampler already completed in this bin')

    except:
        
        if (np.round(N_q/1000)>1):

            print('\nRunning new sampler...\n')
            
            C = np.load('../results/cov_wtheta_G{}_rmin{}_rmax{}_nbins{}_zmin{}zmax{}_Lmax{}_N{}_{}_MC_b{}.npy'.format(
            G_hi, rmin, rmax, nbins, args.zmin, args.zmax, args.Lmax, N_jk, cosmo, b))
            if diag == True:
                C*=np.identity(len(C))
    
            log_det_C, solve = log_det(C, ix_grid)
            
            # initialize our nested sampler
            with Pool(48, loglike, ptform,
                      logl_args=[log_det_C, wp_zbin_i[theta_mask], nbins_mask, solve,
                                 H0_c_dN_dz_norm.value, args.ndim],
                      ptform_args=[args.ndim]) as pool:
                sampler = NestedSampler(pool.loglike, pool.prior_transform, args.ndim, nlive = args.nlive, rstate=np.random.default_rng([int(args.zmin), int(args.zmax), np.abs(int(args.Lmax))]), pool = pool)

                # run the sampler with checkpointing
                sampler.run_nested(dlogz = dlogz)
                
            results = sampler.results
            
            # Print out a summary of the results.
            results.summary()
            print(results.logz[-1], results.logzerr[-1])

            kld = kld_error(results, error='jitter')
            print(np.mean(kld), np.std(kld))
            
            if args.plot:

                print('\nPlotting...\n')

                # https://matplotlib.org/stable/gallery/misc/multipage_pdf.html
                # Create the PdfPages object to which we will save the pages:
                # The with statement makes sure that the PdfPages object is closed properly at
                # the end of the block, even if an Exception occurs.
                with PdfPages('../figures/dyplot{}_ndim{}_dlogz{}_nlive{}.pdf'.format(
                        file, args.ndim, dlogz, args.nlive)) as pdf:

                    fig, axes = dyplot.traceplot(results, 
                                                 show_titles=True,
                                                 trace_cmap='viridis', connect=True,
                                                 connect_highlight=range(5), labels=labels, quantiles = quantile, title_quantiles = quantile)
                    fig.subplots_adjust(hspace=0.7)
                    fig.suptitle(r'$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(args.Lmax, args.zmin, args.zmax), fontsize = 20, y = 1)
                    pdf.savefig()  # saves the current figure into a pdf page
                    plt.close()

                    # initialize figure
                    fig, axes = plt.subplots(args.ndim-1, args.ndim-1, figsize=(args.ndim*5, args.ndim*5))
                    
                    # plot initial run (res1; left)
                    fg, ax = dyplot.cornerpoints(results, cmap='plasma',
                                                 kde=False, labels=labels, fig=(fig, axes))
                    # fig.suptitle(r'$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(args.Lmax, args.zmin, args.zmax), fontsize = 20, y = 0.975, x = 0.575)
                    pdf.savefig()  # saves the current figure into a pdf page
                    plt.close()
                        
                    # # initialize figure
                    fig, axes = plt.subplots(args.ndim, args.ndim, figsize=(args.ndim*5, args.ndim*5))
                    
                    # plot initial run (res1; left)
                    fg, ax = dyplot.cornerplot(results, color=cmap_theta[list(L_bins).index(args.Lmax)], hist2d_kwargs = {'contourf_kwargs': contourf_kwargs, 'contour_kwargs': contour_kwargs}, quantiles_2d = quantile_2d, title_quantiles = quantile, #'blue' #truths=np.zeros(args.ndim)
                                               show_titles=True, #truth_color='black'
                                               max_n_ticks=3, quantiles=quantile, labels=labels, fig=(fig, axes[:, :args.ndim]))
                    # fig.suptitle(r'$M_i\geq{}, {:.1f}<z \leq{:.1f}$'.format(args.Lmax, args.zmin, args.zmax), fontsize = 20, y = 1, x = 0.375)
                    pdf.savefig()  # saves the current figure into a pdf page
                    plt.close()
            
            samples, weights = results.samples, results.importance_weights()
        
            # Get median and 1-sigma (16th, 50th, 84th percentiles) for each parameter
            quantiles = [dyfunc.quantile(samples[:, i], quantile, weights=weights) 
                         for i in range(samples.shape[1])]
        
            for k, (lo, mid, hi) in enumerate(quantiles):
                print(f"{err_labels[k+1]}: {mid:.2f} + {hi-mid:.2f}/-{mid-lo:.2f}")
                np.save('../results/err_{}{}_ndim{}_dlogz{}_nlive{}.npy'.format(err_labels[k+1],
                        file, args.ndim, dlogz, args.nlive), (lo, hi))
            print()
            
            equal_weight_samples = dyfunc.resample_equal(samples, weights)
            
            for k, err_label in enumerate(err_labels[1:]):
                np.save('../results/samples_{}{}_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                        file, args.ndim, dlogz, args.nlive), equal_weight_samples[:,k]) #[ind]
                np.save('../results/median_{}{}_ndim{}_dlogz{}_nlive{}.npy'.format(err_label,
                        file, args.ndim, dlogz, args.nlive), quantiles[k][1])
            
            # Get the maximum likelihood sample from dynesty results
            idx_max = np.argmax(results.logl)
            np.save('../results/best_fit_params{}_ndim{}_dlogz{}_nlive{}.npy'.format(
                        file, args.ndim, dlogz, args.nlive), samples[idx_max])

            # save the results
            with open('../results/results{}_ndim{}_dlogz{}_nlive{}.npy'.format(
                        file, args.ndim, dlogz, args.nlive), 'wb') as fp:
                pickle.dump(results, fp)
    
            print('\nSaved to ../results/*{}_ndim{}_dlogz{}_nlive{}.npy!'.format(
                            file, args.ndim, dlogz, args.nlive))

        else:

            print('\nSkipping this bin, only {} quasars!'.format(N_q))
          
if __name__=='__main__':
    main()
    