import numpy as np
from Corrfunc import mocks
from Corrfunc.utils import convert_3d_counts_to_cf
import healpy as hp
from astropy.table import Table
from healpy.newvisufunc import projview
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import astropy.cosmology.units as cu
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc import theory
from astropy.coordinates import SkyCoord
import scipy.interpolate as interp
from matplotlib import pyplot as plt

# global variables
NSIDE = 64
G_lo = 20.0
fac_stdev = 1.5

# read quaia data
def read(fn_gcatlo, G_lo, fn_sello, NSIDE = NSIDE, mask = 20, plot = False, name_catalog = '$Gaia$-$unWISE$ Quasar Catalog', fac_stdev = 1.45, cmap_map = 'plasma'):
    
    r"""Read in a Quaia catalog (e.g., data, randoms, mocks) given a mask.
    
    Returns the full unmasked catalog, as well as the pixel indices for each object in the masked catalog, the number of objects in the masked catalog, and the mask itself. Transforms RA and DEC coordinates to galactic coordinates for non-data catalogs.
    
    Parameters
    ----------
    fn_gcatlo : str
        The location of the catalog in the user's directory.
    G_lo : {int, float}
        The magnitude cut of the catalog.
    fn_sello : str
        The location of the selection function catalog in the user's directory. 
    NSIDE : int, optional
        The healpix nside parameter, must be a power of 2, less than 2**30. Default is 64.
    mask : {int, float}
        The cut-off for the mask in galactic latitude b. Default is 20. Pass in 0 if unmasked catalog is desired.
        
    Returns
    -------
    tab_gcatlo : astropy table
        The unmasked Quaia catalog.
    pixel_indices_gcatlo_mask : astropy table column
        Pixel indices for each object in the masked catalog.
    N_gcatlo_mask : int
        The number of objects in the masked catalog.
    mask_gcatlo : ndarray
        The mask for the catalog, of length `tab_gcatlo`.
        
    Other Parameters
    ----------------
    plot : bool, optional
        Whether or not to plot the healpy map of the catalog. Default is False.
    name_catalog : str, optional
        The title of the healpy plot. Default is '$Gaia$-$unWISE$ Quasar Catalog'.
    fac_stdev : {int, float}, optional
        The factor of standard deviation for the healpy plot. Default is 1.45.
    cmap_map : str, optional
        The colormap for the healpy plot. Default is 'plasma'.
    """
    
    NPIX = hp.nside2npix(NSIDE)
    tab_gcatlo = Table.read(fn_gcatlo)
    N_gcatlo = len(tab_gcatlo)
    
    print(f"Number of data sources: {N_gcatlo}")
    print(tab_gcatlo.meta)
    print(f"Column names: {tab_gcatlo.columns}")
        
    try:
        mask_gcatlo = np.abs(tab_gcatlo['b'])>=mask
    except:
        c = SkyCoord(ra=tab_gcatlo['ra'].value*u.degree, dec=tab_gcatlo['dec'].value*u.degree)
        mask_gcatlo = np.abs(c.galactic.b.value)>=mask
    
    # make map of quasar number counts
    pixel_indices_gcatlo = hp.ang2pix(NSIDE, tab_gcatlo['ra'][mask_gcatlo], tab_gcatlo['dec'][mask_gcatlo], lonlat=True)
    
    if plot == True:
        
        map_gcatlo = np.bincount(pixel_indices_gcatlo, minlength=NPIX)
        title_gcatlo = rf"{name_catalog}, $G<{G_lo}$ (N={len(tab_gcatlo):,})"
        projview(map_gcatlo, title=title_gcatlo,
                    unit=r"number density per healpixel (deg$^{-2}$)", cmap=cmap_map, coord=['C', 'G'], 
                    min=np.median(map_gcatlo)-fac_stdev*np.std(map_gcatlo), max=np.median(map_gcatlo)+fac_stdev*np.std(map_gcatlo), 
                    norm='log', graticule=True,
                    cbar_ticks=[5, 10, 20]) 
        
        # remove the selection function
        selfunc_lo = hp.fitsfunc.read_map(fn_sello)
        map_selfunc_lo = map_gcatlo/selfunc_lo
        title_gcatlo = rf"{name_catalog}, $G<{G_lo}$ (N={len(tab_gcatlo):,})"
        projview(map_selfunc_lo, title=title_gcatlo,
                    unit=r"number density per healpixel (deg$^{-2}$)", cmap=cmap_map, coord=['C', 'G'], 
                    min=np.nanmedian(map_selfunc_lo)-fac_stdev*np.nanstd(map_selfunc_lo), max=np.nanmedian(map_selfunc_lo)+fac_stdev*np.nanstd(map_selfunc_lo), 
                    norm='log', graticule=True,
                    cbar_ticks=[5, 20, 50]) 
        
    return tab_gcatlo, pixel_indices_gcatlo, len(tab_gcatlo[mask_gcatlo]), mask_gcatlo

# convert z to comoving distance in Mpc/h
def comoving_dist(z, h = 0.6844): # col 2 in fig 7 of https://arxiv.org/pdf/1807.06209

    r"""Transforms redshift to comoving distance.
    
    Returns comoving distances in Mpc/h given a cosmology and redshift. Uses Astropy units and cosmology. Assumes flat Lambda-CDM, with Om0=0.302.
    
    Parameters
    ----------
    z : {int, float, array_like}
        The input redshift.
    h : {int, float}, optional
        Little h.
        
    Returns
    -------
    distance : {int, float, array_like}
        Comoving distance at a given redshift `z` and cosmology.
    """
    
    H0 = h*100 * u.km/u.s/u.Mpc

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.
    cosmo = FlatLambdaCDM(H0=H0, Om0=0.302)
    comoving_r = cosmo.comoving_distance(z)

    # convert from Mpc to Mpc/h
    return (comoving_r*cu.littleh).to(u.Mpc, cu.with_H0(H0))/cu.littleh # equivalent to comoving_d*h 

def recenter(bins):
    
    r"""Computes bin centers.
    
    Returns the centers of a given bin array. Will have length len(`bins`)-1. Useful for plotting the results of `quaia.wp_rp` and `quaia.xi_r` given the input `bins`.
    
    Parameters
    ----------
    bins : array_like
        The input bins.
        
    Returns
    -------
    bins_center : array_like
        The bin centers.
        
    See Also
    --------
    quaia.wp_rp, quaia.xi_r
    """
    
    return 0.5*(bins[1:]+bins[:-1])

def z_dist(tab_gcatlo, tab_randlo, mask_gcatlo, mask_randlo, N_gcatlo_mask, N_randlo_mask, mask = 20, plot = False):
    
    r"""Assigns redshifts to `tab_randlo` given the redshift distribution of `tab_gcatlo`.
    
    Creates a new column in `tab_randlo` for each object's assigned redshift, given the redshift distribution of the masked `tab_gcatlo` catalog. Leaves the masked objects of `tab_randlo` with a redshift of -1.0. A random seed is fixed prior to interpolation of the cumulative redshift distribution.
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The input catalog with the desired redshift distribution.
    tab_randlo : astropy table
        The target catalog that is assigned redshifts.
    mask_gcatlo : array_like
        The mask on `tab_gcatlo`.
    mask_randlo : array_like
        The mask on `tab_randlo`.
    N_gcatlo_mask : int
        The number of objects in `tab_gcatlo`.
    N_randlo_mask : int
        The number of objects in `tab_randlo`.
    mask : {int, float}, optional 
        The cut-off for the mask in galactic latitude b. Default is 20. Pass in 0 if unmasked catalog is desired.
    plot : boolean, optional
        Whether or not to display a histogram of the redshift distributions for `tab_gcatlo` (as given) and `tab_randlo` (as assigned). Default is False
    """
    
    # Compute the empirical cumulative distribution function (ECDF)
    z = np.sort(tab_gcatlo['redshift_quaia'][mask_gcatlo])
    ecdf = np.arange(N_gcatlo_mask) / (N_gcatlo_mask - 1)

    # Create an interpolation function for inverse transform sampling
    inv_ecdf = interp.interp1d(ecdf, z) # have to interpolate, can't do binning method because each point --> cdf probability

    # set random seed: https://builtin.com/data-science/numpy-random-seed
    rng = np.random.default_rng(2023)
    
    if plot == True:
        
        plt.plot(z, ecdf, label = 'ecdf')
        plt.hist(tab_gcatlo['redshift_quaia'][mask_gcatlo], cumulative = True, density = True, bins = 100)
        plt.xlabel('$z$')
        plt.legend()
        plt.show()

    # Generate new random samples from the estimated distribution
    tab_randlo['redshift_quaia_'+str(mask)] = -1.0
    tab_randlo['redshift_quaia_'+str(mask)][mask_randlo] = inv_ecdf(rng.random(size = N_randlo_mask))

# 2D angular clustering w(theta)
def w_theta(tab_gcatlo, tab_randlo, selfunc_lo, pixel_indices_gcatlo, pixel_indices_randlo, N_gcatlo, N_randlo, RR_counts,
            thetabins = np.logspace(np.log10(0.1), np.log10(10.0), 15), nthreads = 8):
    
    # comoving distance
    DD_counts, api_time = mocks.DDtheta_mocks(autocorr = 1, nthreads = nthreads, binfile = thetabins, RA1 = tab_gcatlo['ra'], 
                                          DEC1 = tab_gcatlo['dec'], weights1 = 1/selfunc_lo[pixel_indices_gcatlo],
                                          weight_type='pair_product', c_api_timer = True)
    
    # now measure clustering in random catalog
    DR_counts, api_time = mocks.DDtheta_mocks(autocorr = 0,nthreads = nthreads, binfile = thetabins, RA1 = tab_gcatlo['ra'], 
                                          DEC1 = tab_gcatlo['dec'], weights1 = 1/selfunc_lo[pixel_indices_gcatlo], 
                                          RA2 = tab_randlo['ra'], DEC2 = tab_randlo['dec'], 
                                          weights2 = 1/selfunc_lo[pixel_indices_randlo], weight_type='pair_product', 
                                          c_api_timer = True)
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDtheta_mocks(autocorr = 1, nthreads = nthreads, binfile = thetabins, RA1 = tab_randlo['ra'], 
                                         DEC1 = tab_randlo['dec'], weights1 = 1/selfunc_lo[pixel_indices_randlo], 
                                         weight_type='pair_product', c_api_timer = True)
        RR_counts['thetaavg'] = np.mean([RR_counts['thetamin'], RR_counts['thetamax']], axis = 0)
        RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
    
    # compute bin centers for theta
    DD_counts['thetaavg'] = np.mean([DD_counts['thetamin'], DD_counts['thetamax']], axis = 0)
    DR_counts['thetaavg'] = np.mean([DR_counts['thetamin'], DR_counts['thetamax']], axis = 0)

    # compute weighted pair counts
    DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    return convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)

# 3D projected clustering wp(rp)
def wp_rp(tab_gcatlo, tab_randlo, selfunc_lo, pixel_indices_gcatlo, pixel_indices_randlo, N_gcatlo, N_randlo, RR_counts, rmin = 0.5, rmax = 60.0, 
               nbins = 20, pimax = 40.0, d = 1, mask = 20, nthreads = 8):
    
    # create the bins array
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
    
    # comoving distance 
    DD_counts, api_time = mocks.DDrppi_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'],  # where hubble distance = c/H0 and H0 = 100 km/s/Mpc h
                                         CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = 1/selfunc_lo[pixel_indices_gcatlo],
                                         is_comoving_dist = True, weight_type='pair_product', output_rpavg = True, c_api_timer = True) 
    
    DR_counts, api_time = mocks.DDrppi_mocks(autocorr = 0, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                         CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = 1/selfunc_lo[pixel_indices_gcatlo],
                                         RA2 = tab_randlo['ra'][::d], DEC2 = tab_randlo['dec'][::d],
                                         CZ2 = comoving_dist(tab_randlo['redshift_quaia_'+str(mask)][::d]), 
                                         weights2 = 1/selfunc_lo[pixel_indices_randlo][::d], weight_type='pair_product', 
                                         is_comoving_dist = True, output_rpavg = True, c_api_timer = True)
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDrppi_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_randlo['ra'][::d], DEC1 = tab_randlo['dec'][::d], 
                                         CZ1 = comoving_dist(tab_randlo['redshift_quaia_'+str(mask)][::d]), 
                                         weights1 = 1/selfunc_lo[pixel_indices_randlo][::d], weight_type='pair_product',
                                         is_comoving_dist = True, output_rpavg = True, c_api_timer = True)
        RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
    
    # compute weighted pair counts
    DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    wp = convert_rp_pi_counts_to_wp(N_gcatlo, N_gcatlo, int(N_randlo/d), int(N_randlo/d),
                                DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
    
    # calculate the x axis
    rpavg = [np.sum((DD_counts['rpavg']*DD_counts['npairs'])[DD_counts['rmin']==i])/np.sum(DD_counts['npairs'][DD_counts['rmin']==i]) 
         for i in rbins[:-1]]
    
    return wp, rpavg

# 3d clustering xi(r)
def xi_r(tab_gcatlo, tab_randlo, selfunc_lo, pixel_indices_gcatlo, pixel_indices_randlo, N_gcatlo, N_randlo, RR_counts, correction,
         rbins = np.logspace(np.log10(0.1), np.log10(20.0), 21), mask = 20, nthreads = 8):
    
    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.     
    c = SkyCoord(ra=tab_gcatlo['ra'].value*u.degree, dec=tab_gcatlo['dec'].value*u.degree, distance=comoving_dist(tab_gcatlo['redshift_quaia']))
    X1, Y1, Z1 = c.cartesian.xyz.value - correction 
    DD_counts, api_time = theory.DD(autocorr = 1, nthreads = nthreads, binfile = rbins, periodic = False,
                                    X1 = X1, Y1 = Y1, Z1 = Z1, weights1 = 1/selfunc_lo[pixel_indices_gcatlo],
                                    weight_type='pair_product', output_ravg = True, c_api_timer = True) # cz/H0 = Mpc/h = m/s * km/1000 m / (100 km/s/Mpc h)
    
    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.     
    c = SkyCoord(ra=tab_randlo['ra'], dec=tab_randlo['dec'], distance=comoving_dist(tab_randlo['redshift_quaia_'+str(mask)]))
    X2, Y2, Z2 = c.cartesian.xyz.value - correction
    DR_counts, api_time = theory.DD(autocorr = 0, nthreads = nthreads, binfile = rbins, periodic = False,
                                    X1 = X1, Y1 = Y1, Z1 = Z1, weights1 = 1/selfunc_lo[pixel_indices_gcatlo], 
                                    X2 = X2, Y2 = Y2, Z2 = Z2, weights2 = 1/selfunc_lo[pixel_indices_randlo],
                                    weight_type='pair_product', output_ravg = True, c_api_timer = True)
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = theory.DD(autocorr = 1, nthreads = nthreads, binfile = rbins, periodic = False,
                                X1 = X2, Y1 = Y2, Z1 = Z2, weights1 = 1/selfunc_lo[pixel_indices_randlo], 
                                weight_type='pair_product',output_ravg = True, c_api_timer = True)
        RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
        
    # compute weighted pair counts
    DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    cf = convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)
    
    # make it easier to plot
    ind = np.argsort(DD_counts['ravg'])
    
    return cf[ind], DD_counts['ravg'][ind]