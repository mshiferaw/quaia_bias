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
fac_stdev = 1.5 #1.45

# read quaia data
def read(fn_gcatlo, G_lo, fn_sello, NSIDE = NSIDE, b = 0, mask_type = None, plot = False, 
         name_catalog = '$Gaia$-$unWISE$ Quasar Catalog', fac_stdev = fac_stdev, cmap_map = 'plasma', cbar_ticks = [5, 10, 20], 
         cbar_ticks_selfunc = [5, 20, 50]):
    
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
    b : {int, float}, optional
        The cut-off for the mask in galactic latitude. Default is 0. Leave as 0 if unmasked catalog is desired. The choice does not matter for a selection function-based mask.
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
        
    Returns
    -------
    tab_gcatlo : astropy table
        The unmasked Quaia catalog.
    pixel_indices_gcatlo_mask : astropy table column
        Pixel indices for each object in the masked `tab_gcatlo` catalog.
    N_gcatlo_mask : int
        The number of objects in the masked `tab_gcatlo` catalog.
    mask_gcatlo : ndarray
        The mask for the catalog, of length `tab_gcatlo`.
        
    Other Parameters
    ----------------
    plot : bool, optional
        Whether or not to plot the healpy map of the catalog. Default is False.
    name_catalog : str, optional
        The title of the healpy plot. Default is '$Gaia$-$unWISE$ Quasar Catalog'.
    fac_stdev : {int, float}, optional
        The factor of standard deviation for the healpy plot. Default is 1.5.
    cmap_map : str, optional
        The colormap for the healpy plot. Default is 'plasma'.
    cmap_ticks : array_like, optional
        The ticks of the colorbar for the healpy plot. Default is [5, 10, 20].
    cmap_ticks_selfunc : array_like, optional
        The ticks of the colorbar for the selection function-corrected healpy plot. Default is [5, 20, 50].
    """
    
    NPIX = hp.nside2npix(NSIDE)
    tab_gcatlo = Table.read(fn_gcatlo)
    N_gcatlo = len(tab_gcatlo)
    
    print(f"Number of data sources: {N_gcatlo}")
    print(tab_gcatlo.meta)
    print(f"Column names: {tab_gcatlo.columns}")
            
    try:
        mask_gcatlo = np.abs(tab_gcatlo['b'])>=b
    except:
        c = SkyCoord(ra=tab_gcatlo['ra'].value*u.degree, dec=tab_gcatlo['dec'].value*u.degree)
        mask_gcatlo = np.abs(c.galactic.b.value)>=b

    # make map of quasar number counts
    pixel_indices_gcatlo = hp.ang2pix(NSIDE, tab_gcatlo['ra'][mask_gcatlo], tab_gcatlo['dec'][mask_gcatlo], lonlat=True)
    
    # to do selfunc based mask
    selfunc_lo = hp.fitsfunc.read_map(fn_sello)
    if mask_type == 'selfunc':
        pixel_indices_gcatlo = hp.ang2pix(NSIDE, tab_gcatlo['ra'], tab_gcatlo['dec'], lonlat=True)
        mask_gcatlo = selfunc_lo[pixel_indices_gcatlo]>=0.5
        pixel_indices_gcatlo = pixel_indices_gcatlo[mask_gcatlo]
    
    N_gcatlo_mask = len(tab_gcatlo[mask_gcatlo])
    
    if plot == True:
        
        map_gcatlo = np.bincount(pixel_indices_gcatlo, minlength=NPIX)
        title_gcatlo = rf"{name_catalog}, $G<{G_lo}$ (N={N_gcatlo_mask:,})"
        projview(map_gcatlo, title=title_gcatlo,
                    unit=r"number density per healpixel (deg$^{-2}$)", cmap=cmap_map, coord=['C', 'G'], 
                    min=np.median(map_gcatlo)-fac_stdev*np.std(map_gcatlo), max=np.median(map_gcatlo)+fac_stdev*np.std(map_gcatlo), 
                    norm='log', graticule=True, cbar_ticks=cbar_ticks)
        
        # remove the selection function
        map_selfunc_lo = map_gcatlo/selfunc_lo
        title_gcatlo = rf"{name_catalog}, $G<{G_lo}$ (N={N_gcatlo_mask:,})"
        projview(map_selfunc_lo, title=title_gcatlo,
                    unit=r"number density per healpixel (deg$^{-2}$)", cmap=cmap_map, coord=['C', 'G'], 
                    min=np.nanmedian(map_selfunc_lo)-fac_stdev*np.nanstd(map_selfunc_lo), max=np.nanmedian(map_selfunc_lo)+fac_stdev*np.nanstd(map_selfunc_lo), 
                    norm='log', graticule=True, cbar_ticks=cbar_ticks_selfunc) 
        
    return tab_gcatlo, pixel_indices_gcatlo, N_gcatlo_mask, mask_gcatlo

# convert z to comoving distance in Mpc/h
def comoving_dist(z, h = 0.6844, units = 'Mpc/h'): # col 2 in fig 7 of https://arxiv.org/pdf/1807.06209

    r"""Transforms redshift to comoving distance.
    
    Returns comoving distances in Mpc/h given a cosmology and redshift. Uses Astropy units and cosmology. Assumes flat Lambda-CDM, with Om0=0.302.
    
    Parameters
    ----------
    z : {int, float, array_like}
        The input redshift.
    h : {int, float}, optional
        Little h. Default is 0.6844.
    units : str, optional
        Units to return comoving distance in. Options are 'Mpc/h' and 'pc'. Default is 'Mpc/h'.
        
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
    if units == 'Mpc/h':
        return (comoving_r*cu.littleh).to(u.Mpc, cu.with_H0(H0))/cu.littleh # equivalent to comoving_d*h 
    else:
        return comoving_r.to(u.pc)

def absolute(tab_datalo, G, dust = np.load('../data/maps/map_dust_NSIDE64.npy'), h = 0.6844, NSIDE = NSIDE, 
             k = np.loadtxt('../data/maps/datafile4.txt')):
    
    r"""Transforms apparent to absolute magnitude.
    
    Returns absolute magnitude in the i-band given a cosmology and redshift. Uses Astropy units and cosmology. Assumes flat Lambda-CDM, with Om0=0.302.
    
    Parameters
    ----------
    tab_datalo : astropy table
        The input catalog.
    G : {int, float}
        The apparent G-band magnitude threshold cut of the input catalog.
        
    Returns
    -------
    M : float
        Absolute i-band magnitude at a given redshift `z` and cosmology.
        
    Other Parameters
    ----------------
    dust : array_like, optional
        The dust map. Default is the Quaia dust map.
    h : {int, float}, optional
        Little h. Default is 0.6844.
    NSIDE : int, optional
        The healpix nside parameter, must be a power of 2, less than 2**30. Default is 64.
    """
    
    G_i = {20: 0.09177204284667972, 20.5: 0.10518160949706967} # {20: 0.09178003723144457, 20.5: 0.10519113616943265}
    
    m_g = tab_datalo['phot_g_mean_mag'] # need to convert to i band!
    m_i = m_g - G_i[G]  # 0.05688247070312613  # 0.08305740356445312
    
    H0 = h*100 * u.km/u.s/u.Mpc

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.
    cosmo = FlatLambdaCDM(H0=H0, Om0=0.302)
    d = cosmo.luminosity_distance(tab_datalo['redshift_quaia']).value
    
    pixel_indices_datalo = hp.ang2pix(NSIDE, tab_datalo['ra'], tab_datalo['dec'], lonlat=True)
    A_i = 1.698*dust[pixel_indices_datalo]
    
    k_z = interp.interp1d(k[:,0], k[:,1])
    K_i = k_z(tab_datalo['redshift_quaia'])
    
    tab_datalo['M_i']=m_i-25-5*np.log10(d)-A_i-K_i # d is in Mpc

def recenter(bins):
    
    r"""Computes bin centers.
    
    Returns the centers of a given bin array. Will have length len(`bins`)-1. Useful for plotting the results of `quaia.w_theta` given the input `bins`.
    
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
    quaia.w_theta
    """
    
    return 0.5*(bins[1:]+bins[:-1])

def z_dist(tab_gcatlo, tab_randlo, mask_gcatlo, mask_randlo, N_gcatlo_mask, N_randlo_mask, b = 0, mask_type = None, plot = False, new = False):
    
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
        The number of objects in the masked `tab_gcatlo` catalog.
    N_randlo_mask : int
        The number of objects in the masked `tab_randlo` catalog.
    b : {int, float}, optional
        The cut-off for the mask in galactic latitude. Default is 0. Leave as 0 if unmasked catalog is desired.
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
        
    Returns
    -------
    key : str
        The suffix of the new key created for the random redshifts. Use as ``tab_randlo['redshift_quaia_'+key].``
        
    Other Parameters
    ----------------
    plot : boolean, optional
        Whether or not to display a histogram of the redshift distributions for `tab_gcatlo` (as given) and `tab_randlo` (as assigned). Default is False.
    new : boolean, optional
        Whether or not to return the modified random catalogue `tab_randlo` alongside the suffix of the key for the random redshifts.
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
    if type(mask_type) == str:
        key = mask_type
    else:
        key = str(b)
    tab_randlo['redshift_quaia_'+key] = -1.0
    tab_randlo['redshift_quaia_'+key][mask_randlo] = inv_ecdf(rng.random(size = N_randlo_mask))
    
    if new == True:
        return key, tab_randlo
    else:
        return key

def make_bins(G, bb, prebinned, cut, method = None, tab_gcat = None, tab_gcat_type = None, i = None, n_bins = None, bins = None, z = True, fac_rand = 10, allsky = "", NSIDE = NSIDE, plot = False, mask_type = 'selfunc'):
    
    r"""Computes catalogs binned in redshift or luminosity with a selection-function based mask.
    
    Based on the G threshold cut and the method of binning (Lsplit vs LminLmax), returns the catalog in a specific redshift bin.
    
    Parameters
    ----------
    G : {int, float}
        Threshold cut of input catalog.
    bb : int
        Desired redshift bin.
    prebinned : boolean
        Whether or not the input catalog is prebinned, or if the bins must be computed from the whole input catalog.
    cut : str
        Which dimension to bin along. Options are "z" or "L".
    method : str, optional
        Which binning method to use. Options are "split" and "minmax". Set to None if `prebinned` is True.
    tab_gcat_type : str, optional
        The type of input catalog. Options are "data" and "mock". Default is None. Leave as None if providing a prebinned input catalog via `tab_gcat`.
    n_bins : int, optional
        The number of bins for the "split" bin method. Default is None. Leave as None if using "minmax" method.
    bins : array_like, optional
        The bins using the "minmax" method. Set `n_bins` = None to avoid the "split" method.
    z : boolean, optional
        Whether or not to assign redshifts to the random catalog. Default is True. Set to False if only intending to compute angular clustering via `quaia.w_theta`. Uses `quaia.z_dist`. 
    fac_rand : int, optional
        The factor by which the random catalog exceeds the input catalog. Default is 10. 
        
    Returns
    -------
    tab_datahi_mask_bin0 : array_like
        The binned quasar catalog.
    tab_randhi_mask_bin0 : array_like
        The binned random catalog
    key_bin0 : str
        The suffix of the new key created for the random redshifts. Use as ``tab_randhi_mask_bin0['redshift_quaia_'+key_bin0].``
    bins : array_like
        The final bins used to compute `tab_datahi_mask_bin0` and `tab_randhi_mask_bin0`.
    weights1 : array_like
        Weights for the first set of points. Pass into `quaia.w_theta`, `quaia.wp_rp`, and `quaia.xi_s` to upweight data by the selection function.
        
    Other Parameters
    ----------------
    tab_gcat : array_like, optional
        The input catalog. Default is None, as it will be loaded based on `G` and `n_bins`. Optional to pass in if another input catalog is desired.
    i : int, optional
        The index of the mock to use as the input catalog. Options are any integer from 1-100. Default is None. Set `prebinned` = False and `tab_gcat_type` = 'mock' if using `i`.
    allsky : str, optional
        Suffux for random catalog filename. Default is "". Set to "_allsky" if using allsky random catalog. 
    NSIDE : int, optional
        The healpix nside parameter, must be a power of 2, less than 2**30. Default is 64.
    plot : boolean, optional
        Whether or not to display a histogram of the redshift distributions for the data and random catalogs. Default is False.
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.

    See Also
    --------
    quaia.z_dist
    quaia.w_theta
    quaia.wp_rp
    quaia.xi_s
    """
        
    # select binning method
    if method == 'split':
        fname = '_G{:.1f}_{}split{}bin{}'.format(G, cut, n_bins, bb)
    elif method == 'minmax':
        fname = '_G{:.1f}_{}min{}{}max{}'.format(G, cut, bins[bb], cut, bins[bb+1])
    else:
        print('Leaving tab_gcat_catalog unbinned')
        fname = '_G{:.1f}'.format(G)
        if prebinned == False:
            raise Exception('To leave tab_gcat_catalog unbinned, set prebinned = True.')
    # print(fname)
    
    # for prebinned catalogs (only data)
    if prebinned == True:
        if tab_gcat_type != 'data':
            raise Exception('Only data catalogs are prebinned currently')
        fn_datahi_bin0 = '../data/quaia{}.fits'.format(fname)
        tab_datahi_bin0 = Table.read(fn_datahi_bin0)
        
    # for catalogs that need to be binned
    elif prebinned == False:
        
        if tab_gcat_type == 'mocks':
            if i is None:
                raise Exception('Select which mock to run.')
            if G==20.0:
                G=20 # for mock file naming convention
            fn_datahi = '../Quaia_mock_catalogs-20250211T213139Z-001/Quaia_mock_catalogs/G{}/with_selection_function/mock_catalog_quaia_G{}_mock{}.fits'.format(G, G, i)
            tab_gcat = Table.read(fn_datahi)
        elif tab_gcat_type == 'data':
            fn_datahi = '../data/quaia_G{:.1f}.fits'.format(G)
            tab_gcat = Table.read(fn_datahi)
        else:
            if tab_gcat is None:
                raise Exception('Must provide unbinned tab_gcat catalog if binned == False and tab_gcat_type != mocks or data')
            
        if cut == 'z':
            key = 'redshift_quaia'
        else:
            absolute(tab_gcat, G)
            key = 'M_i'
                    
        if method == 'split':
                z_percentiles = np.linspace(0.0, 100.0, n_bins+1)
                # print(z_percentiles)
                bins = np.percentile(list(tab_gcat[key]), z_percentiles)
                bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
                bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

                # print("bins:", bins)
                # print("n_bins:", n_bins)

        i_bin = (tab_gcat[key] >= bins[bb]) & (tab_gcat[key] < bins[bb+1])
        tab_datahi_bin0 = tab_gcat[i_bin]
        # print("min:", np.min(tab_datahi_bin0[key]))
        # print("max:", np.max(tab_datahi_bin0[key]))
    else:
        raise Exception('prebinned bust be a boolean')
            
    # load selection function
    fn_selhi_bin0 = '../data/maps/selection_function_NSIDE64{}.fits'.format(fname) 
    selfunc_hi_bin0 = hp.fitsfunc.read_map(fn_selhi_bin0)

    # quasar data catalog
    pixel_indices_datahi_bin0 = hp.ang2pix(NSIDE, tab_datahi_bin0['ra'], tab_datahi_bin0['dec'], lonlat=True)
    
    # impose bin and selfunc mask on data
    mask_datahi_bin0 = selfunc_hi_bin0[pixel_indices_datahi_bin0]>=0.5 
    tab_datahi_mask_bin0 = tab_datahi_bin0[mask_datahi_bin0]
    
    # record N in each bin
    N_datahi_mask_bin0 = len(tab_datahi_mask_bin0)
    print('{}: {:.2f}              {}'.format(fname[1:], 1-N_datahi_mask_bin0/len(tab_datahi_bin0), N_datahi_mask_bin0))
    
    # random catalog
    fn_randhi_bin0 = '../data/randoms/random{}{}_{}x.fits'.format(fname, allsky, fac_rand)
    tab_randhi_bin0 = Table.read(fn_randhi_bin0)
        
    # get the pixel indices for objects in each random bin
    pixel_indices_randhi_bin0 = hp.ang2pix(NSIDE, tab_randhi_bin0['ra'], tab_randhi_bin0['dec'], lonlat=True)
    
    # select where the randoms in each bin pass the selfunc mask
    mask_randhi_bin0 = selfunc_hi_bin0[pixel_indices_randhi_bin0]>=0.5
    tab_randhi_mask_bin0 = tab_randhi_bin0[mask_randhi_bin0]

    # assign redshifts to randoms in each bin, mimicking the distribution of data in each bin that pass the selfunc mask
    N_randhi_mask_bin0 = len(tab_randhi_mask_bin0)
    if z == True:
        key_bin0 = z_dist(tab_datahi_mask_bin0, tab_randhi_mask_bin0, np.full(N_datahi_mask_bin0, True), 
                             np.full(N_randhi_mask_bin0, True), N_datahi_mask_bin0, N_randhi_mask_bin0, plot = plot, 
                             mask_type = mask_type+'_bin0')
    else:
        key_bin0 = None
    
    return tab_datahi_mask_bin0, tab_randhi_mask_bin0, key_bin0, bins, 1/selfunc_hi_bin0[pixel_indices_datahi_bin0][mask_datahi_bin0], 1-N_datahi_mask_bin0/len(tab_datahi_bin0)

# 2D angular clustering w(theta)
def w_theta(tab_gcatlo, tab_randlo, selfunc_lo = None, weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, 
            thetabins = np.logspace(np.log10(0.1), np.log10(10.0), 15)):
    
    r"""Computes the angular clustering of `tab_gcatlo`.
    
    Returns the two-point angular correlation function computed as a function of `thetabins` centers. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    nthreads: int, optional
        The number of OpenMP threads to use for computing each set of pair counts. Default is 8.
    thetabins : array_like, optional
        The angular bins. Default is ``np.logspace(np.log10(0.1), np.log10(10.0), 15)``.
        
    Returns
    -------
    wp : array_like
        The angular correlation function.
        
    Other Parameters
    ----------------
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    weights1 : None, optional
        Weights for the first set of points. Default is None if not weighting pair counts by the selection function. 
    weights2 : None, optional
        Weights for the second set of points. Default is None if not weighting pair counts by the selection function.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default is None if pre-computed pair counts are not supplied.
        
    See Also
    --------
    quaia.recenter : Used to compute the bin centers of `thetabins`.
    quaia.read: Used to obtain the pixel indices and mask for `tab_gcatlo` and `tab_randlo`.
    """
    
    if type(weights1) == np.ndarray or type(weights2) == np.ndarray:
        print('weighting DR')
        weight_type12 = 'pair_product'
    else:
        weight_type12 = None
    if type(weights1) == np.ndarray:
        print('weighting DD')
        weight_type1 = 'pair_product'
    else:
        weight_type1 = None
        weights1 = np.ones(len(tab_gcatlo))
    if type(weights2) == np.ndarray:
        print('weighting RR')
        weight_type2 = 'pair_product'
    else:
        weight_type2 = None
        weights2 = np.ones(len(tab_randlo))
        
    # comoving distance
    DD_counts, api_time = mocks.DDtheta_mocks(autocorr = 1, nthreads = nthreads, binfile = thetabins, RA1 = tab_gcatlo['ra'], 
                                          DEC1 = tab_gcatlo['dec'], weights1 = weights1,
                                          weight_type=weight_type1, c_api_timer = True)
    print('DD: {}'.format(api_time))
    
    # now measure clustering in random catalog
    DR_counts, api_time = mocks.DDtheta_mocks(autocorr = 0,nthreads = nthreads, binfile = thetabins, RA1 = tab_gcatlo['ra'], 
                                          DEC1 = tab_gcatlo['dec'], weights1 = weights1, 
                                          RA2 = tab_randlo['ra'], DEC2 = tab_randlo['dec'], 
                                          weights2 = weights2, weight_type=weight_type12, 
                                          c_api_timer = True)
    print('DR: {}'.format(api_time))
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDtheta_mocks(autocorr = 1, nthreads = nthreads, binfile = thetabins, RA1 = tab_randlo['ra'], 
                                         DEC1 = tab_randlo['dec'], weights1 = weights2, 
                                         weight_type=weight_type2, c_api_timer = True)
        print('RR: {}'.format(api_time))
        RR_counts['thetaavg'] = np.mean([RR_counts['thetamin'], RR_counts['thetamax']], axis = 0)
        
        if weight_type2 == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
    
    # compute bin centers for theta
    DD_counts['thetaavg'] = np.mean([DD_counts['thetamin'], DD_counts['thetamax']], axis = 0)
    DR_counts['thetaavg'] = np.mean([DR_counts['thetamin'], DR_counts['thetamax']], axis = 0)

    # compute weighted pair counts
    if weight_type1 == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    if weight_type12 == 'pair_product':
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    N_gcatlo = np.sum(weights1)
    N_randlo = np.sum(weights2)
    return convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)

# 3D projected clustering wp(rp)
def wp_rp(tab_gcatlo, tab_randlo, key, selfunc_lo = None, weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, 
          rbins = np.logspace(np.log10(0.5), np.log10(60.0), 21), nbins = 20, pimax = 40.0, d = 1):
    
    r"""Computes the two-point projected 3D clustering of `tab_gcatlo`.
    
    Returns the projected correlation function computed as a function of projected quasar separation. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    key : str
        The suffix of the key for redshifts assigned to `tab_randlo`, created under the desired masking scheme using quaia.z_dist.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default option is None if pre-computed pair counts are not supplied.
    nthreads: int, optional
        The number of OpenMP threads to use for computing each set of pair counts. Default is 8.
    rbins : array_like, optional
        The projected separation bins. Default is ``np.logspace(np.log10(0.5), np.log10(60.0), 21)``.
    nbins : int, optional
        The number of bins in `rbins`. Default is 20.
    pimax : {int, float}, optional
        The maximum separation along the line-of-sight, in the z-direction. Default is 40.0.
        
    Returns
    -------
    wp : array_like
        The projected correlation function.
        
    Other Parameters
    ----------------
    weights1 : None, optional
        Weights for the first set of points. Default is None if not weighting pair counts by the selection function. 
    weights2 : None, optional
        Weights for the second set of points. Default is None if not weighting pair counts by the selection function.
    d : int, optional
        The factor by which to downsample `tab_randlo` in order to reduce computation time. Default is 1. Increasing `d` comes at the cost of noise in the final clustering measurement.
        
    See Also
    --------
    quaia.z_dist : Used to assign redshifts and a mask to `tab_randlo`.
    quaia.read: Used to obtain the pixel indices and mask for `tab_gcatlo` and `tab_randlo`.
    """
    
    if type(weights1) == np.ndarray or type(weights2) == np.ndarray:
        print('weighting DR')
        weight_type12 = 'pair_product'
    else:
        weight_type12 = None
    if type(weights1) == np.ndarray:
        print('weighting DD')
        weight_type1 = 'pair_product'
    else:
        weight_type1 = None
        weights1 = np.ones(len(tab_gcatlo))
    if type(weights2) == np.ndarray:
        print('weighting RR')
        weight_type2 = 'pair_product'
    else:
        weight_type2 = None
        weights2 = np.ones(len(tab_randlo))
        
    # comoving distance 
    DD_counts, api_time = mocks.DDrppi_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'],  # where hubble distance = c/H0 and H0 = 100 km/s/Mpc h
                                         CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1[::d],
                                         is_comoving_dist = True, weight_type=weight_type1, output_rpavg = True, c_api_timer = True) 
    
    DR_counts, api_time = mocks.DDrppi_mocks(autocorr = 0, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                         CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1[::d],
                                         RA2 = tab_randlo['ra'][::d], DEC2 = tab_randlo['dec'][::d],
                                         CZ2 = comoving_dist(tab_randlo['redshift_quaia_'+key][::d]), 
                                         weights2 = weights2[::d], weight_type=weight_type12,
                                         is_comoving_dist = True, output_rpavg = True, c_api_timer = True)

    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDrppi_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_randlo['ra'][::d], DEC1 = tab_randlo['dec'][::d], 
                                         CZ1 = comoving_dist(tab_randlo['redshift_quaia_'+key][::d]), 
                                         weights1 = weights2[::d], weight_type=weight_type2,
                                         is_comoving_dist = True, output_rpavg = True, c_api_timer = True)

        if weight_type2 == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']

    # compute weighted pair counts
    if weight_type1 == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    if weight_type12 == 'pair_product':
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    N_gcatlo = np.sum(weights1[::d])
    N_randlo = np.sum(weights2[::d])
    wp = convert_rp_pi_counts_to_wp(N_gcatlo, N_gcatlo, int(N_randlo), int(N_randlo),
                                DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
    
    # calculate the x axis
    rpavg = [np.sum((DD_counts['rpavg']*DD_counts['npairs'])[DD_counts['rmin']==i])/np.sum(DD_counts['npairs'][DD_counts['rmin']==i]) 
         for i in rbins[:-1]]
    
    return wp, rpavg

# 3d clustering xi(r) computed with 1 mu bin
def xi_s(tab_gcatlo, tab_randlo, key, selfunc_lo = None, weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, 
         rbins = np.logspace(np.log10(0.1), np.log10(20.0), 21)):
    
    r"""Computes the two-point  3D clustering of `tab_gcatlo`.
    
    Returns the correlation function computed as a function of quasar separation in redshift-space. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    key : str
        The suffix of the key for redshifts assigned to `tab_randlo`, created under the desired masking scheme using quaia.z_dist.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default option is None if pre-computed pair counts are not supplied.
    nthreads: int, optional
        The number of OpenMP threads to use for computing each set of pair counts. Default is 8.
    rbins : array_like, optional
        The separation bins. Default is ``np.logspace(np.log10(0.1), np.log10(20.0), 21)``.
        
    Returns
    -------
    cf : array_like
        The 3D correlation function.
        
    Other Parameters
    ----------------
    weights1 : None, optional
        Weights for the first set of points. Default is None if not weighting pair counts by the selection function. 
    weights2 : None, optional
        Weights for the second set of points. Default is None if not weighting pair counts by the selection function.
        
    See Also
    --------
    quaia.z_dist : Used to assign redshifts and a mask to `tab_randlo`.
    quaia.read: Used to obtain the pixel indices and mask for `tab_gcatlo` and `tab_randlo`.
    """        
        
    if type(weights1) == np.ndarray or type(weights2) == np.ndarray:
        print('weighting DR')
        weight_type12 = 'pair_product'
    else:
        weight_type12 = None
    if type(weights1) == np.ndarray:
        print('weighting DD')
        weight_type1 = 'pair_product'
    else:
        weight_type1 = None
        weights1 = np.ones(len(tab_gcatlo))
    if type(weights2) == np.ndarray:
        print('weighting RR')
        weight_type2 = 'pair_product'
    else:
        weight_type2 = None
        weights2 = np.ones(len(tab_randlo))
    
    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.   
    DD_counts, api_time = mocks.DDsmu_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, binfile = rbins, mu_max = 1, nmu_bins = 1, 
                                            RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                            CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1, 
                                            weight_type = weight_type1, output_savg = True, is_comoving_dist=True, c_api_timer = True)

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.   
    DR_counts, api_time = mocks.DDsmu_mocks(autocorr = 0, cosmology = 2, nthreads = nthreads, binfile = rbins, mu_max = 1, nmu_bins = 1, 
                                            RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                            CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1, 
                                            RA2 = tab_randlo['ra'], DEC2 = tab_randlo['dec'], 
                                            CZ2 = comoving_dist(tab_randlo['redshift_quaia_'+key]), weights2 = weights2, 
                                            weight_type = weight_type12, output_savg = True, is_comoving_dist = True, c_api_timer = True)
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDsmu_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, binfile = rbins, mu_max = 1, nmu_bins = 1,
                                                RA1 = tab_randlo['ra'], DEC1 = tab_randlo['dec'], 
                                                CZ1 = comoving_dist(tab_randlo['redshift_quaia_'+key]), weights1 = weights2, 
                                                weight_type = weight_type2, output_savg = True, is_comoving_dist=True, c_api_timer = True)
        if weight_type2 == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
        
    # compute weighted pair counts
    if weight_type1 == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    if weight_type12 == 'pair_product':
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    N_gcatlo = np.sum(weights1)
    N_randlo = np.sum(weights2)
    cf = convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)
    
    # make it easier to plot
    ind = np.argsort(DD_counts['savg'])
    
    return cf
    
# 3d clustering xi(r)
def xi_r(tab_gcatlo, tab_randlo, key, selfunc_lo = None, weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, 
         rbins = np.logspace(np.log10(0.1), np.log10(20.0), 21), correction = -5205.182232081812):
    
    r"""Computes the two-point  3D clustering of `tab_gcatlo`.
    
    Returns the correlation function computed as a function of quasar separation. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    key : str
        The suffix of the key for redshifts assigned to `tab_randlo`, created under the desired masking scheme using quaia.z_dist.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default option is None if pre-computed pair counts are not supplied.
    nthreads: int, optional
        The number of OpenMP threads to use for computing each set of pair counts. Default is 8.
    rbins : array_like, optional
        The separation bins. Default is ``np.logspace(np.log10(0.1), np.log10(20.0), 21)``.
    correction: {int, float}, optional
        The distance, in Mpc/h, by which to shift the catalog to ensure that it exists in a box of coordinates [0, boxsize].
        
    Returns
    -------
    cf : array_like
        The 3D correlation function.
        
    Other Parameters
    ----------------
    weights1 : None, optional
        Weights for the first set of points. Default is None if not weighting pair counts by the selection function. 
    weights2 : None, optional
        Weights for the second set of points. Default is None if not weighting pair counts by the selection function.
        
    See Also
    --------
    quaia.z_dist : Used to assign redshifts and a mask to `tab_randlo`.
    quaia.read: Used to obtain the pixel indices and mask for `tab_gcatlo` and `tab_randlo`.
    """
    
    if type(weights1) == np.ndarray or type(weights2) == np.ndarray:
        print('weighting DR')
        weight_type12 = 'pair_product'
    else:
        weight_type12 = None
    if type(weights1) == np.ndarray:
        print('weighting DD')
        weight_type1 = 'pair_product'
    else:
        weight_type1 = None
        weights1 = np.ones(len(tab_gcatlo))
    if type(weights2) == np.ndarray:
        print('weighting RR')
        weight_type2 = 'pair_product'
    else:
        weight_type2 = None
        weights2 = np.ones(len(tab_randlo))

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.     
    c = SkyCoord(ra=tab_gcatlo['ra'].value*u.degree, dec=tab_gcatlo['dec'].value*u.degree, distance=comoving_dist(tab_gcatlo['redshift_quaia']))
    X1, Y1, Z1 = c.cartesian.xyz.value - correction 
    DD_counts, api_time = theory.DD(autocorr = 1, nthreads = nthreads, binfile = rbins, periodic = False,
                                    X1 = X1, Y1 = Y1, Z1 = Z1, weights1 = weights1,
                                    weight_type=weight_type1, output_ravg = True, c_api_timer = True) # cz/H0 = Mpc/h = m/s * km/1000 m / (100 km/s/Mpc h)

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.     
    c = SkyCoord(ra=tab_randlo['ra'], dec=tab_randlo['dec'], distance=comoving_dist(tab_randlo['redshift_quaia_'+key]))
    X2, Y2, Z2 = c.cartesian.xyz.value - correction
    DR_counts, api_time = theory.DD(autocorr = 0, nthreads = nthreads, binfile = rbins, periodic = False,
                                    X1 = X1, Y1 = Y1, Z1 = Z1, weights1 = weights1, 
                                    X2 = X2, Y2 = Y2, Z2 = Z2, weights2 = weights2,
                                    weight_type=weight_type12, output_ravg = True, c_api_timer = True)
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = theory.DD(autocorr = 1, nthreads = nthreads, binfile = rbins, periodic = False,
                                X1 = X2, Y1 = Y2, Z1 = Z2, weights1 = weights2, 
                                weight_type=weight_type2, output_ravg = True, c_api_timer = True)
        if weight_type2 == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
        
    # compute weighted pair counts
    if weight_type1 == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
    if weight_type12 == 'pair_product':
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    N_gcatlo = np.sum(weights1)
    N_randlo = np.sum(weights2)
    
    # All the pair counts are done, get the angular correlation function
    cf = convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)
    
    # make it easier to plot
    ind = np.argsort(DD_counts['ravg'])
    
    return cf[ind], DD_counts['ravg'][ind]