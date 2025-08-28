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
    selfunc_lo : array_like
        The selection function. Default is ``None``. Leave as ``None`` if using a galactic latitude-based mask.
        
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
                    # min=0.1, max=np.median(map_gcatlo)+fac_stdev*np.std(map_gcatlo), 
                    min=np.median(map_gcatlo)-fac_stdev*np.std(map_gcatlo), max=np.median(map_gcatlo)+fac_stdev*np.std(map_gcatlo), 
                    norm='log', graticule=True, cbar_ticks=cbar_ticks)
                    # cbar_ticks=[5, 10, 20]) 
        
        # remove the selection function
        map_selfunc_lo = map_gcatlo/selfunc_lo
        title_gcatlo = rf"{name_catalog}, $G<{G_lo}$ (N={N_gcatlo_mask:,})"
        projview(map_selfunc_lo, title=title_gcatlo,
                    unit=r"number density per healpixel (deg$^{-2}$)", cmap=cmap_map, coord=['C', 'G'], 
                    # min=0.1, max=np.nanmedian(map_selfunc_lo)+fac_stdev*np.nanstd(map_selfunc_lo), 
                    min=np.nanmedian(map_selfunc_lo)-fac_stdev*np.nanstd(map_selfunc_lo), max=np.nanmedian(map_selfunc_lo)+fac_stdev*np.nanstd(map_selfunc_lo), 
                    norm='log', graticule=True, cbar_ticks=cbar_ticks_selfunc) 
                    # cbar_ticks=[5, 20, 50]) 
        
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

# # convert z to comoving distance in Mpc/h
# def luminosity_dist(z, h = 0.6844, units = 'Mpc/h'): # col 2 in fig 7 of https://arxiv.org/pdf/1807.06209

#     r"""Transforms redshift to luminosity distance.
    
#     Returns comoving distances in Mpc/h given a cosmology and redshift. Uses Astropy units and cosmology. Assumes flat Lambda-CDM, with Om0=0.302.
    
#     Parameters
#     ----------
#     z : {int, float, array_like}
#         The input redshift.
#     h : {int, float}, optional
#         Little h.
        
#     Returns
#     -------
#     distance : {int, float, array_like}
#         Luminosity distance at a given redshift `z` and cosmology.
#     """
    
#     H0 = h*100 * u.km/u.s/u.Mpc

#     # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.
#     cosmo = FlatLambdaCDM(H0=H0, Om0=0.302)
#     comoving_r = cosmo.luminosity_distance(z)

#     # # convert from Mpc to Mpc/h
#     # if units == 'Mpc/h':
#     #     return (comoving_r*cu.littleh).to(u.Mpc, cu.with_H0(H0))/cu.littleh # equivalent to comoving_d*h 
#     # else:
#     #     return comoving_r.to(u.pc)
#     return comoving_r

def absolute(tab_datalo, dust = np.load('../data/maps/map_dust_NSIDE64.npy'), h = 0.6844, NSIDE = NSIDE, 
             k = np.loadtxt('../data/maps/datafile4.txt')):
    
    r"""Transforms apparent to absolute magnitude.
    
    Returns absolute magnitude in the G-band given a cosmology and redshift. Uses Astropy units and cosmology. Assumes flat Lambda-CDM, with Om0=0.302.
    
    Parameters
    ----------
    tab_datalo : astropy table
        The input catalog.
    dust : array_like, optional
        The dust map. Default is the Quaia dust map.
    h : {int, float}, optional
        Little h. Default is 0.6844.
    NSIDE : int, optional
        The healpix nside parameter, must be a power of 2, less than 2**30. Default is 64.
        
    Returns
    -------
    M : float
        Absolute magnitude at a given redshift `z` and cosmology.
    """
    
    m_i = tab_datalo['phot_g_mean_mag'] # need to convert to i band!
    
    H0 = h*100 * u.km/u.s/u.Mpc

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.
    cosmo = FlatLambdaCDM(H0=H0, Om0=0.302)
    d = cosmo.luminosity_distance(tab_datalo['redshift_quaia']).value
    
    pixel_indices_datalo = hp.ang2pix(NSIDE, tab_datalo['ra'], tab_datalo['dec'], lonlat=True)
    A_i = 1.698*dust[pixel_indices_datalo]
    
    k_z = interp.interp1d(k[:,0], k[:,1])
    K_i = k_z(tab_datalo['redshift_quaia'])
    
    tab_datalo['G']=m_i-25-5*np.log10(d)-A_i-K_i # d is in Mpc

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

# def make_zbins(n_zbins, tab_gcat):
    
#     z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
#     print(z_percentiles)
#     z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
#     z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
#     z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

#     print("zbins:", z_bins)
#     print("n_zbins:", n_zbins)
    
#     return z_bins
    
def zbins(G, bb, i=1, plot = False, n_zbins = 2, NSIDE = NSIDE, data = 'data', mask_type = 'selfunc', tab_gcat = None, tab_datahi_zbin0 = None, z_bins = [0.0, 1.0, 2.0, 3.0, 4.0]):
    
    r"""Computes catalogs binned in redshift with a selection-function based mask.
    
    Based on the G threshold cut and the method of binning (zsplit vs zminzmax), returns the catalog in a specific redshift bin.
    
    Parameters
    ----------
    G : {int, float}
        Threshold cut of input catalog.
    bb : int
        Desired redshift bin.
    i : float, optional
        Mock number. Default is 1.
    n_zbins : int, optional
        The number of redshift bins. Default is 2. Set to None if using zminzmax method.
    data : str, optional
        The type of input catalog. Options are "data" or "mocks".
    tab_gcat : array_like, optional
        The input catalog. Default is None, as it will be loaded based on `G` and `n_zbins`. Optional to pass in if another input catalog is desired.
    tab_datahi_zbin0 : array_like, optional
        The input catalog in the desired redshift bin. Default is None, as it will be calculated based on `G`, `bb`, and `n_zbins`. Optional to pass in a pre-binned input catalog.
    z_bins : array_like, optional
        The redshift bins using the zminzmax method. Default is [0.0, 1.0, 2.0, 3.0, 4.0]. Set n_zbins = None to avoid the zsplit method. 
        
    Returns
    -------
    key : str
        The suffix of the new key created for the random redshifts. Use as ``tab_randlo['redshift_quaia_'+key].``
        
    Other Parameters
    ----------------
    NSIDE : int, optional
        The healpix nside parameter, must be a power of 2, less than 2**30. Default is 64.
    plot : boolean, optional
        Whether or not to display a histogram of the redshift distributions for `tab_datahi_mask_zbin0` and `tab_randhi_mask_zbin0`. Default is False.
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
    """
    
    if data == 'mocks':
        
        # check if unbinned catalog is provided
        if tab_gcat is None:
            fn_datahi = '../Quaia_mock_catalogs-20250211T213139Z-001/Quaia_mock_catalogs/G{}/with_selection_function/mock_catalog_quaia_G{}_mock{}.fits'.format(G, G, i)
            tab_gcat = Table.read(fn_datahi)

        if n_zbins is not None:
            z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
            print(z_percentiles)
            z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
            z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
            z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

            print("zbins:", z_bins)
            print("n_zbins:", n_zbins)

            # z_bins = make_zbins(n_zbins, tab_gcat)

        # for bb in range(n_zbins):
        i_zbin = (tab_gcat['redshift_quaia'] >= z_bins[bb]) & (tab_gcat['redshift_quaia'] < z_bins[bb+1])
        tab_datahi_zbin0 = tab_gcat[i_zbin]
        print("zmin:", np.min(tab_datahi_zbin0['redshift_quaia']))
        print("zmax:", np.max(tab_datahi_zbin0['redshift_quaia']))
    
    elif data == 'data':
        
        if n_zbins == 1:
            fn_datahi_zbin0 = '../data/quaia_G{:.1f}.fits'.format(G)            
            fn_selhi_zbin0 = '../data/maps/selection_function_NSIDE64_G{:.1f}.fits'.format(G)
            fn_randhi_zbin0 = '../data/randoms/random_G{:.1f}_10x.fits'.format(G)
        
        else:
            fn_datahi_zbin0 = '../data/quaia_G{:.1f}_zsplit{}bin{}.fits'.format(G, n_zbins, bb)
            
        tab_datahi_zbin0 = Table.read(fn_datahi_zbin0)
        
    if n_zbins is None:
        
        # parameters
        fn_selhi_zbin0 = '../data/maps/selection_function_NSIDE64_G{:.1f}_zmin{}zmax{}.fits'.format(G, z_bins[bb], z_bins[bb+1])
        fn_randhi_zbin0 = '../data/randoms/random_G{:.1f}_zmin{}zmax{}_10x.fits'.format(G, z_bins[bb], z_bins[bb+1])
        
    elif n_zbins>1:
        # parameters
        fn_selhi_zbin0 = '../data/maps/selection_function_NSIDE64_G{:.1f}_zsplit{}bin{}.fits'.format(G, n_zbins, bb)
        fn_randhi_zbin0 = '../data/randoms/random_G{:.1f}_zsplit{}bin{}_10x.fits'.format(G, n_zbins, bb)
        
    selfunc_hi_zbin0 = hp.fitsfunc.read_map(fn_selhi_zbin0)

    # quasar data catalog
    pixel_indices_datahi_zbin0 = hp.ang2pix(NSIDE, tab_datahi_zbin0['ra'], tab_datahi_zbin0['dec'], lonlat=True)
    
    mask_datahi_zbin0 = selfunc_hi_zbin0[pixel_indices_datahi_zbin0]>=0.5 

    # impose zbin and selfunc mask on data
    tab_datahi_mask_zbin0 = tab_datahi_zbin0[mask_datahi_zbin0]
    
    # record N in each bin
    N_datahi_mask_zbin0 = len(tab_datahi_mask_zbin0)
    
    # random catalog
    tab_randhi_zbin0 = Table.read(fn_randhi_zbin0)
        
    # get the pixel indices for objects in each random zbin
    pixel_indices_randhi_zbin0 = hp.ang2pix(NSIDE, tab_randhi_zbin0['ra'], tab_randhi_zbin0['dec'], lonlat=True)

    # select where the randoms in each zbin pass the selfunc mask
    mask_randhi_zbin0 = selfunc_hi_zbin0[pixel_indices_randhi_zbin0]>=0.5

    # assign redshifts to randoms in each zbin, mimicking the distribution of data in each zbin that pass the selfunc mask
    tab_randhi_mask_zbin0 = tab_randhi_zbin0[mask_randhi_zbin0]
    N_randhi_mask_zbin0 = len(tab_randhi_mask_zbin0)
    key_zbin0 = z_dist(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, np.full(N_datahi_mask_zbin0, True), 
                             np.full(N_randhi_mask_zbin0, True), N_datahi_mask_zbin0, N_randhi_mask_zbin0, plot = plot, 
                             mask_type = mask_type+'_zbin0')
    
    return tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0

def make_zbins(G, bb, prebinned, method = None, i=None, plot = False, n_zbins = None, NSIDE = NSIDE, tab_gcat_type = None, mask_type = 'selfunc', tab_gcat = None, z_bins = None):
    
    r"""Computes catalogs binned in redshift with a selection-function based mask.
    
    Based on the G threshold cut and the method of binning (zsplit vs zminzmax), returns the catalog in a specific redshift bin.
    
    Parameters
    ----------
    G : {int, float}
        Threshold cut of input catalog.
    bb : int
        Desired redshift bin.
    n_zbins : int, optional
        The number of redshift bins. Default is 2. Set to None if using zminzmax method.
    data : str, optional
        The type of input catalog. Options are "data" or "mocks".
    tab_gcat : array_like, optional
        The input catalog. Default is None, as it will be loaded based on `G` and `n_zbins`. Optional to pass in if another input catalog is desired.
    tab_datahi_zbin0 : array_like, optional
        The input catalog in the desired redshift bin. Default is None, as it will be calculated based on `G`, `bb`, and `n_zbins`. Optional to pass in a pre-binned input catalog.
    z_bins : array_like, optional
        The redshift bins using the zminzmax method. Default is [0.0, 1.0, 2.0, 3.0, 4.0]. Set n_zbins = None to avoid the zsplit method. 
        
    Returns
    -------
    key : str
        The suffix of the new key created for the random redshifts. Use as ``tab_randlo['redshift_quaia_'+key].``
        
    Other Parameters
    ----------------
    NSIDE : int, optional
        The healpix nside parameter, must be a power of 2, less than 2**30. Default is 64.
    plot : boolean, optional
        Whether or not to display a histogram of the redshift distributions for `tab_datahi_mask_zbin0` and `tab_randhi_mask_zbin0`. Default is False.
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
    """
        
    # select binning method
    if method == 'zsplit':
        fname = '_G{:.1f}_zsplit{}bin{}'.format(G, n_zbins, bb)
    elif method == 'zminzmax':
        fname = '_G{:.1f}_zmin{}zmax{}'.format(G, z_bins[bb], z_bins[bb+1])
    else:
        print('Leaving tab_gcat_catalog unbinned')
        fname = '_G{:.1f}'.format(G)
        if prebinned == True:
            raise Exception('To do so, set prebinned = False')
        
    # for prebinned catalogs (only data)
    if prebinned == True:
        if tab_gcat_type != 'data':
            raise Exception('Only data catalogs are prebinned currently')
        fn_datahi_zbin0 = '../data/quaia{}.fits'.format(fname)
        tab_datahi_zbin0 = Table.read(fn_datahi_zbin0)
        
    # for catalogs that need to be binned
    elif prebinned == False:
        
        if tab_gcat_type == 'mocks':
            fn_datahi = '../Quaia_mock_catalogs-20250211T213139Z-001/Quaia_mock_catalogs/G{}/with_selection_function/mock_catalog_quaia_G{}_mock{}.fits'.format(G, G, i)
            tab_gcat = Table.read(fn_datahi)
        elif tab_gcat_type == 'data':
            fn_datahi = '../data/quaia_G{}.fits'.format(G)
            tab_gcat = Table.read(fn_datahi)
        else:
            if tab_gcat is None:
                raise Exception('Must provide unbinned tab_gcat catalog if binned == False and tab_gcat_type != mocks or data')
            
        if method == 'zsplit':
                z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
                print(z_percentiles)
                z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
                z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
                z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

                print("zbins:", z_bins)
                print("n_zbins:", n_zbins)

        # for bb in range(n_zbins):
        i_zbin = (tab_gcat['redshift_quaia'] >= z_bins[bb]) & (tab_gcat['redshift_quaia'] < z_bins[bb+1])
        tab_datahi_zbin0 = tab_gcat[i_zbin]
        print("zmin:", np.min(tab_datahi_zbin0['redshift_quaia']))
        print("zmax:", np.max(tab_datahi_zbin0['redshift_quaia']))
            
    # load selection function
    fn_selhi_zbin0 = '../data/maps/selection_function_NSIDE64{}.fits'.format(fname) 
    selfunc_hi_zbin0 = hp.fitsfunc.read_map(fn_selhi_zbin0)

    # quasar data catalog
    pixel_indices_datahi_zbin0 = hp.ang2pix(NSIDE, tab_datahi_zbin0['ra'], tab_datahi_zbin0['dec'], lonlat=True)
    
    mask_datahi_zbin0 = selfunc_hi_zbin0[pixel_indices_datahi_zbin0]>=0.5 

    # impose zbin and selfunc mask on data
    tab_datahi_mask_zbin0 = tab_datahi_zbin0[mask_datahi_zbin0]
    
    # record N in each bin
    N_datahi_mask_zbin0 = len(tab_datahi_mask_zbin0)
    
    # random catalog
    fn_randhi_zbin0 = '../data/randoms/random{}_10x.fits'.format(fname)
    tab_randhi_zbin0 = Table.read(fn_randhi_zbin0)
        
    # get the pixel indices for objects in each random zbin
    pixel_indices_randhi_zbin0 = hp.ang2pix(NSIDE, tab_randhi_zbin0['ra'], tab_randhi_zbin0['dec'], lonlat=True)

    # select where the randoms in each zbin pass the selfunc mask
    mask_randhi_zbin0 = selfunc_hi_zbin0[pixel_indices_randhi_zbin0]>=0.5

    # assign redshifts to randoms in each zbin, mimicking the distribution of data in each zbin that pass the selfunc mask
    tab_randhi_mask_zbin0 = tab_randhi_zbin0[mask_randhi_zbin0]
    N_randhi_mask_zbin0 = len(tab_randhi_mask_zbin0)
    key_zbin0 = z_dist(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, np.full(N_datahi_mask_zbin0, True), 
                             np.full(N_randhi_mask_zbin0, True), N_datahi_mask_zbin0, N_randhi_mask_zbin0, plot = plot, 
                             mask_type = mask_type+'_zbin0')
    
    return tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0

# 2D angular clustering w(theta)
def w_theta(tab_gcatlo, tab_randlo, N_gcatlo, N_randlo, selfunc_lo = None, pixel_indices_gcatlo = None, pixel_indices_randlo = None, weight_type = None,
            weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, thetabins = np.logspace(np.log10(0.1), np.log10(10.0), 15)):
    
    r"""Computes the angular clustering of `tab_gcatlo`.
    
    Returns the two-point angular correlation function computed as a function of `thetabins` centers. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    N_gcatlo : int
        The number of objects in `tab_gcatlo`.
    N_randlo : int
        The number of objects in `tab_randlo`.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    pixel_indices_gcatlo : astropy table column, optional
        Pixel indices for each object in `tab_gcatlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    pixel_indices_randlo : astropy table column, optional 
        Pixel indices for each object in `tab_randlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    weight_type : {None, 'pair_product'}, optional
        Whether or not to weight pair counts by the selection function. Default is None.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default is None if pre-computed pair counts are not supplied.
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
    weights1 : None, optional
        Weights for the first set of points. Default is None if not weighting pair counts by the selection function. 
    weights2 : None, optional
        Weights for the second set of points. Default is None if not weighting pair counts by the selection function.
        
    See Also
    --------
    quaia.recenter : Used to compute the bin centers of `thetabins`.
    quaia.read: Used to obtain the pixel indices and mask for `tab_gcatlo` and `tab_randlo`.
    """
    
    if weight_type == 'pair_product':
        weights1 = 1/selfunc_lo[pixel_indices_gcatlo]
        weights2 = 1/selfunc_lo[pixel_indices_randlo]
        
    # comoving distance
    DD_counts, api_time = mocks.DDtheta_mocks(autocorr = 1, nthreads = nthreads, binfile = thetabins, RA1 = tab_gcatlo['ra'], 
                                          DEC1 = tab_gcatlo['dec'], weights1 = weights1,
                                          weight_type=weight_type, c_api_timer = True)
    
    # now measure clustering in random catalog
    DR_counts, api_time = mocks.DDtheta_mocks(autocorr = 0,nthreads = nthreads, binfile = thetabins, RA1 = tab_gcatlo['ra'], 
                                          DEC1 = tab_gcatlo['dec'], weights1 = weights1, 
                                          RA2 = tab_randlo['ra'], DEC2 = tab_randlo['dec'], 
                                          weights2 = weights2, weight_type=weight_type, 
                                          c_api_timer = True)
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDtheta_mocks(autocorr = 1, nthreads = nthreads, binfile = thetabins, RA1 = tab_randlo['ra'], 
                                         DEC1 = tab_randlo['dec'], weights1 = weights2, 
                                         weight_type=weight_type, c_api_timer = True)
        RR_counts['thetaavg'] = np.mean([RR_counts['thetamin'], RR_counts['thetamax']], axis = 0)
        
        if weight_type == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
    
    # compute bin centers for theta
    DD_counts['thetaavg'] = np.mean([DD_counts['thetamin'], DD_counts['thetamax']], axis = 0)
    DR_counts['thetaavg'] = np.mean([DR_counts['thetamin'], DR_counts['thetamax']], axis = 0)

    # compute weighted pair counts
    if weight_type == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    return convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)

# 3D projected clustering wp(rp)
def wp_rp(tab_gcatlo, tab_randlo, N_gcatlo, N_randlo, key, selfunc_lo = None, pixel_indices_gcatlo = None, pixel_indices_randlo = None, weight_type = None, 
          weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, rbins = np.logspace(np.log10(0.5), np.log10(60.0), 21), nbins = 20, pimax = 40.0, 
          d = 1):
    
    r"""Computes the two-point projected 3D clustering of `tab_gcatlo`.
    
    Returns the projected correlation function computed as a function of projected quasar separation. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    N_gcatlo : int
        The number of objects in `tab_gcatlo`.
    N_randlo : int
        The number of objects in `tab_randlo`.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    pixel_indices_gcatlo : astropy table column, optional
        Pixel indices for each object in `tab_gcatlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    pixel_indices_randlo : astropy table column, optional 
        Pixel indices for each object in `tab_randlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    weight_type : {None, 'pair_product'}, optional
        Whether or not to weight pair counts by the selection function. Default is None.
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
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
        
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
    
    if weight_type == 'pair_product':
        weights1 = 1/selfunc_lo[pixel_indices_gcatlo][::d]
        weights2 = 1/selfunc_lo[pixel_indices_randlo][::d]
        
    # comoving distance 
    DD_counts, api_time = mocks.DDrppi_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'],  # where hubble distance = c/H0 and H0 = 100 km/s/Mpc h
                                         CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1,
                                         is_comoving_dist = True, weight_type=weight_type, output_rpavg = True, c_api_timer = True) 
    
    DR_counts, api_time = mocks.DDrppi_mocks(autocorr = 0, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                         CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1,
                                         RA2 = tab_randlo['ra'][::d], DEC2 = tab_randlo['dec'][::d],
                                         CZ2 = comoving_dist(tab_randlo['redshift_quaia_'+key][::d]), 
                                         weights2 = weights2, weight_type=weight_type,
                                         is_comoving_dist = True, output_rpavg = True, c_api_timer = True)

    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDrppi_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, pimax = pimax, binfile = rbins, 
                                         RA1 = tab_randlo['ra'][::d], DEC1 = tab_randlo['dec'][::d], 
                                         CZ1 = comoving_dist(tab_randlo['redshift_quaia_'+key][::d]), 
                                         weights1 = weights2, weight_type=weight_type,
                                         is_comoving_dist = True, output_rpavg = True, c_api_timer = True)
        if weight_type == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']

    # compute weighted pair counts
    if weight_type == 'pair_product':
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
def xi_r(tab_gcatlo, tab_randlo, N_gcatlo, N_randlo, key, selfunc_lo = None, pixel_indices_gcatlo = None, pixel_indices_randlo = None, weight_type = None,
         weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, rbins = np.logspace(np.log10(0.1), np.log10(20.0), 21), 
         correction = -5205.182232081812):
    
    r"""Computes the two-point  3D clustering of `tab_gcatlo`.
    
    Returns the correlation function computed as a function of quasar separation. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    N_gcatlo : int
        The number of objects in `tab_gcatlo`.
    N_randlo : int
        The number of objects in `tab_randlo`.
    key : str
        The suffix of the key for redshifts assigned to `tab_randlo`, created under the desired masking scheme using quaia.z_dist.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    pixel_indices_gcatlo : astropy table column, optional
        Pixel indices for each object in `tab_gcatlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    pixel_indices_randlo : astropy table column, optional 
        Pixel indices for each object in `tab_randlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    weight_type : {None, 'pair_product'}, optional
        Whether or not to weight pair counts by the selection function. Default is None.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default option is None if pre-computed pair counts are not supplied.
    nthreads: int, optional
        The number of OpenMP threads to use for computing each set of pair counts. Default is 8.
    rbins : array_like, optional
        The separation bins. Default is ``np.logspace(np.log10(0.1), np.log10(20.0), 21)``.
    correction: {int, float}, optional
        The distance, in Mpc/h, by which to shift the catalog to ensure that it exists in a box of coordinates [0, boxsize].
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
        
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
    
    if weight_type == 'pair_product':
        weights1 = 1/selfunc_lo[pixel_indices_gcatlo]
        weights2 = 1/selfunc_lo[pixel_indices_randlo]

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.     
    c = SkyCoord(ra=tab_gcatlo['ra'].value*u.degree, dec=tab_gcatlo['dec'].value*u.degree, distance=comoving_dist(tab_gcatlo['redshift_quaia']))
    X1, Y1, Z1 = c.cartesian.xyz.value - correction 
    DD_counts, api_time = theory.DD(autocorr = 1, nthreads = nthreads, binfile = rbins, periodic = False,
                                    X1 = X1, Y1 = Y1, Z1 = Z1, weights1 = weights1,
                                    weight_type=weight_type, output_ravg = True, c_api_timer = True) # cz/H0 = Mpc/h = m/s * km/1000 m / (100 km/s/Mpc h)

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.     
    c = SkyCoord(ra=tab_randlo['ra'], dec=tab_randlo['dec'], distance=comoving_dist(tab_randlo['redshift_quaia_'+key]))
    X2, Y2, Z2 = c.cartesian.xyz.value - correction
    DR_counts, api_time = theory.DD(autocorr = 0, nthreads = nthreads, binfile = rbins, periodic = False,
                                    X1 = X1, Y1 = Y1, Z1 = Z1, weights1 = weights1, 
                                    X2 = X2, Y2 = Y2, Z2 = Z2, weights2 = weights2,
                                    weight_type=weight_type, output_ravg = True, c_api_timer = True)
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = theory.DD(autocorr = 1, nthreads = nthreads, binfile = rbins, periodic = False,
                                X1 = X2, Y1 = Y2, Z1 = Z2, weights1 = weights2, 
                                weight_type=weight_type, output_ravg = True, c_api_timer = True)
        if weight_type == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
        
    # compute weighted pair counts
    if weight_type == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    cf = convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)
    
    # make it easier to plot
    ind = np.argsort(DD_counts['ravg'])
    
    return cf[ind], DD_counts['ravg'][ind]

# 3d clustering xi(r) computed with 1 mu bin
def xi_s(tab_gcatlo, tab_randlo, N_gcatlo, N_randlo, key, selfunc_lo = None, pixel_indices_gcatlo = None, pixel_indices_randlo = None,
         weight_type = None, weights1 = None, weights2 = None, RR_counts = None, nthreads = 8, rbins = np.logspace(np.log10(0.1), np.log10(20.0), 21)):
    
    r"""Computes the two-point  3D clustering of `tab_gcatlo`.
    
    Returns the correlation function computed as a function of quasar separation in redshift-space. Uses Corrfunc to compute the pair counts first, and then convert pair counts to clustering via the Landy-Szalay estimator. Provides the option to supply pre-computed random-random pair counts to speed up the computation. Leave `RR_counts` as None to avoid this option. 
    
    Parameters
    ----------
    tab_gcatlo : astropy table
        The quasar catalog. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    tab_randlo : astropy table
        The random catalog, used for computing the clustering from pair counts. Pass in a masked catalog as ``tab_randlo[mask_randlo]``. Obtain the mask using quaia.read.
    N_gcatlo : int
        The number of objects in `tab_gcatlo`.
    N_randlo : int
        The number of objects in `tab_randlo`.
    key : str
        The suffix of the key for redshifts assigned to `tab_randlo`, created under the desired masking scheme using quaia.z_dist.
    selfunc_lo : array_like, optional
        The selection function. Default is None if not weighting pair counts by the selection function.
    pixel_indices_gcatlo : astropy table column, optional
        Pixel indices for each object in `tab_gcatlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    pixel_indices_randlo : astropy table column, optional 
        Pixel indices for each object in `tab_randlo`, obtained using quaia.read. Default is None if not weighting pair counts by the selection function.
    weight_type : {None, 'pair_product'}, optional
        Whether or not to weight pair counts by the selection function. Default is None.
    RR_counts : array_like, optional 
        The number of random-random pair counts in `tab_randlo`. Default option is None if pre-computed pair counts are not supplied.
    nthreads: int, optional
        The number of OpenMP threads to use for computing each set of pair counts. Default is 8.
    rbins : array_like, optional
        The separation bins. Default is ``np.logspace(np.log10(0.1), np.log10(20.0), 21)``.
    mask_type : {None, 'selfunc'}, optional
        Choose either a galactic latitude-based mask or a selection function-based mask. Default is None if galactic latitude-based mask is desired. Pass in ``'selfunc'`` if using a selection function-based mask.
        
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
    
    if weight_type == 'pair_product':
        weights1 = 1/selfunc_lo[pixel_indices_gcatlo]
        weights2 = 1/selfunc_lo[pixel_indices_randlo]

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.   
    DD_counts, api_time = mocks.DDsmu_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, binfile = rbins, mu_max = 1, nmu_bins = 1, 
                                            RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                            CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1, 
                                            weight_type = weight_type, output_savg = True, is_comoving_dist=True, c_api_timer = True)

    # obtain r: The comoving distance along the line-of-sight between two objects remains constant with time for objects in the Hubble flow.   
    DR_counts, api_time = mocks.DDsmu_mocks(autocorr = 0, cosmology = 2, nthreads = nthreads, binfile = rbins, mu_max = 1, nmu_bins = 1, 
                                            RA1 = tab_gcatlo['ra'], DEC1 = tab_gcatlo['dec'], 
                                            CZ1 = comoving_dist(tab_gcatlo['redshift_quaia']), weights1 = weights1, 
                                            RA2 = tab_randlo['ra'], DEC2 = tab_randlo['dec'], 
                                            CZ2 = comoving_dist(tab_randlo['redshift_quaia_'+key]), weights2 = weights2, 
                                            weight_type = weight_type, output_savg = True, is_comoving_dist = True, c_api_timer = True)
    
    # now measure clustering in random catalog
    if RR_counts == None:
        
        RR_counts, api_time = mocks.DDsmu_mocks(autocorr = 1, cosmology = 2, nthreads = nthreads, binfile = rbins, mu_max = 1, nmu_bins = 1,
                                                RA1 = tab_randlo['ra'], DEC1 = tab_randlo['dec'], 
                                                CZ1 = comoving_dist(tab_randlo['redshift_quaia_'+key]), weights1 = weights2, 
                                                weight_type = weight_type, output_savg = True, is_comoving_dist=True, c_api_timer = True)
        if weight_type == 'pair_product':
            RR_counts['npairs'] = RR_counts['npairs']*RR_counts['weightavg']
        
    # compute weighted pair counts
    if weight_type == 'pair_product':
        DD_counts['npairs'] = DD_counts['npairs']*DD_counts['weightavg']
        DR_counts['npairs'] = DR_counts['npairs']*DR_counts['weightavg']
    
    # All the pair counts are done, get the angular correlation function
    cf = convert_3d_counts_to_cf(N_gcatlo, N_gcatlo, N_randlo, N_randlo, DD_counts, DR_counts, DR_counts, RR_counts)
    
    # make it easier to plot
    ind = np.argsort(DD_counts['savg'])
    
    # calculate the x axis
    rpavg = [np.sum((DD_counts['savg']*DD_counts['npairs'])[DD_counts['smin']==i])/np.sum(DD_counts['npairs'][DD_counts['smin']==i]) 
         for i in rbins[:-1]]
    
    # return cf[ind], DD_counts['savg'][ind], cf, DD_counts['savg'], rpavg
    return cf