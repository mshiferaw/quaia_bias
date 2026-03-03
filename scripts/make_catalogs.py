import numpy as np
from astropy.table import Table, join
import utils
import sys
import argparse
import quaia

def main():
    overwrite = True

    ### Make catalogs with G-cut and redshifts
    # tag_qspec = ''
    # #tag_cat = '_mags-0.05'
    # tag_cat = ''
    # G_maxs = [20.0, 20.5, 20.6]
    # #G_maxs = [20.6]
    # for G_max in G_maxs:
    #     merge_gaia_spzs_and_cutGmax(G_max=G_max, tag_qspec=tag_qspec, tag_cat=tag_cat, overwrite=overwrite)

    # ### Make public-facing catalog
    #tag_qspec = ''
    #tag_cat = '_mags-0.05'
    #tag_cat = ''
    #G_maxs = [20.0, 20.5]
    #for G_max in G_maxs:
    #    make_public_catalog(G_max=G_max, tag_qspec=tag_qspec, tag_cat=tag_cat, overwrite=overwrite)

    # ### Make redshift-split catalogs
    #G_max = 20.5
    #n_zbins = 4
    #make_redshift_split_catalogs(G_max, n_zbins)

    ### Make redshift-split catalogs for CIB analysis by Giulia
    # G_max = 20.5
    # z_bins = [0, 1.0, 2.3, 5]
    # make_redshift_split_catalogs(G_max, z_bins=z_bins, save_tag='CIB')
    # z_bins = [0,  0.5, 1.0, 1.5, 2.0, 2.5, 5]
    # make_redshift_split_catalogs(G_max, z_bins=z_bins, save_tag='CIB')

    ### Make redshift-split catalogs for autocorr-dutycycle analysis by christina eilers & mariona
    # G_max = 20.5
    #z_bins = [0.0,1.0,2.0,3.0,4.0]
    #make_redshift_split_catalogs(G_max, z_bins=z_bins)
    # z_bins = [2.9,3.5,5.0]
    # make_redshift_split_catalogs(G_max, z_bins=z_bins)
    
    parser = argparse.ArgumentParser()
    parser.add_argument("G", help="magnitude cut", type = float)
    parser.add_argument("--n_zbins", help="number of z bins", type = int)
    parser.add_argument("--z_bins", help="z bins", type = float, nargs = '*')
    parser.add_argument("--n_Lbins", help="number of L bins", type = int)
    parser.add_argument("--L_bins", help="L bins", type = float, nargs = '*')
    parser.add_argument("--L_bolmin", help="L_bol threshold", type = float)
    args = parser.parse_args()
    # if args.n_zbins:
    #     fn_gcat_zbin = make_redshift_split_catalogs(args.G, n_zbins=int(args.n_zbins))
    # elif args.z_bins:
    #     # print(args.z_bins)
    #     fn_gcat_zbin = make_redshift_split_catalogs(args.G, z_bins=args.z_bins)
    # elif args.n_Lbins:
    #     make_luminosity_split_catalogs(args.G, n_Lbins=int(args.n_Lbins))
    # elif args.L_bins:
    #     make_luminosity_split_catalogs(args.G, L_bins=args.L_bins)
    # fn_gcat_zbin = None
    if (args.n_zbins or args.z_bins) and not (args.n_Lbins or args.L_bins or args.L_bolmin):
        make_redshift_split_catalogs(args.G, n_zbins=args.n_zbins, z_bins=args.z_bins)
    elif not (args.n_zbins or args.z_bins) and (args.n_Lbins or args.L_bins):
    # if args.n_Lbins or args.L_bins:
        make_luminosity_split_catalogs(args.G, n_Lbins=args.n_Lbins, L_bins=args.L_bins)
    # joint split
    # if args.n_Lbins or args.L_bins:
    #     make_redshift_luminosity_split_catalogs(fn_gcat_zbin, args.G, n_Lbins=int(args.n_Lbins), L_bins=args.L_bins)
    elif (args.n_zbins or args.z_bins) and (args.n_Lbins or args.L_bins):
        make_redshift_luminosity_split_catalogs(args.G, n_zbins=args.n_zbins, z_bins=args.z_bins, n_Lbins=args.n_Lbins, L_bins=args.L_bins)
    elif (args.n_zbins or args.z_bins) and (args.L_bolmin):
        make_redshift_luminosity_threshold_catalogs(args.G, n_zbins=args.n_zbins, z_bins=args.z_bins, L_bolmin=args.L_bolmin)    

    # zbins = sys.argv[2]
    # print(zbins)
    # print(type(zbins))
    # if type(zbins)==list:
    #     make_redshift_split_catalogs(sys.argv[1], z_bins=zbins)
    # else:
    #     make_redshift_split_catalogs(sys.argv[1], n_zbins=int(zbins))
    
    # make_redshift_split_catalogs(sys.argv[1], n_zbins=int(sys.argv[2]))


def merge_gaia_spzs_and_cutGmax(G_max=20.5, tag_qspec='', tag_cat='', overwrite=False):

    # save name
    fn_gcat = f'../data/catalog_G{G_max}{tag_qspec}{tag_cat}.fits'

    # data paths
    fn_gaia = f'../data/gaia_candidates_clean{tag_qspec}{tag_cat}.fits'
    fn_spz = f'../data/redshift_estimates/redshifts_spz{tag_qspec}{tag_cat}_kNN_K27_std.fits'

    # load data, cut to G_max
    tab_gaia = utils.load_table(fn_gaia)
    tab_gaia = tab_gaia[tab_gaia['phot_g_mean_mag'] < G_max]

    # SPZ-only table
    tab_spz = utils.load_table(fn_spz)
    tab_spz.keep_columns(['source_id', 'redshift_spz', 'redshift_spz_raw', 'redshift_spz_err'])

    tab_gcat = join(tab_gaia, tab_spz, keys='source_id', join_type='inner')
    utils.add_randints_column(tab_gcat)
    tab_gcat.write(fn_gcat, overwrite=overwrite)
    print(f"Wrote table with {len(tab_gcat)} objects to {fn_gcat}")



def make_public_catalog(G_max=20.5, tag_qspec='', tag_cat='', overwrite=False):

    # working catalog
    fn_gcat = f'../data/catalog_G{G_max}{tag_qspec}{tag_cat}.fits'
    # update to final name choice!
    fn_public = f'../data/quaia_G{G_max}{tag_qspec}{tag_cat}.fits'

    tab_gcat = utils.load_table(fn_gcat)

    columns_to_keep = ['source_id', 'unwise_objid', 
                       'redshift_spz', 'redshift_spz_err', 
                       'ra', 'dec', 'l', 'b', 
                       'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 
                       'mag_w1_vg', 'mag_w2_vg', 
                       'pm', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error']

    tab_public = Table()
    tab_public.meta = {'name': '\emph{{Gaia}}--\emph{{unWISE}} Quasar Catalog',
                       'abbrv': 'Quaia'
                       }

    rename_dict = {'redshift_spz': 'redshift_quaia',
                   'redshift_spz_err': 'redshift_quaia_err'
                   }

    for cn in columns_to_keep:
        if cn in rename_dict:
            cn_new = rename_dict[cn]
        else: 
            cn_new = cn
        tab_public[cn_new] = tab_gcat[cn]
        tab_public[cn_new].info.unit = utils.label2unit_dict[cn_new]
        tab_public[cn_new].info.description = utils.label2description_dict[cn_new]
    
    # for tc in tab_public.columns:
    #     print(tc, tab_public[tc].info.unit)
    print(tab_public.columns)
    tab_public.write(fn_public, overwrite=overwrite)
    print(f"Wrote table with {len(tab_public)} objects to {fn_public}")


def make_redshift_split_catalogs(G_max, n_zbins=None, z_bins=None, overwrite=True,
                                 save_tag='', z_split = True):

    assert n_zbins is not None or z_bins is not None, "Either n_zbins or z_bins must be passed!"

    if z_bins is not None and n_zbins is not None:
        print("z_bins passed, ignoring n_zbins")
    if z_bins is not None:
        n_zbins = len(z_bins)-1
        print(f"z_bins: {z_bins}, setting n_zbins={n_zbins}")
        z_split = False

    fn_gcat = f'../data/quaia_G{G_max}.fits'
    tab_gcat = utils.load_table(fn_gcat)

    if z_bins is None:
        n_zbins = int(n_zbins)
        z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
        print(z_percentiles)
        z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
        z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
        z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

    print("zbins:", z_bins)
    print("n_zbins:", n_zbins)

    for bb in range(n_zbins):
        i_zbin = (tab_gcat['redshift_quaia'] >= z_bins[bb]) & (tab_gcat['redshift_quaia'] < z_bins[bb+1])
        tab_gcat_zbin = tab_gcat[i_zbin]
        # if z_bins is None:
        if z_split == True:
            fn_gcat_zbin = f'../data/quaia_G{G_max}_zsplit{n_zbins}bin{bb}{save_tag}.fits'
        else:
            fn_gcat_zbin = f'../data/quaia_G{G_max}_zmin{z_bins[bb]}zmax{z_bins[bb+1]}{save_tag}.fits' 
            # fn_gcat_zbin = f'../data/quaia_G{G_max}_zsplit{n_zbins}bin{bb}{save_tag}.fits'
        tab_gcat_zbin.write(fn_gcat_zbin, overwrite=overwrite)
        print("zmin:", np.min(tab_gcat_zbin['redshift_quaia']))
        print("zmax:", np.max(tab_gcat_zbin['redshift_quaia']))
        print(f"Wrote table with {len(tab_gcat_zbin)} objects to {fn_gcat_zbin}")

    # return fn_gcat_zbin

def make_luminosity_split_catalogs(G_max, n_Lbins=None, L_bins=None, overwrite=True,
                                 save_tag='', L_split = True, fn_gcat = None, prefix = ''):

    assert n_Lbins is not None or L_bins is not None, "Either n_Lbins or L_bins must be passed!"

    if L_bins is not None and n_Lbins is not None:
        print("L_bins passed, ignoring n_Lbins")
    if L_bins is not None:
        n_Lbins = len(L_bins)-1
        print(f"L_bins: {L_bins}, setting n_Lbins={n_Lbins}")
        L_split = False

    if fn_gcat is None:
        fn_gcat = f'../data/quaia_G{G_max}.fits'
    else:
        prefix = fn_gcat[19:-5]
        print(prefix, 'already computed')
    tab_gcat = utils.load_table(fn_gcat)

    quaia.absolute(tab_gcat, G_max)
    if L_bins is None:
        L_percentiles = np.linspace(0.0, 100.0, n_Lbins+1)
        print(L_percentiles)
        L_bins = np.percentile(list(tab_gcat['M_i']), L_percentiles)
        L_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
        L_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

    print("Lbins:", L_bins)
    print("n_Lbins:", n_Lbins)

    for bb in range(n_Lbins):
        i_Lbin = (tab_gcat['M_i'] >= L_bins[bb]) & (tab_gcat['M_i'] < L_bins[bb+1])
        tab_gcat_Lbin = tab_gcat[i_Lbin]
        # if L_bins is None:
        if L_split == True:
            fn_gcat_Lbin = f'../data/quaia_G{G_max}{prefix}_Lsplit{n_Lbins}bin{bb}{save_tag}.fits'
        else:
            fn_gcat_Lbin = f'../data/quaia_G{G_max}{prefix}_Lmin{L_bins[bb]}Lmax{L_bins[bb+1]}{save_tag}.fits' 
        tab_gcat_Lbin.write(fn_gcat_Lbin, overwrite=overwrite)
        print("Lmin:", np.min(tab_gcat_Lbin['M_i']))
        print("Lmax:", np.max(tab_gcat_Lbin['M_i']))
        print(f"Wrote table with {len(tab_gcat_Lbin)} objects to {fn_gcat_Lbin}")
        
# def make_redshift_luminosity_split_catalogs(fn_gcat_zbin, G_max, n_Lbins=None, L_bins=None, overwrite=True, save_tag='', L_split = True):

#     # option 1: write function that does it separately
#     # option 2: merge into one func
#     # option 3: write function that loads it if it's already split in z first
    
#     assert n_Lbins is not None or L_bins is not None, "Either n_Lbins or L_bins must be passed!"

#     if L_bins is not None and n_Lbins is not None:
#         print("L_bins passed, ignoring n_Lbins")
#     if L_bins is not None:
#         n_Lbins = len(L_bins)-1
#         print(f"L_bins: {L_bins}, setting n_Lbins={n_Lbins}")
#         L_split = False

#     # fn_gcat = f'../data/quaia_G{G_max}.fits'
#     tab_gcat = utils.load_table(fn_gcat_zbin)

#     quaia.absolute(tab_gcat, G_max)
#     if L_bins is None:
#         L_percentiles = np.linspace(0.0, 100.0, n_Lbins+1)
#         print(L_percentiles)
#         L_bins = np.percentile(list(tab_gcat['M_i']), L_percentiles)
#         L_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
#         L_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

#     print("Lbins:", L_bins)
#     print("n_Lbins:", n_Lbins)

#     for bb in range(n_Lbins):
#         i_Lbin = (tab_gcat['M_i'] >= L_bins[bb]) & (tab_gcat['M_i'] < L_bins[bb+1])
#         tab_gcat_Lbin = tab_gcat[i_Lbin]
#         if L_split == True:
#             fn_gcat_Lbin = fn_gcat_zbin[:-5]+f'_Lsplit{n_Lbins}bin{bb}{save_tag}.fits'
#         else:
#             fn_gcat_Lbin = fn_gcat_zbin[:-5]+f'_Lmin{L_bins[bb]}Lmax{L_bins[bb+1]}{save_tag}.fits' 
#         tab_gcat_Lbin.write(fn_gcat_Lbin, overwrite=overwrite)
#         print("Lmin:", np.min(tab_gcat_Lbin['M_i']))
#         print("Lmax:", np.max(tab_gcat_Lbin['M_i']))
#         print(f"Wrote table with {len(tab_gcat_Lbin)} objects to {fn_gcat_Lbin}")
        
    # option 1
def make_redshift_luminosity_split_catalogs(G_max, n_zbins=None, z_bins=None, n_Lbins=None, L_bins=None, overwrite=True, save_tag='', z_split = True, L_split = True):
    
    assert n_zbins is not None or z_bins is not None, "Either n_bins or z_bins must be passed!"
    assert n_Lbins is not None or L_bins is not None, "Either n_Lbins or L_bins must be passed!"

    if z_bins is not None and n_zbins is not None:
        print("z_bins passed, ignoring n_zbins")
    if z_bins is not None:
        n_zbins = len(z_bins)-1
        print(f"z_bins: {z_bins}, setting n_zbins={n_zbins}")
        z_split = False
    if L_bins is not None and n_Lbins is not None:
        print("L_bins passed, ignoring n_Lbins")
    if L_bins is not None:
        n_Lbins = len(L_bins)-1
        print(f"L_bins: {L_bins}, setting n_Lbins={n_Lbins}")
        L_split = False

    fn_gcat = f'../data/quaia_G{G_max}.fits'
    tab_gcat = utils.load_table(fn_gcat)
    
    if z_bins is None:
        n_zbins = int(n_zbins)
        z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
        print(z_percentiles)
        z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
        z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
        z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

    print("zbins:", z_bins)
    print("n_zbins:", n_zbins)

    for bb in range(n_zbins):
        i_zbin = (tab_gcat['redshift_quaia'] >= z_bins[bb]) & (tab_gcat['redshift_quaia'] < z_bins[bb+1])
        tab_gcat_zbin = tab_gcat[i_zbin]
        quaia.absolute(tab_gcat_zbin, G_max)
        if z_split == True:
            fn_gcat_zbin = f'_zsplit{n_zbins}bin{bb}{save_tag}'
        else:
            fn_gcat_zbin = f'_zmin{z_bins[bb]}zmax{z_bins[bb+1]}{save_tag}' 
        for bb in range(n_Lbins):
            if L_bins is None:
                n_Lbins = int(n_Lbins)
                L_percentiles = np.linspace(0.0, 100.0, n_Lbins+1)
                print(L_percentiles)
                L_bins = np.percentile(list(tab_gcat['M_i']), L_percentiles)
                L_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
                L_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included
            print("Lbins:", L_bins)
            print("n_Lbins:", n_Lbins)
            i_Lbin = (tab_gcat_zbin['M_i'] >= L_bins[bb]) & (tab_gcat_zbin['M_i'] < L_bins[bb+1])
            tab_gcat_Lbin = tab_gcat_zbin[i_Lbin]
            if L_split == True:
                fn_gcat_Lbin = f'../data/quaia_G{G_max}{fn_gcat_zbin}_Lsplit{n_zbins}bin{bb}{save_tag}.fits'
            else:
                fn_gcat_Lbin = f'../data/quaia_G{G_max}{fn_gcat_zbin}_Lmin{L_bins[bb]}Lmax{L_bins[bb+1]}{save_tag}.fits' 
            tab_gcat_Lbin.write(fn_gcat_Lbin, overwrite=overwrite)
            try:
                print("zmin:", np.min(tab_gcat_Lbin['redshift_quaia']))
                print("zmax:", np.max(tab_gcat_Lbin['redshift_quaia']))
                print("Lmin:", np.min(tab_gcat_Lbin['M_i']))
                print("Lmax:", np.max(tab_gcat_Lbin['M_i']))
            except ValueError:  #raised if `y` is empty.
                pass
            print(f"Wrote table with {len(tab_gcat_Lbin)} objects to {fn_gcat_Lbin}")
            
def make_redshift_luminosity_threshold_catalogs(G_max, n_zbins=None, z_bins=None, n_Lbins=None, L_bins=None, overwrite=True, save_tag='', z_split = True, L_split = True, Lbol_threshold = True, L_bolmin = None, key = 'M_i'):
    
    assert n_zbins is not None or z_bins is not None, "Either n_bins or z_bins must be passed!"
    assert n_Lbins is not None or L_bins is not None or L_bolmin is not None, "Either n_Lbins or L_bins or L_bolmin must be passed!"

    if z_bins is not None and n_zbins is not None:
        print("z_bins passed, ignoring n_zbins")
    if z_bins is not None:
        n_zbins = len(z_bins)-1
        print(f"z_bins: {z_bins}, setting n_zbins={n_zbins}")
        z_split = False
    if L_bins is not None and n_Lbins is not None:
        print("L_bins passed, ignoring n_Lbins")
    if L_bins is not None:
        n_Lbins = len(L_bins)-1
        print(f"L_bins: {L_bins}, setting n_Lbins={n_Lbins}")
        L_split = False
    if L_bolmin is not None:
        n_Lbins = 1
        L_split = False
        L_bol_threshold = True

    fn_gcat = f'../data/quaia_G{G_max}.fits'
    tab_gcat = utils.load_table(fn_gcat)
    
    if z_bins is None:
        n_zbins = int(n_zbins)
        z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
        print(z_percentiles)
        z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
        z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
        z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

    print("zbins:", z_bins)
    print("n_zbins:", n_zbins)

    for bb in range(n_zbins):
        i_zbin = (tab_gcat['redshift_quaia'] >= z_bins[bb]) & (tab_gcat['redshift_quaia'] < z_bins[bb+1])
        tab_gcat_zbin = tab_gcat[i_zbin]
        quaia.absolute(tab_gcat_zbin, G_max)
        quaia.L_bol(tab_gcat_zbin)
        if z_split == True:
            fn_gcat_zbin = f'_zsplit{n_zbins}bin{bb}{save_tag}'
        else:
            fn_gcat_zbin = f'_zmin{z_bins[bb]}zmax{z_bins[bb+1]}{save_tag}' 
        for bb in range(n_Lbins):
            if L_bolmin is not None:
                L_bins = [10**L_bolmin, np.inf]
                key = 'L_bol'
            elif L_bins is None:
                n_Lbins = int(n_Lbins)
                L_percentiles = np.linspace(0.0, 100.0, n_Lbins+1)
                print(L_percentiles)
                L_bins = np.percentile(list(tab_gcat['M_i']), L_percentiles)
                L_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
                L_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included 
            print("Lbins:", L_bins)
            print("n_Lbins:", n_Lbins)
            print(key, tab_gcat_zbin[key])
            i_Lbin = (tab_gcat_zbin[key] >= L_bins[bb]) & (tab_gcat_zbin[key] < L_bins[bb+1])
            tab_gcat_Lbin = tab_gcat_zbin[i_Lbin]
            if L_split == True:
                fn_gcat_Lbin = f'../data/quaia_G{G_max}{fn_gcat_zbin}_Lsplit{n_Lbins}bin{bb}{save_tag}.fits'
            elif L_bol_threshold == True:
                fn_gcat_Lbin = f'../data/quaia_G{G_max}{fn_gcat_zbin}_Lbolmin{L_bolmin}{save_tag}.fits'
            else:
                fn_gcat_Lbin = f'../data/quaia_G{G_max}{fn_gcat_zbin}_Lmin{L_bins[bb]}Lmax{L_bins[bb+1]}{save_tag}.fits' 
            tab_gcat_Lbin.write(fn_gcat_Lbin, overwrite=overwrite)
            try:
                print("zmin:", np.min(tab_gcat_Lbin['redshift_quaia']))
                print("zmax:", np.max(tab_gcat_Lbin['redshift_quaia']))
                print("Lmin:", np.min(tab_gcat_Lbin[key]))
                print("Lmax:", np.max(tab_gcat_Lbin[key]))
            except ValueError:  #raised if `y` is empty.
                pass
            print(f"Wrote table with {len(tab_gcat_Lbin)} objects to {fn_gcat_Lbin}")
            
#     # option 2
# def make_redshift_luminosity_split_catalogs(G_max, n_zbins=None, z_bins=None, n_Lbins=None, L_bins=None, overwrite=True, save_tag='', z_split = True, L_split = True):
    
#     assert n_zbins is not None or z_bins is not None, "Either n_bins or z_bins must be passed!"
#     assert n_Lbins is not None or L_bins is not None, "Either n_Lbins or L_bins must be passed!"

#     if z_bins is not None and n_zbins is not None:
#         print("z_bins passed, ignoring n_zbins")
#     if z_bins is not None:
#         n_zbins = len(z_bins)-1
#         print(f"z_bins: {z_bins}, setting n_zbins={n_zbins}")
#         z_split = False
#     if L_bins is not None and n_Lbins is not None:
#         print("L_bins passed, ignoring n_Lbins")
#     if L_bins is not None:
#         n_Lbins = len(L_bins)-1
#         print(f"L_bins: {L_bins}, setting n_Lbins={n_Lbins}")
#         L_split = False

#     fn_gcat = f'../data/quaia_G{G_max}.fits'
#     tab_gcat = utils.load_table(fn_gcat)
    
#     if z_bins is None:
#         n_zbins = int(n_zbins)
#         z_percentiles = np.linspace(0.0, 100.0, n_zbins+1)
#         print(z_percentiles)
#         z_bins = np.percentile(list(tab_gcat['redshift_quaia']), z_percentiles)
#         z_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
#         z_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included

#     print("zbins:", z_bins)
#     print("n_zbins:", n_zbins)

#     for bb in range(n_zbins):
#         i_zbin = (tab_gcat['redshift_quaia'] >= z_bins[bb]) & (tab_gcat['redshift_quaia'] < z_bins[bb+1])
#         tab_gcat_zbin = tab_gcat[i_zbin]
#         quaia.absolute(tab_gcat_zbin, G_max)
#         for bb in range(n_Lbins):
#             if L_bins is None:
#                 L_percentiles = np.linspace(0.0, 100.0, n_Lbins+1)
#                 print(L_percentiles)
#                 L_bins = np.percentile(list(tab_gcat['M_i']), L_percentiles)
#                 L_bins[-1] += 0.01 # add a bit to maximum bin to make sure the highest-z source gets included
#                 L_bins[0] -= 0.01 # add a bit to minimum bin to make sure the lowest-z source gets included
#             i_Lbin = (tab_gcat_zbin['M_i'] >= L_bins[bb]) & (tab_gcat_zbin['M_i'] < L_bins[bb+1])
#             tab_gcat_Lbin = tab_gcat[i_Lbin]
#             if z_split == True:
#                 fn_gcat_zbin = f'../data/quaia_G{G_max}_zsplit{n_zbins}bin{bb}{save_tag}.fits'
#             else:
#                 fn_gcat_zbin = f'../data/quaia_G{G_max}_zmin{z_bins[bb]}zmax{z_bins[bb+1]}{save_tag}.fits' 
#             tab_gcat_zbin.write(fn_gcat_zbin, overwrite=overwrite)
#             print("zmin:", np.min(tab_gcat_zbin['redshift_quaia']))
#             print("zmax:", np.max(tab_gcat_zbin['redshift_quaia']))
#             print(f"Wrote table with {len(tab_gcat_zbin)} objects to {fn_gcat_zbin}")

if __name__=='__main__':
    main()
