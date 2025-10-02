import numpy as np
import quaia
import sys
    
def main():
    
    nthreads = 18
    
    G = float(sys.argv[1])
    zbin = sys.argv[2]
    fac_rand = 10
    tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, _ = quaia.make_bins(G, zbin, True, 'z', method = 'split', n_bins = 2, tab_gcat_type = 'data', z = True, fac_rand = fac_rand, mask_type = None, b = 20)
    
    # Create the bins array
    rmin = 0.1 # 0.1 # start higher
    rmax = 180.0 # 20.0 #50.0
    nbins = 20 #90 #20
    sbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    
    cf_datahi_zbin0 = quaia.xi_s(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = sbins)
    np.save('../results/xi_G{:.1f}_zsplit2bin{}_b20_{}x'.format(G, zbin, fac_rand), cf_datahi_zbin0)
    
    # Create the bins array
    rmin = 0.15 # 0.5 # 0.1 # start higher
    rmax = 240 # 120 #85 # 60.0 # 20.0 
    nbins = 20 
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    pimax = 80 # 40 # 5000.0 # 40.0
    
    wp_datahi_zbin0, rpavg_datahi_zbin0 = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = rbins,
                                                                   nbins = nbins, pimax = pimax)
    np.save('../results/wp_G{:.1f}_zsplit2bin{}_b20_{}x'.format(G, zbin, fac_rand), wp_datahi_zbin0)
    np.save('../results/rpavg_G{:.1f}_zsplit2bin{}_b20_{}x'.format(G, zbin, fac_rand), rpavg_datahi_zbin0)
    
    # no binning
    if zbin=='0':
        tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, _ = quaia.make_bins(G, zbin, True, 'z', tab_gcat_type = 'data', z = True, fac_rand = fac_rand, mask_type = None, b = 20)
        cf_datahi_zbin0 = quaia.xi_s(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = sbins)
        np.save('../results/xi_G{:.1f}_b20_{}x'.format(G, fac_rand), cf_datahi_zbin0)

        wp_datahi_zbin0, rpavg_datahi_zbin0 = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = rbins,
                                                                       nbins = nbins, pimax = pimax)
        np.save('../results/wp_G{:.1f}_b20_{}x'.format(G, fac_rand), wp_datahi_zbin0)
        np.save('../results/rpavg_G{:.1f}_b20_{}x'.format(G, fac_rand), rpavg_datahi_zbin0)
    

if __name__=='__main__':
    main()