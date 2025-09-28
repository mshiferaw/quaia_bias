import numpy as np
import quaia
import sys

def main():
    
    nthreads = 20
    z_array = [0,1]
    # thetabins = np.logspace(np.log10(0.1), np.log10(35), 15) #(np.log10(0.15), np.log10(35), 15)

    # Create the bins array
    rmin = 0.15 # 0.5 # 0.1 # start higher
    rmax = 240 # 120 #85 # 60.0 # 20.0 
    nbins = 20 
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    pimax = 80 # 40 # 5000.0 # 40.0
    
    G = float(sys.argv[1])
    zbin = int(sys.argv[2])
    i = int(sys.argv[3])
    fac_rand = 10

    tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _, _, _ = quaia.make_bins(G, zbin, False, 'z', tab_gcat_type = 'mocks', method = 'split', n_bins = 2, fac_rand = fac_rand, i = i+1, mask_type = None, b = 20)
    
    wp_datahi_zbin0, rpavg_datahi_zbin0 = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = rbins,
                                                                   nbins = nbins, pimax = pimax)
    np.save('../results/wp_G{:.1f}_zsplit2bin{}_b20_mock{}_{}x'.format(G, zbin, i+1, fac_rand), wp_datahi_zbin0)
    np.save('../results/rpavg_G{:.1f}_zsplit2bin{}_b20_mock{}_{}x'.format(G, zbin, i+1, fac_rand), rpavg_datahi_zbin0)
    
    # Create the bins array
    rmin = 0.1 # 0.1 # start higher
    rmax = 180.0 # 20.0 #50.0
    nbins = 20 #90 #20
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    
    cf_datahi_zbin0 = quaia.xi_s(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = rbins)
    np.save('../results/xi_G{:.1f}_zsplit2bin{}_b20_mock{}_{}x'.format(G, zbin, i+1, fac_rand), cf_datahi_zbin0)
                
if __name__=='__main__':
    main()