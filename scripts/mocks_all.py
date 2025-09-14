import numpy as np
import quaia

def main():
    
    nthreads = 19
    z_array = [0,1]
    thetabins = np.logspace(np.log10(0.1), np.log10(35), 15) #(np.log10(0.15), np.log10(35), 15)

    # Create the bins array
    rmin = 0.15 # 0.5 # 0.1 # start higher
    rmax = 240 # 120 #85 # 60.0 # 20.0 
    nbins = 20 
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    pimax = 80 # 40 # 5000.0 # 40.0
    
    for G in [20, 20.5]:
        for i in range(100):
                tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _ = quaia.make_bins(G, 0, True, 'z', tab_gcat_type = 'mocks', fac_rand = 25, i = i+1)
                wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
                np.save('../results/wtheta_G{:.1f}_mock{}_25x'.format(G, i+1), wtheta_datahi_zbin0)

                wp_datahi_zbin0, rpavg_datahi_zbin0 = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = rbins,
                                                                   nbins = nbins, pimax = pimax)
                np.save('../results/wp_G{:.1f}_mock{}_25x'.format(G, i+1), wp_datahi_zbin0)
                np.save('../results/rpavg_G{:.1f}_mock{}_25x'.format(G, i+1), rpavg_datahi_zbin0)
                
if __name__=='__main__':
    main()