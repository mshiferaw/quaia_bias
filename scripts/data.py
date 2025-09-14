import numpy as np
import quaia
    
def main():
    
    nthreads = 18
    z_array = [0,1]
    thetabins = np.logspace(np.log10(0.1), np.log10(35), 15) #(np.log10(0.15), np.log10(35), 15)

    # Create the bins array
    rmin = 0.15 # 0.5 # 0.1 # start higher
    rmax = 240 # 120 #85 # 60.0 # 20.0 
    nbins = 20 
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1) # Mpc/h https://github.com/manodeep/Corrfunc/issues/202
    pimax = 80 # 40 # 5000.0 # 40.0

    for G in [20, 20.5]:
        for zbin in z_array:

            # tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0 = quaia.zbins(G, zbin)
            # wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
            #                                                    N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
            # print(G, N_datahi_mask_zbin0)
            # np.save('../results/wtheta_G{:.1f}_zsplit2bin{}'.format(G, zbin), wtheta_datahi_zbin0)
            
            # #more debugging
            # tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0 = quaia.make_zbins(G, zbin, False, method = 'zsplit', n_zbins = 2, tab_gcat_type = 'data')
            # wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
            #                                                    N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
            # print(G, N_datahi_mask_zbin0)
            # np.save('../results/wtheta_G{:.1f}_zsplit2bin{}_debug'.format(G, zbin), wtheta_datahi_zbin0)
            
            # #more debugging
            # tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0, z_bins = quaia.make_zbins(G, zbin, False, method = 'zsplit', n_zbins = 2, tab_gcat_type = 'data', z = False)
            # wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
            #                                                    N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
            # print(G, N_datahi_mask_zbin0)
            # np.save('../results/wtheta_G{:.1f}_zsplit2bin{}_debug2'.format(G, zbin), wtheta_datahi_zbin0)
            
            # tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0, z_bins = quaia.make_zbins(G, zbin, True, method = 'zsplit', n_zbins = 2, tab_gcat_type = 'data', z = False)
            # wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
            #                                                    N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
            # print(G, N_datahi_mask_zbin0, len(tab_datahi_mask_zbin0))
            # np.save('../results/wtheta_G{:.1f}_zsplit2bin{}_debug3'.format(G, zbin), wtheta_datahi_zbin0)

            tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, _, _ = quaia.make_bins(G, zbin, True, method = 'split', n_bins = 2, tab_gcat_type = 'data', z = False, fac_rand = 25)
            wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
            print(G, N_datahi_mask_zbin0, len(tab_datahi_mask_zbin0))
            np.save('../results/wtheta_G{:.1f}_zsplit2bin{}_25x'.format(G, zbin), wtheta_datahi_zbin0)

            wp_datahi_zbin0, rpavg_datahi_zbin0 = quaia.wp_rp(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, key_zbin0, nthreads = nthreads, rbins = rbins,
                                                           nbins = nbins, pimax = pimax)
            np.save('../results/wp_G{:.1f}_zsplit2bin{}_25x'.format(G, zbin), wp_datahi_zbin0)
            np.save('../results/rpavg_G{:.1f}_zsplit2bin{}_25x'.format(G, zbin), rpavg_datahi_zbin0)
            
#         # for debugging purposes...
#         tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0 = quaia.zbins(G, 0, n_zbins = 1)
#         wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
#                                                            N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
#         np.save('../results/wtheta_G{:.1f}'.format(G), wtheta_datahi_zbin0)
        
#         tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0 = quaia.zbins(G, 0, n_zbins = 1, tab_gcat = )
#         wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
#                                                            N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
#         np.save('../results/wtheta_G{:.1f}'.format(G), wtheta_datahi_zbin0)

if __name__=='__main__':
    main()