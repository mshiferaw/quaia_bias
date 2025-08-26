import numpy as np
import quaia

def main():
    
    nthreads = 18
    z_array = [0,1]
    thetabins = np.logspace(np.log10(0.1), np.log10(35), 15) #(np.log10(0.15), np.log10(35), 15)
    
    for G in [20, 20.5]:
        for i in range(100):
            for zbin in z_array:

                tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, N_randhi_mask_zbin0, key_zbin0 = quaia.zbins(G, zbin, i=i+1, data = 'mocks')
                wtheta_datahi_zbin0 = quaia.w_theta(tab_datahi_mask_zbin0, tab_randhi_mask_zbin0, N_datahi_mask_zbin0, 
                                                                   N_randhi_mask_zbin0, nthreads = nthreads, thetabins = thetabins)
                np.save('../results/wtheta_G{}_zsplit2bin{}_mock{}'.format(G, zbin, i+1), wtheta_datahi_zbin0)
                
if __name__=='__main__':
    main()