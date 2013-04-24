import sys
#sys.path.append('/user/bostroem/science/programs')
sys.path.append('/Users/bostroem/science/programs')
import normalize_spec
import pyfits
from matplotlib import pyplot
from stsci.convolve import boxcar
import pdb
import numpy as np

def make_plots(fuv, ccd1, ccd2, ccd3, ccd4, ax1_ylim, ax2_ylim):
    #Open data files
    tbdata1 = pyfits.getdata(ccd1, 1)
    tbdata2 = pyfits.getdata(ccd2, 1)
    tbdata3 = pyfits.getdata(ccd3, 1)
    tbdata4 = pyfits.getdata(ccd4, 1)
    tbdata5 = pyfits.getdata(fuv, 1)

    wl1 = tbdata1['wavelength'].ravel()
    net1 = tbdata1['net'].ravel()
    dq = tbdata1['dq'].ravel()
    wl1 = wl1[dq&512 != 512]
    net1 = net1[dq&512 != 512]

    wl2 = tbdata2['wavelength'].ravel()
    net2 = tbdata2['net'].ravel()
    dq = tbdata2['dq'].ravel()
    wl2 = wl2[dq&512 != 512]
    net2 = net2[dq&512 != 512]

    wl3 = tbdata3['wavelength'].ravel()
    net3 = tbdata3['net'].ravel()
    dq = tbdata3['dq'].ravel()
    wl3 = wl3[dq&512 != 512]
    net3 = net3[dq&512 != 512]

    wl4 = tbdata4['wavelength'].ravel()
    net4 = tbdata4['net'].ravel()
    dq = tbdata4['dq'].ravel()
    wl4 = wl4[dq&512 != 512]
    net4 = net4[dq&512 != 512]

    wl5 = tbdata5['wavelength'].ravel()
    net5 = tbdata5['net'].ravel()    
    dq = tbdata5['dq'].ravel()
    wl5 = wl5[dq&512 != 512]
    net5 = net5[dq&512 != 512]
    

    #normalize spectra
    net1 = normalize_spec.manually_remove_CR(wl1, net1)
    wl1, norm_spec1, cont1 = normalize_spec.normalize_spectrum(wl1, net1)
    norm_spec2 = normalize_spec.manually_remove_CR(wl2, net2)
    wl2, norm_spec2, cont2 = normalize_spec.normalize_spectrum(wl2, net2)
    norm_spec3 = normalize_spec.manually_remove_CR(wl3, net3)
    wl3, norm_spec3, cont3 = normalize_spec.normalize_spectrum(wl3, net3)
    norm_spec4 = normalize_spec.manually_remove_CR(wl4, net4)
    wl4, norm_spec4, cont4 = normalize_spec.normalize_spectrum(wl4, net4)
    wl5, norm_spec5, cont5 = normalize_spec.normalize_spectrum(wl5, net5, kernel_size = 401)



    #create plotting objects
    fig = pyplot.figure(1)
    fig.clf()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    ax2.plot(wl1, norm_spec1, 'r')
    ax2.plot(wl2, norm_spec2, 'g')
    ax2.plot(wl3, norm_spec3, 'm')
    ax2.plot(wl4, norm_spec4, 'b')

    ax3.plot(wl1, boxcar(norm_spec1, (3,)), 'r')
    ax3.plot(wl2, boxcar(norm_spec2, (3,)), 'g')
    ax3.plot(wl3, boxcar(norm_spec3, (3,)), 'm')
    ax3.plot(wl4, boxcar(norm_spec4, (3,)), 'b')

    ax1.plot(wl5, norm_spec5, 'k')

    ax2.set_xlim(3790, 4841)
    ax3.set_xlim(3790, 4841)
    ax1.set_ylim(ax1_ylim[0], ax1_ylim[1])
    ax2.set_ylim(ax2_ylim[0], ax2_ylim[1])
    ax3.set_ylim(ax2_ylim[0], ax2_ylim[1])

    ax1_xticks = np.arange(1150, 1750, 50)
    ax2_xticks = np.arange(3800, 4900, 50)
    ax1.set_xticks(ax1_xticks)
    ax2.set_xticks(ax2_xticks)
    ax3.set_xticks(ax2_xticks)
    
    ax1.grid(b = 'on', color = 'gray', ls = ':')
    ax2.grid(b = 'on', color = 'gray', ls = ':')
    ax2.grid(b = 'on', color = 'gray', ls = ':')

    #Label Plots
    ax2.set_title('Optical Spectrum H9, R136b')
    ax1.set_title('FUV Spectrum H9, R136b')
    ax3.set_title('Optical Spectrum H9, R136b R~4000')

    ax3.set_xlabel('Wavelength')
    ax2.set_ylabel('Relative Intensity')

    return ax1, ax2, ax3, fig


if __name__ == "__main__":
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_3936_combined_img_loc545.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4194_combined_img_loc495.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4451_combined_img_loc545.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4706_combined_img_loc542.fits'
    fuv = '/Users/bostroem/Dropbox/R136/mama/R136b_SE8loc459.fits'


    #line_dict = {'name': ['N IV 3480', 'Si IV 4089', 'Si IV 4116', 'N IV 4058', 'He II 4200', 'He I 4471', 'He II 4542', 'N III 4636', 'N III 4642', 'He II 4686', 'H-Beta', 'H-Gamma', 'H-Delta'], \
    #        'rest_wavelength': [3480, 4089, 4116, 4058, 4200, 4471, 4542, 4636, 4642, 4686, 4861.33, 4340.47, 4101.74]}

    line_dict = {'name':['N II 3995', 'He I 4009' ,'He I+II 4026', 'N II 4041/44', 'Ni IV 4058', 'C III 4068/69/70',\
         'O II 4070/72/76', 'Si IV 4089', 'N III 4097', 'H-Delta 4102','Si IV 4116' ,'He I 4121','He I 4144','C III 4187','He II 4200',\
         'O II 4254','C II 4267','O II 4276/85','O II 4317/20','H-Gamma 4340','O II 4349','O II 4367', 'N III 4379' ,'He I 4388',\
         'O II 4415/17','DIB is 4429','He I 4471','Mg II 4481' ,'S IV 4486','DIB is 4502','S IV 4504','N III 4511/15','He II 4542','Si III 4553/68/75',\
         'O II 4591/4596' ,'N V 4604/20','Si IV 4631', 'N III 4535/414/42','C III 4647/50/51' ,'O II 4662/76' ,'He II 4686', 'O II 4699/705',\
         'He I 4713','DIB is 4727','DIB is 4762','DIB is 4780','H-Beta 4861','DIB is 4881','N III 4905', 'He I 4921'], \
        'rest_wavelength':[[3995], [4009], [4026], [4041, 4044], [4058], [4068, 4069, 4070], \
                [4070, 4072, 4076], [4089], [4097], [4102], [4116], [4121], [4144], [4187], [4200],\
                 [4254], [4267], [4276, 4285], [4317, 4320], [4340], [4349], [4367], [4379], [4388], \
                [4415, 4417], [4429], [4471], [4481], [4486], [4502], [4504], [4511, 4515], [4542], [4553, 4568, 4575], \
                [4591, 4596], [4604, 4620], [4631], [4634, 4641, 4642], [4647, 4650, 4651], [4662, 4676], [4686], [4699, 4705], \
                [4713], [4727],[4762], [4780], [4861], [4881], [4905], [4921]]}

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, ax1_ylim = [-0.1, 2.0], ax2_ylim = [0.5, 2.0])

    pdb.set_trace()

    pyplot.savefig('/user/bostroem/science/2013_greece/H9_R136b_cr_remove_spec.pdf')

    line_dict = normalize_spec.apply_redshift(line_dict, 278.0E3)
    normalize_spec.mark_spectrum(line_dict, ax2)
    normalize_spec.mark_spectrum(line_dict, ax3)

    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H9_R136b_cr_remove_spec_line_labels.pdf')