import sys
#sys.path.append('/user/bostroem/science/programs')
sys.path.append('/Users/bostroem/science/programs')
import normalize_spec
import pyfits
from matplotlib import pyplot
from stsci.convolve import boxcar
import pdb
import numpy as np
from optparse import OptionParser
import os

def make_plots(fuv, ccd1, ccd2, ccd3, ccd4, reject_cr_flag, norm_spec_flag, ax1_ylim, ax2_ylim, title, bin_flag = None):
    check_existance_of_files(ccd1, ccd2, ccd3, ccd4, fuv)
    if norm_spec_flag == True:
        tbdata5 = pyfits.getdata(fuv, 1)
        wl5 = tbdata5['wavelength'].ravel()
        net5 = tbdata5['flux'].ravel()    
        #dq = tbdata5['dq'].ravel()
        #wl5 = wl5[dq&512 != 512]
        #net5 = net5[dq&512 != 512]
        #cr_reject
        if reject_cr_flag:
            #Open data files
            tbdata1 = pyfits.getdata(ccd1, 1)
            tbdata2 = pyfits.getdata(ccd2, 1)
            tbdata3 = pyfits.getdata(ccd3, 1)
            tbdata4 = pyfits.getdata(ccd4, 1)
    
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
    
        else:
            tbdata1 = pyfits.getdata(ccd1.replace('.fits', '_crj.fits'), 1)
            wl1 = tbdata1['wavelength'].ravel()
            net1 = tbdata1['net'].ravel()
            tbdata2 = pyfits.getdata(ccd2.replace('.fits', '_crj.fits'), 1)
            wl2 = tbdata2['wavelength'].ravel()
            net2 = tbdata2['net'].ravel()
            tbdata3 = pyfits.getdata(ccd3.replace('.fits', '_crj.fits'), 1)
            wl3 = tbdata3['wavelength'].ravel()
            net3 = tbdata3['net'].ravel()
            tbdata4 = pyfits.getdata(ccd4.replace('.fits', '_crj.fits'), 1)
            wl4 = tbdata4['wavelength'].ravel()
            net4 = tbdata4['net'].ravel()
    
        net1 = normalize_spec.manually_remove_CR(wl1, net1, flt = '/user/bostroem/science/combined_dithered_images/%s_%i_combined_img.fits' %(pyfits.getval(ccd1, 'targname', 0)[5:], pyfits.getval(ccd1, 'cenwave', 0)))
        write_fits_file(ccd1.replace('.fits', '_crj.fits'), wl1, net1, pyfits.getheader(ccd1, 0), spec_name = 'net')
    
        norm_spec2 = normalize_spec.manually_remove_CR(wl2, net2, flt = '/user/bostroem/science/combined_dithered_images/%s_%i_combined_img.fits' %(pyfits.getval(ccd2, 'targname', 0)[5:], pyfits.getval(ccd2, 'cenwave', 0)))
        write_fits_file(ccd2.replace('.fits', '_crj.fits'), wl2, net2, pyfits.getheader(ccd2, 0), spec_name = 'net')
    
        norm_spec3 = normalize_spec.manually_remove_CR(wl3, net3, flt = '/user/bostroem/science/combined_dithered_images/%s_%i_combined_img.fits' %(pyfits.getval(ccd3, 'targname', 0)[5:], pyfits.getval(ccd3, 'cenwave', 0)))
        write_fits_file(ccd3.replace('.fits', '_crj.fits'), wl3, net3, pyfits.getheader(ccd3, 0), spec_name = 'net')
    
        norm_spec4 = normalize_spec.manually_remove_CR(wl4, net4, flt = '/user/bostroem/science/combined_dithered_images/%s_%i_combined_img.fits' %(pyfits.getval(ccd4, 'targname', 0)[5:], pyfits.getval(ccd4, 'cenwave', 0)))
        write_fits_file(ccd4.replace('.fits', '_crj.fits'), wl4, net4, pyfits.getheader(ccd4, 0), spec_name = 'net')
    
    
    
    
        if bin_flag:
            wl1, net1 = bin_data(wl1, net1, bin_size = 3)
            wl2, net2 = bin_data(wl2, net2, bin_size = 3)
            wl3, net3 = bin_data(wl3, net3, bin_size = 3)
            wl4, net4 = bin_data(wl4, net4, bin_size = 3)

        #normalize spectra


        wl1, norm_spec1, cont1 = normalize_spec.normalize_spectrum(wl1, net1)
        write_fits_file(ccd1.replace('.fits', '_norm.fits'), wl1, norm_spec1, pyfits.getheader(ccd1, 0), spec_name = 'net')
    
        wl2, norm_spec2, cont2 = normalize_spec.normalize_spectrum(wl2, net2)
        write_fits_file(ccd2.replace('.fits', '_norm.fits'), wl2, norm_spec2, pyfits.getheader(ccd2, 0), spec_name = 'net')
    
        wl3, norm_spec3, cont3 = normalize_spec.normalize_spectrum(wl3, net3)
        write_fits_file(ccd3.replace('.fits', '_norm.fits'), wl3, norm_spec3, pyfits.getheader(ccd3, 0), spec_name = 'net')
    
        wl4, norm_spec4, cont4 = normalize_spec.normalize_spectrum(wl4, net4)
        write_fits_file(ccd4.replace('.fits', '_norm.fits'), wl4, norm_spec4, pyfits.getheader(ccd4, 0), spec_name = 'net')
    
        wl5, norm_spec5, cont5 = normalize_spec.normalize_spectrum(wl5, net5, kernel_size = 401)
        write_fits_file(fuv .replace('.fits', '_norm.fits'), wl5, norm_spec5, pyfits.getheader(fuv, 0))

    else:
        tbdata1 = pyfits.getdata(ccd1.replace('.fits', '_norm.fits'), 1)
        wl1 = tbdata1['wavelength']
        norm_spec1 = tbdata1['net']
        tbdata2 = pyfits.getdata(ccd2.replace('.fits', '_norm.fits'), 1)
        wl2 = tbdata2['wavelength']
        norm_spec2 = tbdata2['net']
        tbdata3 = pyfits.getdata(ccd3.replace('.fits', '_norm.fits'), 1)
        wl3 = tbdata3['wavelength']
        norm_spec3 = tbdata3['net']
        tbdata4 = pyfits.getdata(ccd4.replace('.fits', '_norm.fits'), 1)
        wl4 = tbdata4['wavelength']
        norm_spec4 = tbdata4['net']
        tbdata5 = pyfits.getdata(fuv.replace('.fits', '_norm.fits'), 1)
        wl5 = tbdata5['wavelength']
        norm_spec5 = tbdata5['flux']


    #create plotting objects
    fig = pyplot.figure(6, figsize = [26, 15])
    fig.clf()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    ax2.plot(wl1, norm_spec1, 'k')
    ax2.plot(wl2, norm_spec2, 'b')
    ax2.plot(wl3, norm_spec3, 'k')
    ax2.plot(wl4, norm_spec4, 'b')

    ax3.plot(wl1, boxcar(norm_spec1, (3,)), 'k')
    ax3.plot(wl2, boxcar(norm_spec2, (3,)), 'b')
    ax3.plot(wl3, boxcar(norm_spec3, (3,)), 'k')
    ax3.plot(wl4, boxcar(norm_spec4, (3,)), 'b')

    ax1.plot(wl5, norm_spec5, 'k')

    ax2.set_xlim(3900, 4750)
    ax3.set_xlim(3900, 4750)
    ax1.set_xlim(1150, 1700)
    ax1.set_ylim(ax1_ylim[0], ax1_ylim[1])
    ax2.set_ylim(ax2_ylim[0], ax2_ylim[1])
    ax3.set_ylim(ax2_ylim[0], ax2_ylim[1])

    ax1_xticks = np.arange(1150, 1750, 50)
    ax2_xticks = np.arange(3900, 4750, 50)
    ax1.set_xticks(ax1_xticks)
    ax2.set_xticks(ax2_xticks)
    ax3.set_xticks(ax2_xticks)
    
    ax1.grid(b = 'on', color = 'gray', ls = ':')
    ax2.grid(b = 'on', color = 'gray', ls = ':')
    ax3.grid(b = 'on', color = 'gray', ls = ':')

    #Label Plots
    ax2.set_title('Optical Spectrum %s' %(title))
    ax1.set_title('FUV Spectrum %s' %(title))
    ax3.set_title('Optical Spectrum %s R~4000' %(title))

    ax3.set_xlabel('Wavelength')
    ax2.set_ylabel('Relative Intensity')
    pyplot.draw()
    return ax1, ax2, ax3, fig


def check_existance_of_files(ccd1, ccd2, ccd3, ccd4, fuv):
    if not os.path.exists(ccd1):
        raise IOError('%s does not exist' %(ccd1))
    if not os.path.exists(ccd2):
        raise IOError('%s does not exist' %(ccd2))
    if not os.path.exists(ccd3):
        raise IOError('%s does not exist' %(ccd3))
    if not os.path.exists(ccd4):
        raise IOError('%s does not exist' %(ccd4))
    if not os.path.exists(fuv):
        raise IOError('%s does not exist' %(fuv))

def write_fits_file(filename, wl, spec, hdr, spec_name = 'flux'):
    c1 = pyfits.Column(name = 'wavelength', format = 'D', unit = 'Ang', array = wl)
    c2 = pyfits.Column(name = spec_name, format = 'D', array = spec)
    tbhdu = pyfits.new_table([c1, c2])
    hdu = pyfits.PrimaryHDU()
    thdulist = pyfits.HDUList([hdu, tbhdu])
    thdulist.writeto(filename, clobber = True)
    
def label_spectrum(ax1, ax2, ax3):

    NII_dict = {'name':['N II 3995', 'N II 4041/44'],\
                'rest_wavelength':[[3995],  [4041, 4044]]}

    NIII_dict = {'name':['N III 4097','N III 4379' ,'N III 4511/15','N III 4634/14/42','N III 4905'],\
                'rest_wavelength':[[4097],[4379],[4511, 4515],[4634, 4641, 4642],[4905]]}
    NIV_dict = {'name':['N IV 4058'], \
                'rest_wavelength': [[4058]]}

    NV_dict = {'name':['N V 1238/42','N V 4604/20'],\
                'rest_wavelength':[[1238, 1242],[4604, 4620]]}

    HeI_dict = {'name':['He I 4009' ,'He I 4121','He I 4144','He I 4388', 'He I 4471','He I 4713','He I 4921'],\
                'rest_wavelength':[[4009], [4121],[4144],  [4388],[4471],[4713],[4921]]}

    HeI_II_dict = {'name':['He I+II 4026'],\
                'rest_wavelength':[[4026]]}

    HeII_dict = {'name':['He II 1640', 'He II 4200','He II 4542','He II 4686'],\
                'rest_wavelength':[ [1640], [4200],[4542], [4686]]}

    NiIV_dict = {'name':['Ni IV 4058'],\
                'rest_wavelength':[[4058]]} 

    CIV_dict = {'name':['C IV 1548/51'], \
                'rest_wavelength':[[1548, 1551]]}

    CIII_dict = {'name':['C III 1176','C III 4068/69/70','C III 4187','C III 4647/50/51'],\
                'rest_wavelength':[ [1176],[4068, 4069, 4070],[4187],[4647, 4650, 4651]]}

    CII_dict = {'name':['C II 4267'],\
                'rest_wavelength':[[4267]]}

    OII_dict = {'name':['O II 4070/72/76','O II 4254','O II 4276/85','O II 4317/20', 'O II 4349','O II 4367','O II 4415/17','O II 4591/4596' ,'O II 4662/76' ,'O II 4699/705'],\
                'rest_wavelength':[[4070, 4072, 4076],[4254],[4276, 4285], [4317, 4320],[4349], [4367], [4415, 4417],[4591, 4596],[4662, 4676],[4699, 4705]]}

    OV_dict = {'name':['O V 1371'],\
                'rest_wavelength':[[1371]]}

    SiIV_dict = {'name':['Si IV 1393/402','Si IV 4089','Si IV 4116'],\
                'rest_wavelength':[[1391, 1402], [4089], [4116]]}

    SiIII_dict = {'name':['Si III 4553/68/75'],\
                'rest_wavelength':[[4553, 4568, 4575]]}
    SV_dict = {'name': ['SV 1502'], \
                'rest_wavelength':[[1502]]}

    MgII_dict = {'name':['Mg II 4481'],\
                'rest_wavelength':[[4481]]}

    SIV_dict = {'name':['S IV 4486','S IV 4504'],\
                'rest_wavelength':[[4486],[4504]]}
    
    balmer_dict = {'name':[r'H $\epsilon$ 3970', r'H $\delta$ 4102',r'H $\gamma$ 4340',r'H $\beta$ 4861'],\
                'rest_wavelength':[[3970], [4102], [4340],[4861]]}

    dib_dict = {'name':['DIB is 4429','DIB is 4502','DIB is 4727','DIB is 4762','DIB is 4780','DIB is 4881'],\
                'rest_wavelength':[[4429],  [4502], [4727],[4762], [4780],[4881]]}

    ISM_dict = {'name': ['Ly $\\alpha$ 1215', 'Si II 1260', 'OI 1302',  'C II 1335', 'Si II 1526', 'Fe II 1608', 'Al II 1671', 'Ca II H 3968', 'Ca II K 3933'], \
                'rest_wavelength': [[1215.6], [1260.4], [1302.2], [1334.5], [1526.7], [1608.5], [1670.8], [3968], [3933]]}

    ISM_dict_shifted = {'name': ['Si II 1304'], \
                'rest_wavelength': [[1304.4]]}


    line_dict = normalize_spec.apply_redshift(HeI_dict, 268.0E3)
    normalize_spec.mark_spectrum(HeI_dict, ax1, 'r')
    normalize_spec.mark_spectrum(HeI_dict, ax2, 'r')
    normalize_spec.mark_spectrum(HeI_dict, ax3, 'r')

    line_dict = normalize_spec.apply_redshift(HeII_dict, 268.0E3)
    normalize_spec.mark_spectrum(HeII_dict, ax1, 'r', offset = 10.0)
    normalize_spec.mark_spectrum(HeII_dict, ax2, 'r', offset = 10.0)
    normalize_spec.mark_spectrum(HeII_dict, ax3, 'r', offset = 10.0)

    line_dict = normalize_spec.apply_redshift(HeI_II_dict, 268.0E3)
    normalize_spec.mark_spectrum(HeI_II_dict, ax1, 'r')
    normalize_spec.mark_spectrum(HeI_II_dict, ax2, 'r')
    normalize_spec.mark_spectrum(HeI_II_dict, ax3, 'r')

    line_dict = normalize_spec.apply_redshift(NIII_dict, 268.0E3)
    normalize_spec.mark_spectrum(NIII_dict, ax1, 'purple')
    normalize_spec.mark_spectrum(NIII_dict, ax2, 'purple')
    normalize_spec.mark_spectrum(NIII_dict, ax3, 'purple')


    line_dict = normalize_spec.apply_redshift(NV_dict, 268.0E3)
    normalize_spec.mark_spectrum(NV_dict, ax1, 'purple')
    normalize_spec.mark_spectrum(NV_dict, ax2, 'purple')
    normalize_spec.mark_spectrum(NV_dict, ax3, 'purple')

    line_dict = normalize_spec.apply_redshift(NIV_dict, 268.0E3)
    normalize_spec.mark_spectrum(NIV_dict, ax1, 'purple')
    normalize_spec.mark_spectrum(NIV_dict, ax2, 'purple')
    normalize_spec.mark_spectrum(NIV_dict, ax3, 'purple')

    line_dict = normalize_spec.apply_redshift(CIII_dict, 268.0E3)
    normalize_spec.mark_spectrum(CIII_dict, ax1, 'g', offset = 10.0)
    normalize_spec.mark_spectrum(CIII_dict, ax2, 'g', offset = 10.0)
    normalize_spec.mark_spectrum(CIII_dict, ax3, 'g', offset = 10.0)

    line_dict = normalize_spec.apply_redshift(SiIV_dict, 268.0E3)
    normalize_spec.mark_spectrum(SiIV_dict, ax1, 'c')
    normalize_spec.mark_spectrum(SiIV_dict, ax2, 'c')
    normalize_spec.mark_spectrum(SiIV_dict, ax3, 'c')

    line_dict = normalize_spec.apply_redshift(balmer_dict, 268.0E3)
    normalize_spec.mark_spectrum(balmer_dict, ax1, 'brown', offset = 9.0)
    normalize_spec.mark_spectrum(balmer_dict, ax2, 'brown', offset = 9.0)
    normalize_spec.mark_spectrum(balmer_dict, ax3, 'brown', offset = 9.0)

    line_dict = normalize_spec.apply_redshift(CIV_dict, 268.0E3)
    normalize_spec.mark_spectrum(CIV_dict, ax1, 'g')
    normalize_spec.mark_spectrum(CIV_dict, ax2, 'g')
    normalize_spec.mark_spectrum(CIV_dict, ax3, 'g')

    line_dict = normalize_spec.apply_redshift(OV_dict, 268.0E3)
    normalize_spec.mark_spectrum(OV_dict, ax1, '#FF66CC')
    normalize_spec.mark_spectrum(OV_dict, ax2, '#FF66CC')
    normalize_spec.mark_spectrum(OV_dict, ax3, '#FF66CC')

    line_dict = normalize_spec.apply_redshift(ISM_dict, 268.0E3)
    normalize_spec.mark_spectrum(ISM_dict, ax1, 'gray')
    normalize_spec.mark_spectrum(ISM_dict, ax2, 'gray')
    normalize_spec.mark_spectrum(ISM_dict, ax3, 'gray')

    line_dict = normalize_spec.apply_redshift(ISM_dict_shifted, 268.0E3)
    normalize_spec.mark_spectrum(ISM_dict_shifted, ax1, 'gray', offset = 5.0)
    normalize_spec.mark_spectrum(ISM_dict_shifted, ax2, 'gray', offset = 5.0)
    normalize_spec.mark_spectrum(ISM_dict_shifted, ax3, 'gray', offset = 5.0)

    line_dict = normalize_spec.apply_redshift(SV_dict, 268.0E3)
    normalize_spec.mark_spectrum(SV_dict, ax1, 'm')
    normalize_spec.mark_spectrum(SV_dict, ax2, 'm')
    normalize_spec.mark_spectrum(SV_dict, ax3, 'm')

    return ax1, ax2, ax3

def bin_data(wl, net, bin_size = 3):
    array_len = 3*len(wl[bin_size - 1::bin_size])  #In case the len(wl) is not evenly divisible by bin_size
    bin_wl = wl[0:array_len:bin_size]
    bin_net = net[0:array_len:bin_size]
    for i in range(1, bin_size):
        bin_wl += wl[i:array_len:bin_size]
        bin_net += net[i:array_len:bin_size]
        #pdb.set_trace()
    bin_wl = bin_wl / float(bin_size)
    return bin_wl, bin_net

def make_r136b_spec(options):
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_3936_combined_img_loc545.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4194_combined_img_loc496.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4451_combined_img_loc546.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4706_combined_img_loc543.fits'
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0005_G140L.fits'
    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 2.3], [0.6, 2.1], 'H9, R136b')
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H9_R136b_cr_remove_spec.pdf')
    pyplot.close()
    pyplot.close()
    pyplot.close()
    pyplot.close()
    pyplot.close()
    pyplot.close()

def make_H36_spec(options):
    #ccd1 = '/Users/bostroem/dropbox/R136/ccd/H36_3936SE3loc547.fits.gz'
    #ccd2 = '/Users/bostroem/dropbox/R136/ccd/H36_4194SE3loc495.fits.gz'
    #ccd3 = '/Users/bostroem/dropbox/R136/ccd/H36_4451SE3loc548.fits.gz'
    #ccd4 = '/Users/bostroem/dropbox/R136/ccd/H36_4706SE3loc544.fits.gz'
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0015_G140L.fits'
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_3936_combined_img_loc546.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4194_combined_img_loc494.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4451_combined_img_loc547.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4706_combined_img_loc544.fits'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 2.3], [0.6, 1.4], 'H36', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H36_cr_remove_spec.pdf')

def make_r136a7_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0007_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/R136a7_3936NW3loc517.fits.gz'
    ccd2 = '/Users/bostroem/Dropbox/R136/ccd/R136a7_4194NW3loc526.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/R136a7_4451NW3loc517.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/R136a7_4706NW3loc515.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 2.3], [0.6, 1.4], 'R136a7', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/r136a7_cr_remove_spec.pdf')

def make_r136a6_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0004_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/R136a6_3936SE1loc506.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4194SE1loc535.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4451SE1loc507.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4706SE1loc504.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'R136a6', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/R136a6_cr_remove_spec.pdf')

def make_H50_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0015_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H50_3936SE3loc528.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H50_4194SE3loc514.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H50_4451SE3loc529.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H50_4706SE3loc525.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H50', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H50_cr_remove_spec.pdf') 

def make_H40_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0020_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H40_3936NW8loc500.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H40_4194NW8loc542.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H40_4451NW8loc500.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H40_4706NW8loc499.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H40', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H40_cr_remove_spec.pdf') 

def make_H55_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0023_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H55_3936SE9loc516.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H55_4194SE9loc525.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H55_4451SE9loc517.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H55_4706SE9loc514.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H55', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H55_cr_remove_spec.pdf')

def make_H62_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0024_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H62_3936SE3loc513.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H62_4194SE3loc533.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H62_4451SE3loc514.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H62_4706SE3loc510.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H62', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H62_cr_remove_spec.pdf')

def make_H58_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0024_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H58_3936NW1loc531.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H58_4194NW1loc510.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H58_4451NW1loc532.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H58_4706NW1loc529.fits'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H58', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H58_cr_remove_spec.pdf') 


def make_H48_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0026_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H48_3936SE2loc542.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H48_4194SE2loc500.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H48_4451SE2loc542.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H48_4706SE2loc539.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H48', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H48_cr_remove_spec.pdf')


def make_H47_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0027_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H47_3936NW3loc489.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H47_4194NW3loc554.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H47_4451NW3loc489.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H47_4706NW3loc487.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H47', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H47_cr_remove_spec.pdf') 

def make_H35_spec(options):
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0012_G140L.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H35_3936NW5loc528.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H35_4194NW5loc515.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H35_4451NW5loc528.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H35_4706NW5loc526.fits.gz'

    ax1, ax2, ax3, fig = make_plots(fuv, ccd1, ccd2, ccd3, ccd4, options.reject_cr_flag, options.norm_spec_flag, [-0.1, 1.8], [0.6, 1.4], 'H35', bin_flag = options.bin_flag)
    ax1, ax2, ax3 = label_spectrum(ax1, ax2, ax3)
    pyplot.draw()
    pdb.set_trace()
    pyplot.savefig('/user/bostroem/science/2013_greece/H35_cr_remove_spec.pdf') 

if __name__ == "__main__":

    #If --norm is not supllied norm_spec_flag is set to True = user wants to normalize the spectrum
    parser = OptionParser()
    parser.add_option('--crrej', dest = 'reject_cr_flag', action = 'store_false', help = 'Use previously produce CR rejected file', default = True)
    parser.add_option('--norm', dest = 'norm_spec_flag', action = 'store_false', help = 'Use previously normalized spectrum file', default = True)    
    parser.add_option('--bin', dest = 'bin_flag', action = 'store_true', help = 'Bin data', default = False)
    (options, args) = parser.parse_args()

    ############## H9/R136 b################
    #make_r136b_spec(options)
    ############## H36################
    #make_H36_spec(options)
    ############## R136/a6################
    #make_r136a6_spec(options)
    ############## H62 ################
    #make_H62_spec(options)
    ############## H50 ################
    #make_H50_spec(options)
    ############## H48 ################
    #make_H48_spec(options)
    ############## H55 ################
    #make_H55_spec(options)
    ############## H47 ################
    #make_H47_spec(options)
    ############## H40 ################
    #make_H40_spec(options)
    ############## H35 ################
    #make_H35_spec(options)
    ############## R136/a7################
    #make_r136a7_spec(options)
    ############## H58 ################
    make_H58_spec(options)




    
