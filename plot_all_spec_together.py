import pyfits
from matplotlib import pyplot
import pdb
from degrade import degrader
import normalize_spec
import numpy as np

def plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, offset, ax_ccd, ax_fuv):

    tbdata1 = pyfits.getdata(ccd1, 1)
    print tbdata1.names
    wl1 = tbdata1['wavelength'].ravel()
    net1 = tbdata1['net'].ravel()
    net1_4000 = degrader(wl1, net1, 7000., 4000., quick = 5)
    tbdata2 = pyfits.getdata(ccd2, 1)
    wl2 = tbdata2['wavelength'].ravel()
    net2 = tbdata2['net'].ravel()
    net2_4000 = degrader(wl2, net2,7000., 4000., quick = 5)
    tbdata3 = pyfits.getdata(ccd3, 1)
    wl3 = tbdata3['wavelength'].ravel()
    net3 = tbdata3['net'].ravel()
    net3_4000 = degrader(wl3, net3, 7000., 4000., quick = 5)
    tbdata4 = pyfits.getdata(ccd4, 1)
    wl4 = tbdata4['wavelength'].ravel()
    net4 = tbdata4['net'].ravel()
    net4_4000 = degrader(wl4, net4, 7000., 4000., quick = 5)
    tbdata5 = pyfits.getdata(fuv, 1)
    wl5 = tbdata5['wavelength'].ravel()
    net5 = tbdata5['flux'].ravel()
    
    ax_fuv.plot(wl5, net5+offset, c = 'k')
    ax_fuv.set_xlim(1150, 1700)
    ax_fuv.set_ylim(-0.1, 2.3+offset)

    ax_ccd.plot(wl1, net1_4000+offset, c = 'k')
    ax_ccd.plot(wl2, net2_4000+offset, c = 'b')
    ax_ccd.plot(wl3, net3_4000+offset, c = 'k')
    ax_ccd.plot(wl4, net4_4000+offset, c = 'b')
    ax_ccd.set_xlim(3900, 4750)
    ax_ccd.set_ylim(-0.1, 2.3+offset)
    return ax_ccd, ax_fuv

def label_spectrum(ax_fuv, ax_ccd):

    NII_dict = {'name':['N II 3995', 'N II 4041/44'],\
                'rest_wavelength':[[3995],  [4041, 4044]]}

    NIII_dict = {'name':['N III 4097','N III 4379' ,'N III 4511/15','N III 4634/40/42','N III 4905'],\
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
    normalize_spec.mark_spectrum(HeI_dict, ax_fuv, 'r', offset = 8)
    normalize_spec.mark_spectrum(HeI_dict, ax_ccd, 'r', offset = 8)


    line_dict = normalize_spec.apply_redshift(HeII_dict, 268.0E3)
    normalize_spec.mark_spectrum(HeII_dict, ax_fuv, 'r', offset = 11)
    normalize_spec.mark_spectrum(HeII_dict, ax_ccd, 'r', offset = 11)

    line_dict = normalize_spec.apply_redshift(HeI_II_dict, 268.0E3)
    normalize_spec.mark_spectrum(HeI_II_dict, ax_fuv, 'r')
    normalize_spec.mark_spectrum(HeI_II_dict, ax_ccd, 'r')

    line_dict = normalize_spec.apply_redshift(NIII_dict, 268.0E3)
    normalize_spec.mark_spectrum(NIII_dict, ax_fuv, 'purple')
    normalize_spec.mark_spectrum(NIII_dict, ax_ccd, 'purple')

    line_dict = normalize_spec.apply_redshift(NV_dict, 268.0E3)
    normalize_spec.mark_spectrum(NV_dict, ax_fuv, 'purple')
    normalize_spec.mark_spectrum(NV_dict, ax_ccd, 'purple')

    line_dict = normalize_spec.apply_redshift(NIV_dict, 268.0E3)
    normalize_spec.mark_spectrum(NIV_dict, ax_fuv, 'purple')
    normalize_spec.mark_spectrum(NIV_dict, ax_ccd, 'purple')


    line_dict = normalize_spec.apply_redshift(CIII_dict, 268.0E3)
    normalize_spec.mark_spectrum(CIII_dict, ax_fuv, 'g', offset = 12)
    normalize_spec.mark_spectrum(CIII_dict, ax_ccd, 'g', offset = 12)

    line_dict = normalize_spec.apply_redshift(SiIV_dict, 268.0E3)
    normalize_spec.mark_spectrum(SiIV_dict, ax_fuv, 'c')
    normalize_spec.mark_spectrum(SiIV_dict, ax_ccd, 'c')

    line_dict = normalize_spec.apply_redshift(balmer_dict, 268.0E3)
    normalize_spec.mark_spectrum(balmer_dict, ax_fuv, 'brown', offset = 9.0)
    normalize_spec.mark_spectrum(balmer_dict, ax_ccd, 'brown', offset = 9.0)

    line_dict = normalize_spec.apply_redshift(CIV_dict, 268.0E3)
    normalize_spec.mark_spectrum(CIV_dict, ax_fuv, 'g')
    normalize_spec.mark_spectrum(CIV_dict, ax_ccd, 'g')

    line_dict = normalize_spec.apply_redshift(OV_dict, 268.0E3)
    normalize_spec.mark_spectrum(OV_dict, ax_fuv, '#FF66CC')
    normalize_spec.mark_spectrum(OV_dict, ax_ccd, '#FF66CC')

    line_dict = normalize_spec.apply_redshift(ISM_dict, 268.0E3)
    normalize_spec.mark_spectrum(ISM_dict, ax_fuv, 'gray')
    normalize_spec.mark_spectrum(ISM_dict, ax_ccd, 'gray')

    line_dict = normalize_spec.apply_redshift(ISM_dict_shifted, 268.0E3)
    normalize_spec.mark_spectrum(ISM_dict_shifted, ax_fuv, 'gray', offset = 15.0)
    normalize_spec.mark_spectrum(ISM_dict_shifted, ax_ccd, 'gray', offset = 15.0)

    line_dict = normalize_spec.apply_redshift(SV_dict, 268.0E3)
    normalize_spec.mark_spectrum(SV_dict, ax_fuv, 'm')
    normalize_spec.mark_spectrum(SV_dict, ax_ccd, 'm')

    return ax_fuv, ax_ccd


if __name__ == "__main__":
    fig = pyplot.figure(1)
    ax_fuv = fig.add_subplot(1,2,1)
    ax_ccd = fig.add_subplot(1,2,2)
    ax_fuv.set_position([0.05, 0.1, 0.35, 0.8])
    ax_ccd.set_position([0.425, 0.1, 0.55, 0.8])
    ax_fuv.set_title('FUV Spectra')
    ax_ccd.set_title('Optical Spectra')

    loc_list = np.arange(14)*0.75 - 0.5

    #H58
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0024_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H58_3936NW1loc531_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H58_4194NW1loc510_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H58_4451NW1loc532_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H58_4706NW1loc529_norm.fits'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[0], ax_ccd, ax_fuv)


    #r136a7
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0007_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/R136a7_3936NW3loc517_norm.fits.gz'
    ccd2 = '/Users/bostroem/Dropbox/R136/ccd/R136a7_4194NW3loc526_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/R136a7_4451NW3loc517_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/R136a7_4706NW3loc515_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[1], ax_ccd, ax_fuv)


    #H35
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0012_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H35_3936NW5loc528_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H35_4194NW5loc515_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H35_4451NW5loc528_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H35_4706NW5loc526_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[2], ax_ccd, ax_fuv)


    #H40
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0020_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H40_3936NW8loc500_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H40_4194NW8loc542_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H40_4451NW8loc500_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H40_4706NW8loc499_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[3], ax_ccd, ax_fuv)


    #H47
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0027_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H47_3936NW3loc489_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H47_4194NW3loc554_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H47_4451NW3loc489_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H47_4706NW3loc487_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[4], ax_ccd, ax_fuv)


    #H55
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0023_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H55_3936SE9loc516_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H55_4194SE9loc525_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H55_4451SE9loc517_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H55_4706SE9loc514_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[5], ax_ccd, ax_fuv)


    #H48
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0026_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H48_3936SE2loc542_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H48_4194SE2loc500_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H48_4451SE2loc542_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H48_4706SE2loc539_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[6], ax_ccd, ax_fuv)



    #H50
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0015_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H50_3936SE3loc528_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H50_4194SE3loc514_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H50_4451SE3loc529_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H50_4706SE3loc525_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[7], ax_ccd, ax_fuv)



    #H62
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0024_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H62_3936SE3loc513_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H62_4194SE3loc533_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H62_4451SE3loc514_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H62_4706SE3loc510_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[8], ax_ccd, ax_fuv)

    #r136a6
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0004_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/R136a6_3936SE1loc506_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4194SE1loc535_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4451SE1loc507_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4706SE1loc504_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[9], ax_ccd, ax_fuv)




    #H36
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0015_G140L_norm.fits'
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_3936_combined_img_loc546_norm.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4194_combined_img_loc494_norm.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4451_combined_img_loc547_norm.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4706_combined_img_loc544_norm.fits'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[10], ax_ccd, ax_fuv)


    #r136b
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_3936_combined_img_loc545_norm.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4194_combined_img_loc496_norm.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4451_combined_img_loc546_norm.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4706_combined_img_loc543_norm.fits'
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0005_G140L_norm.fits'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, loc_list[11], ax_ccd, ax_fuv)

    ax_fuv.set_yticks(loc_list - 0.5)
    ax_ccd.set_yticks(loc_list - 0.5)
    ax_ccd.set_ylim(-0.5, np.max(loc_list) + 2.0)
    ax_fuv.set_ylim(-0.50, np.max(loc_list)+ 2.0)
    yticklabels = ax_fuv.get_yticklabels()
    ccd_yticklabels = ax_ccd.get_yticklabels()
    for i in range(len(ccd_yticklabels)):
        ccd_yticklabels[i] = ''
    ax_ccd.set_yticklabels(ccd_yticklabels)

    label_list = ['H58', 'R136a7', 'H35', 'H40', 'H47', 'H55', 'H48', 'H50', 'H62', 'R136a6', 'H36', 'R136b']
    for i in range(len(yticklabels)):
        yticklabels[i] = ''
    for i, star_name in enumerate(label_list):
        yticklabels[i+2] = star_name
        ccd_yticklabels[i+2] = star_name
    ax_fuv.set_yticklabels(yticklabels)
    ax_ccd.minorticks_on()
    ax_fuv.minorticks_on()
    ax_fuv.set_xlabel('Wavelength ($\AA$)')
    ax_ccd.set_xlabel('Wavelength ($\AA$)')
    ax_fuv.set_ylabel('Relative Intensity')
    ax_fuv.tick_params(axis = 'y', which = 'minor', left = 'off')
    ax_ccd.tick_params(axis = 'y', which = 'minor', left = 'off')
    ax_fuv.tick_params(axis = 'y', which = 'minor', right = 'off')
    ax_ccd.tick_params(axis = 'y', which = 'minor', right = 'off')

    ax_fuv, ax_ccd = label_spectrum(ax_fuv, ax_ccd)
    pdb.set_trace()



























    '''


    #H62
    fuv = '/user/bostroem/science/multispec/R136_G140L/star0024_G140L_norm.fits'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/H62_3936SE3_3936_loc513_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/H62_4194SE3loc533_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/H62_4451SE3loc514_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/H62_4706SE3loc510_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, 1.0, ax_ccd, ax_fuv)

    #R136a6
    fuv = '//user/bostroem/science/multispec/R136_G140L/star0004_G140L_norm.fits'
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_3936_combined_img_loc546_norm.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4194_combined_img_loc494_norm.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4451_combined_img_loc547_norm.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE3_4706_combined_img_loc544_norm.fits'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, 1.0, ax_ccd, ax_fuv)

    #H36
    fuv = '/Users/bostroem/dropbox/R136/mama/R136a6_SE1loc378_norm.fits.gz'
    ccd1 = '/Users/bostroem/dropbox/R136/ccd/R136a6_3936SE1loc506_norm.fits.gz'
    ccd2 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4194SE1loc535_norm.fits.gz'
    ccd3 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4451SE1loc507_norm.fits.gz'
    ccd4 = '/Users/bostroem/dropbox/R136/ccd/R136a6_4706SE1loc504_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, 1.5, ax_ccd, ax_fuv)

    #R136b
    ccd1 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_3936_combined_img_loc545_norm.fits'
    ccd2 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4194_combined_img_loc496_norm.fits'
    ccd3 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4451_combined_img_loc546_norm.fits'
    ccd4 = '/user/bostroem/science/12465_otfr20121109/ccd/SE8_4706_combined_img_loc543_norm.fits'
    fuv = '/Users/bostroem/dropbox/R136/mama/R136b_SE8loc459_norm.fits.gz'
    ax_ccd, ax_fuv = plot_spec(ccd1, ccd2, ccd3, ccd4, fuv, 0, ax_ccd, ax_fuv)

    '''