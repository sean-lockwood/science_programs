import scipy
from scipy import signal
from scipy.signal import medfilt

from matplotlib import pyplot
import numpy as np
import pdb
from astropy import constants
import numpy as np

def normalize_spectrum(wl, net, kernel_size = 201):
    fig = pyplot.figure(figsize = [22, 7])
    ax = fig.add_subplot(1,1,1)
    wl = wl.flatten()
    net = net.flatten()
    ax.plot(wl, net, 'b')
    background_wl = np.empty((0,))
    background_net = np.empty((0,))
    raw_input('Zoom in on spectrum and press enter to continue ')
    select_another_region = 'y'
    print 'Select background regions '
    while select_another_region != 'n':
        pt1, pt2 = pyplot.ginput(n = 2, timeout = -1)
        x1, y1 = pt1
        x2, y2 = pt2
        background_indx = np.where((wl < max(x1, x2)) & (wl > min(x1, x2)))
        background_wl = np.append(background_wl, wl[background_indx])
        background_net = np.append(background_net, net[background_indx])
        pyplot.plot([wl[background_indx[0][0]], wl[background_indx[0][-1]]], [net[background_indx[0][0]], net[background_indx[0][-1]]], 'm--|', lw = 3)

        select_another_region = raw_input('Would you like to select another background regions? (y), n ')
    background_wl = background_wl.flatten()
    background_net = background_net.flatten()
    sort_indx = np.argsort(background_wl)
    background_wl = background_wl[sort_indx]
    background_net = background_net[sort_indx]
    poly_coeff = np.polyfit(background_wl, background_net, 3)
    continuum = np.polyval(poly_coeff, wl)
    ax.plot(wl, continuum, 'c', lw = 2)
    pyplot.draw()
    return wl, net/continuum, continuum

def apply_redshift(line_dict, v):
    #velocity should be in m/s
    redshift_list = []
    for ilist in line_dict['rest_wavelength']:
        elem_list = []
        for iline in ilist:
            elem_list.append(iline* (1.0 + v/constants.c.value))
        redshift_list.append(elem_list)
    line_dict['redshift_wavelength'] = redshift_list
    return line_dict

def mark_spectrum(line_dict, ax):
    for label, lines in zip(line_dict['name'], line_dict['redshift_wavelength']):
        for value in lines:
            ax.axvline(value, color = 'k', linestyle = '--')
        ymin, ymax = ax.get_ylim()
        t = ax.text(lines[0], ymax - (ymax - ymin)/8.0, label, horizontalalignment = 'right',  rotation = 'vertical', alpha = 0.5)
        #t.set_fontweight('bold')
        #t.set_backgroundcolor('w')
        #t.set_rotation('vertical')
    return ax
        
def manually_remove_CR(wl, spec, flt = None):
    fig = pyplot.figure(figsize = [22, 14])
    ax = fig.add_subplot(2,1,1)
    if flt:
        ax2 = fig.add_subplot(2, 1, 2)
        img = pyfits.getdata(flt, 1)
        ax2.imshow(img, interpolation = 'nearest')
    ax.plot(wl, spec)
    pyplot.draw()
    remove_CR_flag = raw_input('Would you like to remove any cosmic-rays? (y), n ')
    while remove_CR_flag != 'n':
        raw_input('Zoom in on cosmic ray and press ENTER when finished ')
        print 'Click on cosmic ray on plot'
        x, y = pyplot.ginput(n = 1, timeout = -1)[0]
        rm_indx = np.where(abs(wl - x) == np.min(abs(wl - x)))[0]
        spec[rm_indx] = (spec[rm_indx - 1] + spec[rm_indx + 1])/2.0
        ax.cla()
        ax.plot(wl, spec)
        pyplot.draw()
        remove_CR_flag = raw_input('Would you like to remove any cosmic-rays? (y), n ')
    pyplot.close()
    return spec





line_dict = {'name': ['N IV 3480', 'Si IV 4089', 'Si IV 4116', 'N IV 4058', 'He II 4200', 'He I 4471', 'He II 4542', 'N III 4636', 'N III 4642', 'He II 4686', 'H-Beta', 'H-Gamma', 'H-Delta'], \
            'rest_wavelength': [[3480], [4089], [4116], [4058], [4200], [4471], [4542], [4636], [4642], [4686], [4861.33], [4340.47], [4101.74]]}




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
