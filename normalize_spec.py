import scipy
from scipy import signal
from scipy.signal import medfilt

from matplotlib import pyplot
import numpy as np
import pdb
from astropy import constants

def normalize_spectrum(wl, net, kernel_size = 201):
    fig = pyplot.figure()
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

    interpol_back = np.interp(wl, background_wl, background_net)
    continuum = medfilt(interpol_back, kernel_size = 69)
    ax.plot(wl, continuum, 'c', lw = 2)
    ax.plot(wl, np.polyval(np.polyfit(background_wl, background_net, 4), wl), 'r', lw = 2)
    pyplot.draw()
    return wl, net/continuum, continuum

def apply_redshift(line_dict, v):
    #velocity should be in m/s
    redshift_list = []
    for iline in line_dict['rest_wavelength']:
        redshift_list.append(iline* (1.0 + v/constants.c.value))
    line_dict['redshift_wavelength'] = redshift_list
    return line_dict

def mark_spectrum(line_dict, ax):
    for label, value in zip(line_dict['name'], line_dict['redshift_wavelength']):
        ax.axvline(value)
        ymin, ymax = ax.get_ylim()
        t = ax.text(value - 8, ax.get_ylim()[1], label, horizontalalignment = 'right',  rotation = 'vertical', alpha = 0.3)
        #t.set_fontweight('bold')
        #t.set_backgroundcolor('w')
        #t.set_rotation('vertical')
    return ax
        

line_dict = {'name': ['N IV 3480', 'Si IV 4089', 'Si IV 4116', 'N IV 4058', 'He II 4200', 'He I 4471', 'He II 4542', 'N III 4636', 'N III 4642', 'He II 4686', 'H-Beta', 'H-Gamma', 'H-Delta'], \
            'rest_wavelength': [3480, 4089, 4116, 4058, 4200, 4471, 4542, 4636, 4642, 4686, 4861.33, 4340.47, 4101.74]}

H-Beta 4861.33
H-Gamma 4340.47
H-Delta 4101.74 

