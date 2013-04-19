import scipy
from scipy import signal
from scipy.signal import medfilt

from matplotlib import pyplot
import numpy as np
import pdb

def normalize_spectrum(wl, net, kernel_size = 201):
    wl = wl.flatten()
    net = net.flatten()
    pyplot.plot(wl, net, 'b')
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
    pyplot.plot(wl, continuum, 'c')
    return wl, net - continuum, continuum


        