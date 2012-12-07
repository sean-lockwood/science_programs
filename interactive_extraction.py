import pyfits
import numpy as np
import sys
import math
from matplotlib import pyplot
import pdb
import os
import shutil

import pyraf
from pyraf import iraf
from iraf import stsdas,hst_calib,stis,x1d


def collapse_spectrum(img):
    collapsed_img = np.sum(img, axis = 1)/1024
    return collapsed_img

def get_background_locations(fig1, ax1, pix_num, collapsed_img, filename, cenwave):
    background_collapsed = np.empty((0,))
    background_pix = np.empty((0,))
    finished_flag = 'n'
    if os.path.exists(filename.replace('.fits', '_c%i_background_regions.txt' %(cenwave))):
        keep_file_flag = raw_input('Background file already exists: Overwrite file (w) or Append to file (a)')
    else:
        keep_file_flag = 'w'
    ofile = open(filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), keep_file_flag)
    while finished_flag != 'y':
        print 'Select background points for cross-dispersion background subtraction'
        temp_region = fig1.ginput(n = 2, timeout = -1)
        back_reg = [temp_region[0][0], temp_region[1][0]]
        l1 = pyplot.plot(back_reg, np.interp(np.array(back_reg), pix_num, collapsed_img), 'r|-', markersize = 5)
        pyplot.draw()
        keep_flag = raw_input('Keep points? (y), n ')
        if keep_flag != 'n':
            ofile.write('%i\t %i \n' %(int(math.floor(temp_region[0][0])),int(math.ceil(temp_region[1][0]))))
            background_collapsed = np.append(background_collapsed, collapsed_img[int(math.floor(temp_region[0][0])):int(math.ceil(temp_region[1][0]))])
            background_pix = np.append(background_pix, pix_num[int(math.floor(temp_region[0][0])):int(math.ceil(temp_region[1][0]))])
        else:
            l1[0].set_visible(False)
            pyplot.draw()
        finished_flag = raw_input('Finished entering points? (n), y ')
    ofile.close()
    sort_indx = np.argsort(background_pix)
    background_pix = background_pix[sort_indx]
    background_collapsed = background_collapsed[sort_indx]

    return fig1, ax1, background_pix, background_collapsed

def fit_background(ax1, background_pix, background_collapsed, pix_num):
    deg = 5
    change_deg_flag = 'y'
    leg_text = ['%i' %(deg)]
    leg_lines = []
    while change_deg_flag != 'n':
        coeff = np.polyfit(background_pix, background_collapsed, deg)
        fit = np.polyval(coeff, pix_num)
        l1 = ax1.plot(pix_num, fit, lw = 2)
        leg_lines.append(l1[0])
        pyplot.legend(leg_lines, leg_text, loc = 'best')
        pyplot.draw()
        change_deg_flag = raw_input('Replot with a different degree polynomial? (y), n ')
        if change_deg_flag != 'n':
            
            l1[0].set_linestyle('--')
            #l1[0].set_visible(False)
            pyplot.draw()

            try:
                deg = int(raw_input('Enter degree of polynomial you wish to fit: '))
            except ValueError:
                print 'The value you entered cannot be converted to an integer, try again'
                deg = int(raw_input('Enter degree of polynomial you wish to fit: '))

            leg_text.append('%i' %(deg))
    return ax1, fit

def subtract_2D_cross_disp_background(img, fit):
    two_D_background = np.rot90(np.rot90(np.rot90(np.tile(fit, (1024, 1)))))
    return img - two_D_background


def write_temp_subtracted_file(subtracted_img, filename):
    shutil.copyfile(filename, filename.replace('.fits', 'sub.fits'))
    ofile = pyfits.open(filename.replace('.fits', 'sub.fits'), mode = 'update')
    ofile[1].data = subtracted_img
    ofile.flush()
    ofile.close()

def select_extraction_location(fig1, c):
    raw_input('Zoom in on spectrum, hit Enter when finished')
    print 'Select spectrum location'
    keep_flag = 'n'
    while keep_flag == 'n':
        coords = fig1.ginput(n = 1, timeout = -1)
        ln = pyplot.plot([coords[0][0]], [coords[0][1]], '%s+' %(c))
        keep_flag = raw_input('Keep selection? (y), n ')
        if keep_flag == 'n':
            ln[0].set_visible(False)
            pyplot.draw()

    return coords[0][0]
    

def select_extraction_box_size(fig1, extrlocy, c):
    print 'Select extraction box'
    keep_flag = 'n'
    while keep_flag == 'n':
        coords = fig1.ginput(n = 2, timeout = -1)
        half_box = min(abs(coords[0][0] - extrlocy), abs(coords[1][0] - extrlocy))
        lb = pyplot.plot([extrlocy - half_box, extrlocy + half_box],[coords[0][1], coords[0][1]] , '%s|-' %(c), markersize = 10)
        pyplot.draw()
        keep_flag = raw_input('Keep points? (y), n ')
        if keep_flag == 'n':
            lb[0].set_visible(False)
            pyplot.draw()
    return half_box*2

def extract_spectrum(extrlocy, extract_box_size, background_loc1, background_size1, background_loc2, background_size2, filename, c):
    cenwave = pyfits.getval(filename, 'cenwave', 0)
    if os.path.exists(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i_%i.fits' %(int(extrlocy), cenwave)))):
        os.remove(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i_%i.fits' %(int(extrlocy), cenwave))))
        
    iraf.stsdas.hst_calib.stis.x1d(filename.replace('.fits', 'sub.fits'), output = filename.replace('.fits', '_loc%i_%i.fits' %(int(extrlocy), cenwave)), \
                                       a2center = extrlocy, extrsize = extract_box_size, maxsrch = 0, bk1offst = background_loc1 - extrlocy, \
                                       bk2offst = background_loc2 - extrlocy, bk1size = background_size1, bk2size = background_size2)
 

if __name__ == "__main__":
    colors = ['r', 'g', 'c', 'k', 'm']

    pyplot.ion()
    os.environ['oref'] = '/Users/bostroem/science/oref/'
    filename = sys.argv[1]
    try:
        ext = int(sys.argv[2])
    except:
        ext = 1
    find_background_flag = raw_input('Select the background regions? y, (n)')
    fig1 = pyplot.figure(1)
    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.set_xlabel('Pixel Number')
    ax1.set_ylabel('Total Counts Summed in the Dispersion Direction')
    img = pyfits.getdata(filename, ext)
    cenwave = pyfits.getval(filename, 'cenwave', 0)
    print type(cenwave)
    collapsed_img = collapse_spectrum(img)
    pix_num = range(len(collapsed_img))
    ax1.plot(pix_num, collapsed_img)
    if find_background_flag == 'y':
        disp_2d = raw_input('Would you like to display the 2D image? y, (n)')
        if disp_2d == 'y':
            fig2 = pyplot.figure(2)
            ax2 = fig2.add_subplot(1,1,1)
            ax2.imshow(np.rot90(np.log10(img)), interpolation = 'nearest', origin = 'lower')
        raw_input('Zoom in on desired spectrum location. Press Enter when finished')
        fig1, ax1, background_pix, background_collapsed = get_background_locations(fig1, ax1, pix_num, collapsed_img, filename, cenwave)

    else:
        coords = np.genfromtxt(filename.replace('.fits','_c%i_background_regions.txt'%(cenwave)))
        background_collapsed = np.empty((0,))
        background_pix = np.empty((0,))
        for region in coords:
            background_pix = np.append(background_pix, pix_num[int(region[0]):int(region[1])])
            background_collapsed = np.append(background_collapsed, collapsed_img[int(region[0]):int(region[1])])
        sort_indx = np.argsort(background_pix)
        background_pix = background_pix[sort_indx]
        background_collapsed = background_collapsed[sort_indx]

    xlims = ax1.get_xlim()
    ylims = ax1.get_ylim()
    ax1, fit = fit_background(ax1, background_pix, background_collapsed, pix_num)
    subtracted_img = subtract_2D_cross_disp_background(img, fit)
    write_temp_subtracted_file(subtracted_img, filename)
    sub_collapsed_img = collapse_spectrum(subtracted_img)
    ax1.cla()
    ax1.plot(pix_num, sub_collapsed_img)
    another_spectrum_flag = 'y'
    i = 0
    while another_spectrum_flag != 'n':
        c = colors[i%len(colors)]
        ax1.set_xlim(xlims[0], xlims[1])
        ax1.set_ylim(ylims[0], ylims[1])
        pyplot.draw()
        print 'Identify spectrum location'
        extrlocy = select_extraction_location(fig1, c)
        print 'Identify extraction box'
        extract_box_size = select_extraction_box_size(fig1, extrlocy, c)
        print 'Identify left background center'
        background_loc1 = select_extraction_location(fig1, c)
        print 'Idnetify left background box'
        background_size1 = select_extraction_box_size(fig1, background_loc1, c)
        print 'Identify right background center'
        background_loc2 = select_extraction_location(fig1, c)
        print 'Identify right background box'
        background_size2 = select_extraction_box_size(fig1, background_loc2, c)
        extract_spectrum(extrlocy, extract_box_size, background_loc1, background_size1, background_loc2, background_size2, filename, c)
        another_spectrum_flag = raw_input('Extract another spectrum?')
        i = i + 1
    os.remove(filename.replace('.fits', 'sub.fits'))




    raw_input('Enter')

