import pyfits
import numpy as np
import sys
import math
from matplotlib import pyplot
import pdb
import os
import shutil
from optparse import OptionParser


import pyraf
from pyraf import iraf
from iraf import stsdas,hst_calib,stis,x1d


def collapse_spectrum(img, num_cols):
    '''Collapse 2D spectrum in the dispersion direction to create a cross dispersion profile'''
    collapsed_img = np.sum(img[:, int(512 - num_cols/2.0):int(512 + num_cols/2.0)], axis = 1)/num_cols
    return collapsed_img

def get_background_locations(fig1, ax1, pix_num, collapsed_img, filename, cenwave, first_time = True):
    '''Interactively select the regions to be fit as background '''
    finished_flag = 'n'
    background_collapsed = np.empty((0,))
    background_pix = np.empty((0,))
    #no background regions defined
    if first_time is True:
        if os.path.exists(filename.replace('.fits', '_c%i_background_regions.txt' %(cenwave))):
            keep_file_flag = raw_input('Background file already exists: Overwrite file (w) or Append to file (a)')
        else:
            keep_file_flag = 'w'
        ofile = open(filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), keep_file_flag)
    #Background regions have been defined - adding or removing some of them
    if first_time is False:
        remove_flag = 'y'
        #Get initial indx values from background file
        start_indx, end_indx  = np.genfromtxt(filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), unpack = True)
        #REMOVE BACKGROUND REGIONS
        while remove_flag != 'n':
            l2 = [] #list for line objects (this is so we can make them invisible later)
            l3 = [] #list of text objects (this is so we can make them invisible later)
            #Draw background regions as currently defined
            for i, sindx, eindx in zip(range(len(start_indx)), start_indx, end_indx):
                temp_indx = np.where((pix_num <= eindx) & (pix_num >= sindx))[0]
                temp_pix = pix_num[temp_indx]
                temp_bkg = collapsed_img[temp_indx]
                l2.append(pyplot.plot([pix_num[temp_indx[0]], pix_num[temp_indx[-1] + 1]], [collapsed_img[temp_indx[0]], collapsed_img[temp_indx[-1]+1]], 'mo--', markersize = 5))
                l3.append(pyplot.text(np.mean(temp_pix), np.max(temp_bkg), str(i)))
            remove_flag = raw_input('Would you like to remove any background regions? (y), n ')
            if remove_flag == 'n':
                break
            try:
                remove_indx = int(raw_input('Enter the number of the region you wish to remove: '))
            except ValueError:
                print 'You must enter something that can be converted into an integer'
                remove_indx = int(raw_input('Enter the number of the region you wish to remove: '))
            #remove indx from start_indx and end_indx
            start_indx = np.append(start_indx[0:remove_indx], start_indx[remove_indx+1:])
            end_indx = np.append(end_indx[0:remove_indx], end_indx[remove_indx+1:])
            #delete lines from plot
            for iline, itxt in zip(l2, l3):
                iline[0].set_visible(False)
                itxt.set_visible(False)
        ofile = open(filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), 'w')
        for sindx, eindx in zip(start_indx, end_indx):
            temp_indx = np.where((pix_num <= eindx) & (pix_num >= sindx))[0]
            background_collapsed = np.append(background_collapsed, collapsed_img[temp_indx[0]:temp_indx[-1]+1])
            background_pix = np.append(background_pix, pix_num[temp_indx[0]:temp_indx[-1]+1])
            ofile.write('%i\t %i \n' %(int(sindx),int(eindx)))
        finished_flag = raw_input('Finished entering points? (n), y ')
    #ADD BACKGROUND REGIONS
    while finished_flag != 'y':
        print 'Select background points for cross-dispersion background subtraction'
        temp_region = fig1.ginput(n = 2, timeout = -1)
        back_reg = [temp_region[0][0], temp_region[1][0]]
        l1 = pyplot.plot(back_reg, np.interp(np.array(back_reg), pix_num, collapsed_img), 'r|-', markersize = 5)
        pyplot.draw()
        keep_flag = raw_input('Keep points? (y), n ')
        if keep_flag != 'n':
            ofile.write('%i\t %i \n' %(int(math.floor(temp_region[0][0])),int(math.ceil(temp_region[1][0]))))
            background_collapsed = np.append(background_collapsed, collapsed_img[int(math.floor(temp_region[0][0])):int(math.ceil(temp_region[1][0]))+1])
            background_pix = np.append(background_pix, pix_num[int(math.floor(temp_region[0][0])):int(math.ceil(temp_region[1][0]))+1])
        else:
            l1[0].set_visible(False)
            pyplot.draw()
        finished_flag = raw_input('Finished entering points? (n), y ')
    ofile.close()
    sort_indx = np.argsort(background_pix)
    background_pix = background_pix[sort_indx]
    background_collapsed = background_collapsed[sort_indx]

    return fig1, ax1, background_pix, background_collapsed

def fit_background(fig1, ax1, background_pix, background_collapsed, pix_num):
    '''Interactively fit background regions with a polynomial '''
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
        change_deg_flag = raw_input('Replot with a different degree polynomial? (y), n -or- redefine background, r ')
        if (change_deg_flag != 'n') and (change_deg_flag != 'r'):
            l1[0].set_linestyle('--')
            pyplot.draw()
            try:
                deg = int(raw_input('Enter degree of polynomial you wish to fit: '))
            except ValueError:
                print 'The value you entered cannot be converted to an integer, try again'
                deg = int(raw_input('Enter degree of polynomial you wish to fit: '))

            leg_text.append('%i' %(deg))
        elif change_deg_flag == 'r':
            for l in leg_lines:
                l.set_visible(False)
            return None, None
    return ax1, fit

def subtract_2D_cross_disp_background(img, fit):
    '''Create a 2D image from the 1D fit to the background in the cross dispersion direction and subtract it from the image '''
    two_D_background = np.rot90(np.rot90(np.rot90(np.tile(fit, (1024, 1)))))
    return img - two_D_background

def write_temp_subtracted_file(subtracted_img, filename):
    '''Write the subtracted image to a temporary file for calstis to extract from. This file is erased at the end of execution '''
    shutil.copyfile(filename, filename.replace('.fits', 'sub.fits'))
    ofile = pyfits.open(filename.replace('.fits', 'sub.fits'), mode = 'update')
    ofile[1].data = subtracted_img
    ofile.flush()
    ofile.close()

def select_extraction_location(fig1, c):
    '''Interactvely select the center of an extraction or background region '''
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
    '''Interactively select the box size of an extraction or background region '''
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

def extract_spectrum(extrlocy, extract_box_size, background_loc1, background_size1, background_loc2, background_size2, filename, c, backcorr_option, bksmode_option):
    '''Extract the spectrum using the calstis task x1d '''
    if os.path.exists(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy))))):
        os.remove(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))))
        
    iraf.stsdas.hst_calib.stis.x1d(filename.replace('.fits', 'sub.fits'), output = filename.replace('.fits', '_loc%i.fits' %(int(extrlocy))), \
                                       a2center = extrlocy + 1, extrsize = extract_box_size, maxsrch = 0, bk1offst = background_loc1 - extrlocy, \
                                       bk2offst = background_loc2 - extrlocy, bk1size = background_size1, bk2size = background_size2, backcorr = backcorr_option, bksmode = bksmode_option)  #a2center is 1 indexed
    add_background_into_x1d(filename, fit, extrlocy)

def add_background_into_x1d(filename, fit, extrlocy):
    '''Add the fitted background back into the gross and background columns of the extracted spectrum'''
    ofile = pyfits.open(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))), mode = 'update')
    ofile[1].data['gross'][:] = ofile[1].data['gross'][:] + fit[int(round(ofile[1].data['a2center'])) - 1] #a2center is 1 indexed
    ofile[1].data['background'][:] = ofile[1].data['background'][:] + fit[int(round(ofile[1].data['a2center'])) - 1] #a2center is 1 indexed
    ofile.flush()
    ofile.close()

def confirm_extraction_location(filename, img, extrlocy):
    fig3 = pyplot.figure(3)
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.imshow(np.log10(img), interpolation = 'nearest', origin = 'lower', cmap = 'gray')
    tbdata =  pyfits.getdata(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))),  1)
    extrlocy = tbdata['extrlocy'].ravel() - 1.0
    extrsize = tbdata['extrsize'].ravel() 
    bk1size = tbdata['bk1size'].ravel()
    bk2size = tbdata['bk2size'].ravel()
    bk1offset = tbdata['bk1offst'].ravel()
    bk2offset = tbdata['bk2offst'].ravel()
    pix = np.arange(len(extrlocy))
    pyplot.plot(pix, extrlocy, 'r', lw = 3)
    pyplot.plot(pix, extrlocy - 0.5*extrsize, 'r--', lw = 2)
    pyplot.plot(pix, extrlocy + 0.5*extrsize, 'r--', lw = 2)
    pyplot.plot(pix, extrlocy + bk1offset, color = '#0022FF', lw = 2)
    pyplot.plot(pix, extrlocy + bk1offset + 0.5*bk1size, color = '#0022FF', ls = '--', lw = 2)
    pyplot.plot(pix, extrlocy + bk1offset - 0.5 * bk1size, color = '#0022FF', ls = '--', lw = 2)
    pyplot.plot(pix, extrlocy + bk2offset, color = '#0022FF', lw = 2)
    pyplot.plot(pix, extrlocy + bk2offset + 0.5*bk2size, color = '#0022FF', ls = '--', lw = 2)
    pyplot.plot(pix, extrlocy + bk2offset - 0.5*bk2size, color = '#0022FF', ls = '--', lw = 2)
    raw_input('Press enter to close and continue')
    pyplot.close(fig3)
    
    
if __name__ == "__main__":

    #Define colors for extracting more than one spectrum
    colors = ['r', 'g', 'c', 'k', 'm']

    pyplot.ion()  #turn plotting on 
    #os.environ['oref'] = '/grp/hst/cdbs/oref/' #set oref environment variable to point to reference file location
    os.environ['oref'] = '/Users/bostroem/science/oref/'
    parser = OptionParser()
    parser.add_option('--backcorr', dest = 'backcorr', help = 'Enter perform (default) or omit to perform or omit the background subtraction in CalSTIS x1d', default = 'perform')
    parser.add_option('--ncol', dest = 'num_cols', type = 'float', help = 'Number of columns summed when examining the cross-dispersion profile', default = 50)
    parser.add_option('--backsmooth', dest = 'bksmode', help = 'Background smoothing mode: off, median, average', default = 'off')
    
    (options, args) = parser.parse_args()
    
    #read in command line arguments
    filename = sys.argv[1]
    try:
        ext = int(sys.argv[2])
    except:
        ext = 1
    cenwave = pyfits.getval(filename, 'cenwave', 0)
    if os.path.exists(os.path.join(os.getcwd(), filename.replace('.fits','_c%i_background_regions.txt'%(cenwave)))):
        find_background_flag = raw_input('Select the background regions? y, (n) ') #read background from file or interactively define?
    else:
        find_background_flag = 'y'

    #set up figure and plot cross-dispersion profile
    fig1 = pyplot.figure(1)
    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.set_xlabel('Pixel Number')
    ax1.set_ylabel('Total Counts Summed in the Dispersion Direction')
    #pdb.set_trace()
    img = pyfits.getdata(filename, ext)
    
    collapsed_img = collapse_spectrum(img, options.num_cols)
    pix_num = np.arange(len(collapsed_img))
    ax1.plot(pix_num, collapsed_img)
    
    if find_background_flag == 'y':  #interactvely define background
        #Display 2D image
        disp_2d = raw_input('Would you like to display the 2D image? y, (n) ')
        if disp_2d == 'y':
            fig2 = pyplot.figure(2)
            ax2 = fig2.add_subplot(1,1,1)
            ax2.imshow(np.rot90(np.log10(img)), interpolation = 'nearest', origin = 'lower')
        raw_input('Zoom in on desired spectrum location. Press Enter when finished')
        
        #Define background regions in cross dispersion profile
        fig1, ax1, background_pix, background_collapsed = get_background_locations(fig1, ax1, pix_num, collapsed_img, filename, cenwave)

    else:  #read background regions from file
        coords = np.genfromtxt(filename.replace('.fits','_c%i_background_regions.txt'%(cenwave)))
        background_collapsed = np.empty((0,))
        background_pix = np.empty((0,))
        for region in coords:
            background_pix = np.append(background_pix, pix_num[int(region[0]):int(region[1])])
            background_collapsed = np.append(background_collapsed, collapsed_img[int(region[0]):int(region[1])])
        sort_indx = np.argsort(background_pix)
        background_pix = background_pix[sort_indx]
        background_collapsed = background_collapsed[sort_indx]

    #Get limits now to scale plot when plotting the background subtracted data (things can go crazy at the edges)
    xlims = ax1.get_xlim()
    ylims = ax1.get_ylim()
    
    #Fit a polynomial to the background
    ax1_temp, fit = fit_background(fig1, ax1, background_pix, background_collapsed, pix_num)
    while ax1_temp is None:
        fig1, ax1, background_pix, background_collapsed = get_background_locations(fig1, ax1, pix_num, collapsed_img, filename, cenwave, first_time = False)
        ax1_temp, fit = fit_background(fig1, ax1, background_pix, background_collapsed, pix_num)
    ax1 = ax1_temp
    #Subtract the background (in the cross dispersion direction) from the image
    subtracted_img = subtract_2D_cross_disp_background(img, fit)
    #Write a temporary file of the subtracted image
    write_temp_subtracted_file(subtracted_img, filename)
    #collapse the subtracted image in the cross-dispersion direction
    sub_collapsed_img = collapse_spectrum(subtracted_img, options.num_cols)
    ax1.cla()
    ax1.plot(pix_num, sub_collapsed_img) #Plot the cross dispersion profile of the subtracted image
    ax1.axhline(0, ls = '--', color = 'k')
    another_spectrum_flag = 'y'  #this flag denotes wanting to extract more than one spectrum from a subtracted image
    i = 0
    while another_spectrum_flag != 'n':  #for each spectrum
        c = colors[i%len(colors)] #Set color for spectrum extraction regions
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

        extract_spectrum(extrlocy, extract_box_size, background_loc1, background_size1, background_loc2, background_size2, filename, c, options.backcorr, options.bksmode)
        confirm_flag = raw_input('Would you like to confirm the location of your extraction on a 2D image? (y), n ')
        if confirm_flag != 'n':
            confirm_extraction_location(filename, img, extrlocy)
        another_spectrum_flag = raw_input('Extract another spectrum? ')
        i = i + 1
    os.remove(filename.replace('.fits', 'sub.fits'))



#Add log file which records degree polynomial fit, extraction locations and box sizes

