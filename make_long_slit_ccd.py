import pyfits
from matplotlib import pyplot
import glob
import numpy as np
import pdb
import math


def plot_ax(img, orientat,  ax,targname, vmin, vmax, ymin, ymax, offset):
    ax.set_title(targname, fontsize = 9)
    l = ax.imshow(img, cmap = 'bone', vmin = vmin, vmax = vmax, aspect = 'auto', interpolation = 'nearest')
    ax.set_ylim(ymin, ymax, offset)
    if orientat < 0:
        ax.cla()
        ax.set_title(targname, fontsize = 9)
        ax.imshow(np.rot90(img, k = 2), cmap = 'bone', vmin = vmin, vmax = vmax, aspect = 'auto', interpolation = 'nearest')
        ax.set_ylim(ymin - offset, ymax - offset)
    ax.set_xticks([25])
    ax.set_xticklabels(['%4.1f' %(orientat)], fontsize = 7)
    return ax, l

def make_long_slit_img(flist, title, ymin, ymax, vmin, vmax, offset):
    fig = pyplot.figure(figsize = [8, 17])
    ax1 = fig.add_subplot(1, 17, 1)
    ax2 = fig.add_subplot(1, 17, 2)
    ax3 = fig.add_subplot(1, 17, 3)
    ax4 = fig.add_subplot(1, 17, 4)
    ax5 = fig.add_subplot(1, 17, 5)
    ax6 = fig.add_subplot(1, 17, 6)
    ax7 = fig.add_subplot(1, 17, 7)
    ax8 = fig.add_subplot(1, 17, 8)
    ax9 = fig.add_subplot(1, 17, 9)
    ax10 = fig.add_subplot(1, 17, 10)
    ax11 = fig.add_subplot(1, 17, 11)
    ax12 = fig.add_subplot(1, 17, 12)
    ax13 = fig.add_subplot(1, 17, 13)
    ax14 = fig.add_subplot(1, 17, 14)
    ax15 = fig.add_subplot(1, 17, 15)
    ax16 = fig.add_subplot(1, 17, 16)
    ax17 = fig.add_subplot(1, 17, 17)
    
    
    for ifile in flist:
        print ifile
        #Get the middle 50 columns of the image
        img = pyfits.getdata(ifile, 1)[:, 487:537]
        orientat = pyfits.getval(ifile, 'orientat', 1)
        if 'SE9' in ifile:
            ax1, l1 = plot_ax(img, orientat,  ax1, 'SE9', vmin, vmax, ymin, ymax, offset)
        if  'SE8' in ifile:
            ax2, l2 = plot_ax(img, orientat, ax2, 'SE8', vmin, vmax, ymin, ymax, offset)
        if  'SE7' in ifile:
            ax3, l3 = plot_ax(img, orientat,  ax3, 'SE7', vmin, vmax, ymin, ymax, offset)
        if  'SE6' in ifile:
            ax4, l4 = plot_ax(img, orientat,  ax4, 'SE6', vmin, vmax, ymin, ymax, offset)
        if  'SE5' in ifile:
            ax5, l5 = plot_ax(img, orientat,  ax5, 'SE5', vmin, vmax, ymin, ymax, offset)
        if  'SE4' in ifile:
            ax6, l6 = plot_ax(img, orientat,  ax6, 'SE4', vmin, vmax, ymin, ymax, offset)
        if  'SE3' in ifile:
            ax7, l7 = plot_ax(img, orientat,  ax7, 'SE3', vmin, vmax, ymin, ymax, offset)
        if  'SE2' in ifile:
            ax8, l8 = plot_ax(img, orientat,  ax8, 'SE2', vmin, vmax, ymin, ymax, offset)
        if  'SE1' in ifile:
             ax9, l9 = plot_ax(img, orientat,  ax9, 'SE1', vmin, vmax, ymin, ymax, offset)
        if  'NW1' in ifile:
            ax10, l10 = plot_ax(img, orientat,  ax10, 'NW1', vmin, vmax, ymin, ymax, offset)
        if  'NW2' in ifile:
            ax11, l11 = plot_ax(img, orientat,  ax11, 'NW2', vmin, vmax, ymin, ymax, offset)
        if  'NW3' in ifile:
            ax12, l12 = plot_ax(img, orientat,  ax12, 'NW3', vmin, vmax, ymin, ymax, offset)
        if  'NW4' in ifile:
            ax13, l13 = plot_ax(img, orientat,  ax13, 'NW4', vmin, vmax, ymin, ymax, offset)
        if  'NW5' in ifile:
            ax14, l14 = plot_ax(img, orientat,  ax14, 'NW5', vmin, vmax, ymin, ymax, offset)
        if  'NW6' in ifile:
            ax15, l15 = plot_ax(img, orientat,  ax15, 'NW6', vmin, vmax, ymin, ymax, offset)
        if  'NW7' in ifile:
            ax16, l16 = plot_ax(img, orientat,  ax16, 'NW7', vmin, vmax, ymin, ymax, offset)
        if  'NW8' in ifile:
            ax17, l17 = plot_ax(img, orientat,  ax17, 'NW8', vmin, vmax, ymin, ymax, offset)

    ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17]
    for ax in ax_list[1:]:
        ax.set_yticks([])
    #ax1.set_xticks([])
    
    fig.suptitle(title)
    return ax_list, fig

#View center of cluster
ymin = 450
ymax = 580
#set vmax and vmin
vmax = 1000
vmin = 0

#orientat1 = ['SE9', 'SE8', 'SE7', 'SE6', 'SE5', 'SE4', 'SE2', 'SE1', 'NW1']
#orientat2 = ['NW2', 'NW3', 'NW4', 'NW5', 'NW6', 'NW7', 'NW8']

#cenwaves = ['3963', '4194', '4451', '4706']

#3963: SE9:NW1 @64, NW2:NW8 @-115
#4194: SE9:NW1 @ -115, NW2:NW8 @64
#4451: SE9:NW1 @64, NW2:NW8 @-115
#4706: SE9:NW1 @64, NW2:NW8 @-115
#6581: all @-115

#derived this for #3963, will have to see if it works for the other gratings: 1031 - yloc_img + 13

flist = glob.glob('/Users/bostroem/science/12465_otfr20121109/ccd/???_3936_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_3936_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_3936_combined_img.fits.gz')
title = 'Long Slit Image G430M 3963'
ax_list3963, fig3963 = make_long_slit_img(flist, title, ymin, ymax, vmin, vmax, 20)
fig3963.savefig('/Users/bostroem/science/'+'_'.join(title.split())+'.pdf')


flist = glob.glob('/Users/bostroem/science/12465_otfr20121109/ccd/???_4194_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_4194_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_4194_combined_img.fits.gz')
title = 'Long Slit Image G430M 4194'
ax_list4194, fig4194 = make_long_slit_img(flist, title, ymin, ymax, vmin, vmax, 7)
fig4194.savefig('/Users/bostroem/science/'+'_'.join(title.split())+'.pdf')

flist = glob.glob('/Users/bostroem/science/12465_otfr20121109/ccd/???_4451_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_4451_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_4451_combined_img.fits.gz')
title = 'Long Slit Image G430M 4451'
ax_list4451, fig4451 = make_long_slit_img(flist, title, ymin, ymax, vmin, vmax, 20)
fig4451.savefig('/Users/bostroem/science/'+'_'.join(title.split())+'.pdf')

flist = glob.glob('/Users/bostroem/science/12465_otfr20121109/ccd/???_4706_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_4706_combined_img.fits') + \
    glob.glob('/Users/bostroem/science/12465_otfr20130503/ccd/???_4706_combined_img.fits.gz')
title = 'Long Slit Image G430M 4706'
ax_list4706, fig4706 = make_long_slit_img(flist, title, ymin, ymax, vmin, vmax, 15)
fig4706.savefig('/Users/bostroem/science/'+'_'.join(title.split())+'.pdf')





