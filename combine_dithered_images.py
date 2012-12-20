import pyfits
import numpy as np
import os
import glob
import c_correlate 
from optparse import OptionParser

#!!!!!!!!!!!!!!!!!!!!!!
#CORRECTNESS OF DQ AND ERR ARRAYS HAVE NOT BEEN CONFIRMED ALTHOUGH THEY SHOULD FOLLOW THE SAME
#RULES AS THE IMAGE EXTENSION
#!!!!!!!!!!!!!!!!!!!!!!


def make_declination_dict(flist):
    '''
    Create a dictionary with declination as keys and filenames as lists to each key
    This is to match 2 dithered images together

    Create a dictionary with declination as keys and filenames as lists to each key
    This is to match 2 dithered images together

    Called from: Main
    Calls to: Nothing

    Input:
        flist: list of files to make dictionary out of
    Outputs:
        Nothing
    Returns:
        dec_dict: a dictionary with declination as keys and filenames as lists to each key
    '''
    dec_dict = {} 
    for ifile in flist:
        targ_dec = pyfits.getval(ifile, 'dec_targ', 0)
        if dec_dict.has_key(targ_dec):
            dec_dict[targ_dec].append(ifile)
        else:
            dec_dict[targ_dec] = [ifile]
    return dec_dict

def combine_dithered_images(dec_dict, targ_dec, offset):
    '''
    The 12465 data is dithered along the slit direction. This program finds the shifts
    of the dithered images by cross-correlating a middle column and combines them.

    Input:
        dec_dict: dictionary with declination as keys and a list of filenames of data taken at 
            the that declination as values
        targ_dec: decclination of target
        offset: dither offset in pixels
    Outputs:
        None
    Returns:
        new_img: combined image array
    '''

    #Read in 2 images
    ofile1 = pyfits.open(dec_dict[targ_dec][0])
    ofile2 = pyfits.open(dec_dict[targ_dec][1])
    hdr0 = ofile1[0].header 
    hdr1 = ofile1[1].header
    hdr2 = ofile1[2].header
    hdr3 = ofile1[3].header
    img1 = ofile1[1].data
    err1 = ofile1[2].data
    dq1 = ofile1[3].data
    img2 = ofile2[1].data
    err2 = ofile2[2].data
    dq2 = ofile2[3].data

    cenwave = hdr0['cenwave']
    exptime2 = pyfits.getval(dec_dict[targ_dec][1], 'texptime', 0)
    hdr0['texptime'] = hdr0['texptime'] + exptime2  #exptime is only in the Primary HDU
    if offset is None:
        #Get middle columns of each image
        col1 = np.sum(img1[:, 350:650], axis = 1)
        col2 = np.sum(img2[:, 350:650], axis = 1)
        #Calculate offset
        lag, corr_mat = c_correlate.c_corr(col1, col2)
        max_indx = np.argmax(corr_mat)
        print 'Found offset of %i' %(lag[max_indx])
        max_lag = lag[max_indx]
        if max_lag == 0:
            lag = np.append(lag[:max_indx], lag[max_indx+1:])
            corr_mat = np.append(corr_mat[:max_indx], corr_mat[max_indx+1:])
            max_indx = np.argmax(corr_mat)
            print 'Found new offset of %i' %(lag[max_indx])
    else:
        max_lag = offset
        print 'Using user input offset of %f' %(offset)
    if max_lag > 0:
        start1 = 0
        end1 = np.shape(img1)[0]
        start2 = max_lag
        end2 = np.shape(img1)[0] + max_lag
        new_img, new_err, new_dq = create_new_combined_arrays(img1, img2, err1, err2, dq1, dq2, start1, end1, start2, end2)

    else:
        start1 = 0
        end1 = np.shape(img1)[0]
        start2 = np.abs(max_lag)
        end2 = np.shape(img1)[0]+np.abs(max_lag)
        new_img, new_err, new_dq = create_new_combined_arrays(img2, img1, err2, err1, dq2, dq1, start1, end1, start2, end2)

    #Write fits file
    make_combined_fits(hdr0, hdr1, hdr2, hdr3, new_img, new_err, new_dq, cenwave)
 
    return new_img

def create_new_combined_arrays(top_img, bottom_img, top_err, bottom_err, top_dq, bottom_dq, start1, end1, start2, end2):
        new_img = np.zeros((end2, end1))
        new_img[start1:end1, :] = new_img[start1:end1, :] + bottom_img
        new_img[start2:end2, :] = new_img[start2:end2, :] + top_img

        new_err1 = np.zeros((end2, end1))
        new_err2 = np.zeros((end2, end1))
        new_err1[start1:end1, :] = new_err1[start1:end1, :] + bottom_err
        new_err2[start2:end2, :] = new_err2[start2:end2, :] + top_err
        new_err = np.sqrt(new_err1**2 + new_err2**2)

        new_dq1 = np.zeros((end2, end1))
        new_dq2 = np.zeros((end2, end1))
        new_dq1[start1:end1, :] = new_dq1[start1:end1, :] + bottom_dq
        new_dq2[start2:end2, :] = new_dq2[start2:end2, :] + top_dq
        new_dq1 = np.int_(new_dq1)
        new_dq2 = np.int_(new_dq2)
        new_dq = np.bitwise_or(new_dq1, new_dq2)
        return new_img, new_err, new_dq


def make_combined_fits(hdr0, hdr1, hdr2, hdr3, new_img, new_err, new_dq, cenwave):
    '''
    write combined image to file
    Called from:
        combine_dithered_images
    Calls to:
        Nothing
    Inputs:
        hdr: fits header object
        img: image array
        cenwave: central wavelength for filename
    Outputs:
        TARGNAME_(CENWAVE)_combined_img.fits files of the combined dithered images
    Returns:
        Nothing

    '''
    hdu0 = pyfits.PrimaryHDU()
    hdu0.header = hdr0
    print '#############################%s######################' %(hdr0['targname'])
    hdu1 = pyfits.ImageHDU(data=new_img, name='SCI', header = hdr1)
    hdu2 = pyfits.ImageHDU(data=new_err, name='ERR', header = hdr2)
    hdu3 = pyfits.ImageHDU(data = np.int16(new_dq), name = 'DQ', header = hdr3)
    hdulist = pyfits.HDUList([hdu0, hdu1, hdu2, hdu3])
    if cenwave:
        hdulist.writeto('%s_%i_combined_img.fits' %(hdr0['targname'], hdr0['cenwave']), clobber = True)
    else:
        hdulist.writeto('%s_combined_img.fits' %(hdr0['targname']), clobber = True)



if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--offset', dest = 'offset', type = 'int', help = 'Set the dithered offset in pixels')
    (options, args) = parser.parse_args()
    
    #for the FUV data
    idir = '/user/bostroem/science/12465_otfr20120425/mama/'
    os.chdir(idir)
    flist = glob.glob('obrc04???_flt.fits')+glob.glob('obrc05???_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.offset)

