import os
import shutil
import glob
import sys
import pyfits
import numpy as np
import pdb
import time





def apply_herringbone_correction(input_dir, input_flist, output_dir = None):
    '''
    This function will copy files from an input folder to the ORIG folder in the herringbone
    correction folder. It will then run the herringbone correction script in batch mode on those files
    and copy the final result from the CLEAN folder back to the original folder
    '''
    if not output_dir:
        output_dir = input_dir
    #Clean out herringbone directory
    cur_dir = os.getcwd()
    os.chdir(herringbone_dir)
    flist = glob.glob('PSUB/*')+glob.glob('CLEAN/*') + glob.glob('TABLES/*') + glob.glob('ORIG/*') + glob.glob('FAILED/*')
    for ifile in flist:
        os.remove(ifile)    
    #Copy files to the ORIG directory in the herringbone_correction folder
    for ifile in input_flist:
        shutil.copy(os.path.join(input_dir, ifile), 'herringbone_correction/ORIG/%s' %(ifile))
    #Run autofilet
    os.chdir(herringbone_dir)
    batch_run_autofilet.autofilet()
    batch_run_autofilet.rebuild_files(interactive = True)
    #Move original raw files to the orig folder and corrected files to input_dir
    os.chdir(os.path.join(cur_dir, input_dir))
    os.mkdir('orig')
    for ifile in input_flist:
        shutil.copy(ifile, os.path.join(cur_dir, input_dir, 'orig/%s' %(ifile))) #copy uncorrected to orig
        shutil.copy(os.path.join(herringbone_dir, 'CLEAN/%s' %(ifile)), os.path.join(cur_dir, output_dir, ifile)) #copy corrected to input_dir
    os.chdir(cur_dir)
    print 'HERRINGBONE CORRECTION APPLIED '
    print '\tinput directory = %s' %(input_dir)
    print '\toutput directory = %s' %(output_dir)

def create_reference_file(input_dir, input_flist):
    '''
    This function will point Justin's refstis scripts to a folder containing darks or biases
    and have it create a superdark or super bias. The script will automatically try to create darks
    and biases for with all files in a directory. Biases and darks should be kept in separate 
    folders. This script automatically divided observations into individual weeks
    '''
    date = []
    for ifile in input_flist:
        date.append(pyfits.getval(os.path.join(input_dir, ifile), 'texpstrt', 0))
        date.append(pyfits.getval(os.path.join(input_dir, ifile),'texpend', 0))
    date = np.array(date)
    month_begin = np.min(date)
    month_end = np.max(date)
    REFSTIS_wrapper.separate_obs(input_dir, month_begin, month_end)
    start = time.time()
    REFSTIS_wrapper.make_ref_files(input_dir)
    end = time.time()
    print 'RUNTIME = %f minutes' %((end - start)/60.0) 
    #Copy reference file from each week folder to the reffiles folder
    os.path.walk(input_dir, copy_reference_files_to_reffiles, '')   

def copy_reference_files_to_reffiles(junk, dirname, junk2):
    '''
    This function copies super bias and super dark reference files to the reffiles folder
    '''
    flist = glob.glob(os.path.join(dirname, 'ref*.fits')) + glob.glob(os.path.join(dirname, 'week*.fits'))
    if len(flist) > 0:
        for ifile in flist:
            shutil.copy(ifile, os.path.join(cur_dir, 'reffiles', ifile.split('/')[-1]))
            print 'Copying %s to reffiles folder' %(ifile)

def calibrate_data_with_cte_correction(input_flist, header_tup = None):
    '''
    This function will:
    1. copy relevant reference files
    2. update relevant header keywords
    3. Start running calstis with BIASCORR, BLEVCORR, and DQICORR
    4. Run the CTE correction code
    5. Finish running calstis on the CTE corrected files with CTECORR turned off
    '''
    copy_relevant_reference_files(reffile_list, output_dir)
    update_header_keywords(input_flist, header_tup)
    run_calstis_part1(input_flist)
    run_cte_correction_code(input_flist)
    run_calstis_part2(input_flist)

def determine_correct_reference_files(input_dir, input_list):
    for ifile in flist:
        filename = os.path.join(input_dir, ifile)
        date = pyfits.getval(ifile, 'expstart', 0)
        

def copy_relevant_reference_files(reffile_list, output_dir):
    '''
    This file copies all reference files in reffile_list from reffiles to output_dir
    '''
    

def update_header_keywords(input_flist, header_tup):
    '''
    This function takes every header (keyword, value) tuple and updates every file in 
    input_flist
    '''
    pass

def run_calstis_part1(input_flist):
    '''
    This function runs BIASCORR, BLEVCORR and DQICORR on the files in input_flist
    '''
    pass

def run_cte_correction_code(input_flist):
    '''
    This function will run the CTE correction code StisPixCteCorr.py on each file in 
    input_flist
    '''
    pass

def run_calstis_part2(input_flist):
    '''
    This function finishes running calstis on the CTE corrected data. Note: CTECORR should be
    set to OMIT
    '''
    pass

if __name__ == "__main__":

    #Set paths
    #Eventually this can be done with argparse
    #set cur_dir, refstis_dir, 
    #What do I use programs for?

    cur_dir = '/user/bostroem/science/cte/2012_04'
    refstis_dir = '/user/bostroem/science/programs/refstis/'
    herringbone_dir = os.path.join(cur_dir, 'herringbone_correction')
    sys.path.append('/user/bostroem/science/programs')
    sys.path.append(herringbone_dir)
    sys.path.append(refstis_dir)
    import REFSTIS_wrapper
    import batch_run_autofilet
    os.chdir(cur_dir)

    #Herringbone correct Bias files
    os.chdir('bias')
    flist_bias = glob.glob('*raw.fits')
    os.chdir(cur_dir)
    ##apply_herringbone_correction('bias', flist_bias)

    
    #Create Superbias
    #create_reference_file('bias', flist_bias)
    
    #Herringbone correct Dark files
    os.chdir('dark')
    flist_dark = glob.glob('*raw.fits')
    os.chdir(cur_dir)
    apply_herringbone_correction('dark', flist_dark)
    '''pwd
    #CTE correct and calibrate Dark files with new superbias
    calibrate_data_with_cte_correction(flist_dark, header_tup = [('pcttab', 'stis_apr_2012_pcte.fits'), ('biasfile', 'superbias.fits')])
    #Create Superdark
    create_reference_file('dark')
    
    #Herringbone correct science data
    os.chdir('sci_data')
    flist_sci = glob.glob('*raw.fits')
    os.chdir(cur_dir)
    apply_herringbone_correction('sci_data', flist_sci)
    #CTE correcte and calibrate science data with new superbias and superdark
    calibrate_data_with_cte_correction(flist_sci, header_tup = [('pcttab', 'stis_apr_2012_pcte.fits'), ('biasfile', 'superbias.fits'), ('darkfile', 'superdark.fits')])
    '''