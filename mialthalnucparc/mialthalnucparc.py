"""Main module."""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import bids
from bids import BIDSLayout
import nibabel as nib
from glob import glob
from operator import itemgetter
from datetime import datetime
import argparse
import csv
import json
import subprocess
import scipy.ndimage as sc
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from concurrent.futures import as_completed
import time
import copy
from threading import Lock
from importlib import reload

from templateflow import api
from templateflow import api as tflow
import clabtoolkit.misctools as cltmisc
import clabtoolkit.freesurfertools as cltfree
import clabtoolkit.parcellationtools as cltparc
from rich.progress import Progress

def _build_args_parser():

    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)

    from argparse import ArgumentParser

    p = argparse.ArgumentParser(formatter_class=cltmisc.SmartFormatter, description='\n Help \n')

    requiredNamed = p.add_argument_group('Required arguments')

    requiredNamed.add_argument('--bidsdir', '-b', action='store', required=False, metavar='BIDSDIR', type=str, nargs=1,
                                help="R| BIDs dataset folder. \n"
                                "\n")
    requiredNamed.add_argument('--derivdir', '-d', action='store', required=False, metavar='DERIVDIR', type=str, nargs=1,
                                help="R| BIDs derivative folder containing the derivatives folder. \n"
                                "\n",
                                default='None')
    
    requiredNamed.add_argument('--nthreads', '-n', action='store', required=False, metavar='NTHREADS', type=str, nargs=1,
                                help="R| Number of processes to run in parallel (default= Number of cores - 4). \n", default=['auto'])

    requiredNamed.add_argument('--t1s', '-t', action='store', required=False, metavar='T1FILE', type=str, nargs=1,
                                help="R| File containing the basename of the NIFTI images that will be ran. \n"
                                    "   This file is useful to tun Chimera, only, on certain T1s in case of multiple T1s \n"
                                    " for the same session.\n"
                                    " Example of this file: \n"
                                    "   sub-00001_ses-0001_run-2 \n"
                                    "   sub-00001_ses-0003_run-1\n"
                                    "   sub-00001_ses-post_acq-mprage\n"
                                    " \n", default=None)
    requiredNamed.add_argument('--force', '-f', action='store_true', required=False,
                                help="R| Overwrite the results. \n"
                                    "\n")
    p.add_argument('--verbose', '-v', action='store', required=False,
                    type=int, nargs=1,
                    help='verbosity level: 1=low; 2=debug')

    args = p.parse_args()

    if args.derivdir is None or args.bidsdir is None:
        print('--bidsdir and --derivdir are REQUIRED arguments')
        sys.exit()

    bids_dir = args.bidsdir[0]
    deriv_dir = args.derivdir[0]

    if not os.path.isdir(bids_dir):
        print("\n")
        print("Please, supply a valid BIDs directory.")
        p.print_help()
        sys.exit()

    return p

def _compute_abased_thal_parc(t1, deriv_dir, conf_file:str = None, thal_tsv:str = None, prob_thresh:float = 0.05):
    """
    Compute atlas-based thalamic parcellation. 
    
    Parameters:
    ----------
    t1 : str
        T1 image file
        
    vol_tparc : str
        Output thalamic parcellation file
        
    deriv_dir : str
        Derivative directory
        
    pathcad : str
        Path to the CAD directory
        
    fullid : str
        Full subject id
        
    aseg_nii : str
        ASEG file
        
    out_str : str
        Output string
        
    Returns:
    --------
    mial_thalparc : str
        MIAL thalamic parcellation file
        
    
    """
    
    global layout

    # Detecting the entities of the T1 image
    ent_dict = layout.parse_file_entities(t1)
    
    if 'session' in ent_dict.keys():
        pattern_fullid = "sub-{subject}_ses-{session}_run-{run}"
        path_cad       = "sub-" + ent_dict["subject"] + os.path.sep + "ses-" + ent_dict["session"]
    else:
        pattern_fullid = "sub-{subject}_run-{run}"
        path_cad       = "sub-" + ent_dict["subject"]
    
    
    anat_dir = os.path.dirname(t1)
    t1_name = os.path.basename(t1)
    t1name_ent_dict = t1_name.split('_')[:-1]
    full_id = "_".join(t1name_ent_dict)

    
    ######## -- Reading the configuration dictionary  ------------ #
    # Get the absolute of this file
    cwd = os.path.dirname(os.path.abspath(__file__))
    cwd = os.path.dirname(cwd)
    if conf_file is None:
        pipe_json = os.path.join(cwd, 'config', 'pipe_config.json')
    else:
        if os.path.isfile(conf_file):
            pipe_json = conf_file
        else:
            raise ValueError("Please, provide a valid JSON file containing the pipeline configuration dictionary.")
    
    with open(pipe_json) as f:
        config_dict = json.load(f)
        
    # Reading the tsv with the parcellation information
    
    if thal_tsv is None:
        thal_tsv = os.path.join(cwd, 'config', 'Thalamus.tsv')
    else:
        if not os.path.isfile(thal_tsv):
            raise ValueError("Please, provide a valid TSV file containing the thalamic nuclei information.")
        
    # Reading the tsv file containing the thalamic nuclei information (index, name, color)
    temp_df = pd.read_csv(thal_tsv, sep='\t')
    
    # Get the hemispheres
    st_hemi = temp_df["hemi"].tolist()
    
    # Get the unique hemispheres
    st_hemi = np.unique(st_hemi).tolist()
    
    info_dict = {} # Create a dictionary for each hemisphere
    for hemi in st_hemi:
        
        # Create a sub dataframe for each hemisphere
        sub_df3 = temp_df.loc[temp_df["hemi"] == hemi]
        
        # Get the indexes of the regions
        indexes = sub_df3["index"].tolist()
        
        # Get the names of the regions
        names = sub_df3["name"].tolist()
        
        # Get the colors of the regions
        colors = sub_df3["color"].tolist()
        
        # Create a dictionary for each hemisphere
        temp_dict = {"index": indexes, "name": names, "color": colors}
        
        # Add the dictionary to the info_dict
        info_dict[hemi] = temp_dict
    
    
    ######## ------------ Selecting the templates  ------------ #
    if config_dict["templates"]["reference"]["tool"] == "templateflow":
        ###########################################################
        ################ Downloading the templates ################
        ###########################################################
        
        # Setting templateflow home directory
        tflow_dir = config_dict["packages"]["templateflow"]["home_dir"]
        
        if tflow_dir == "local":
            
            tflow_dir = os.environ.get('TEMPLATEFLOW_HOME')
            
            if tflow_dir is None:
                # Setting the templateflow home directory in the same directory as the script
                temp_dir = os.path.dirname(os.path.realpath(__file__))
                
                # Select the directory before 
                temp_dir = os.path.dirname(temp_dir)
                tflow_dir = os.path.join(temp_dir, 'templateflow')
        
        # Create the directory if it does not exist using the library Path
        tflow_dir = Path(tflow_dir)
        
        # If the directory does not exist create the directory and if it fails because it does not have write access send an error
        try:
            tflow_dir.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            print("The TemplateFlow directory does not have write access.")
            sys.exit()
            
        if os.path.isdir(tflow_dir):        
            # Setting the templateflow home directory
            os.environ["TEMPLATEFLOW_HOME"] = tflow_dir.as_posix()
            reload(api)
            
        else:
            print("The TemplateFlow directory does not exist.")
            sys.exit()
            
        # Getting the templates
        # Reference space
        temp_cad = config_dict["templates"]["reference"]["space"]
        t1_temp = tflow.get(temp_cad, desc=None, resolution=1, suffix='T1w', extension='nii.gz')
        
        # Getting the thalamic nuclei spams 
        atlas_cad = config_dict["templates"]["spams"]["atlas"]
        thal_spam = tflow.get(temp_cad, desc=None, resolution=1,atlas=atlas_cad, suffix='probseg', extension='nii.gz')
        
    else:
        t1_temp = config_dict["templates"]["reference"]["space"]
        if not os.path.isfile(t1_temp):
            print("The template file does not exist.")
            sys.exit()
        
        thal_spam = config_dict["templates"]["spams"]["atlas"]
        if not os.path.isfile(thal_spam):
            print("The thalamic atlas file does not exist.")
            sys.exit()
        
        temp_cad = "CustomSpace"
        atlas_cad = "CustomParc"
    
    # thal_spam = Path('/media/COSAS/Yasser/Useful_templates/thalamic_nuclei_MIALatlas/Thalamus_Nuclei-HCP-4DSPAMs.nii.gz')
    
    ######## -- Registration to the template space  ------------ #
    # Creating spatial transformation folder
    stransf_dir = Path(os.path.join(deriv_dir, config_dict["outputs"]["transforms"], path_cad, 'anat'))
    
    # If the directory does not exist create the directory and if it fails because it does not have write access send an error
    try:
        stransf_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        print("The directory to store the spatial transformations does not have write access.")
        sys.exit()
            
    # If 'space' is a substring in any of the entities, remove it
    temp_ent = [s for s in t1name_ent_dict if 'space' not in s]
    
    # Base ID for the transformation
    xfm_id = '_'.join(temp_ent)
    
    # Spatial transformation files
    defFile = os.path.join(str(stransf_dir), xfm_id + '_space-' + temp_cad + '_desc-t12mni_')
    xfm_affine = os.path.join(str(stransf_dir), xfm_id + '_space-' + temp_cad + '_desc-t12mni_0GenericAffine.mat')
    xfm_nl = os.path.join(str(stransf_dir), xfm_id + '_space-' + temp_cad + '_desc-t12mni_Warp.nii.gz')
    xfm_invnl = os.path.join(str(stransf_dir), xfm_id + '_space-' + temp_cad + '_desc-t12mni_1InverseWarp.nii.gz')
    
    
    if not os.path.isfile(xfm_invnl):
        # Registration to MNI template
        
        cmd_bashargs = ['antsRegistrationSyNQuick.sh', '-d', '3', '-f', t1_temp, '-m', t1, '-t', 's',
                        '-o', defFile]
        
        cont_tech = config_dict["packages"]["ants"]["con_tech"]
        cont_image = config_dict["packages"]["ants"]["container"]
        cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
    
    ######## -- Applying the transform ------------ #
    thalparc_dir = Path(os.path.join(deriv_dir, config_dict["outputs"]["thalparc"], path_cad, 'anat'))
    
    temp_ent = [s for s in t1name_ent_dict if 'atlas' not in s]
    altas_id = '_'.join(temp_ent)
    
    thalparc_maxprob = os.path.join(thalparc_dir, altas_id + '_atlas-' + atlas_cad + '_dseg.nii.gz')
    thalparc_lutfile = os.path.join(thalparc_dir, altas_id + '_atlas-' + atlas_cad + '_dseg.lut')
    thalparc_tsvfile = os.path.join(thalparc_dir, altas_id + '_atlas-' + atlas_cad + '_dseg.tsv')
    thalparc_spam    = os.path.join(thalparc_dir, altas_id + '_atlas-' + atlas_cad + '_probseg.nii.gz')
    
    # If the directory does not exist create the directory and if it fails because it does not have write access send an error
    try:
        thalparc_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        print("The directory to store the spatial transformations does not have write access.")
        sys.exit()
    
    if not os.path.isfile(thalparc_spam):
        # Applying spatial transform
        cmd_bashargs = ['antsApplyTransforms', '-d', '3', '-e', '3', '-i', thal_spam,
                        '-o', thalparc_spam, '-r', t1, '-t', xfm_invnl,
                        '-t','[' + xfm_affine + ',1]', '-n', 'Linear']
    
        cont_tech = config_dict["packages"]["ants"]["con_tech"]
        cont_image = config_dict["packages"]["ants"]["container"]
        cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
        
        # Removing the Warped images
        iwarped = os.path.join(str(stransf_dir), xfm_id + '_space-' + temp_cad + '_desc-t12mni_1InverseWarped.nii.gz')
        if os.path.isfile(iwarped):
            os.remove(iwarped)
            
        warped = os.path.join(str(stransf_dir), xfm_id + '_space-' + temp_cad + '_desc-t12mni_Warped.nii.gz')
        if os.path.isfile(warped):
            os.remove(warped)
        
    ######## -- Running FreeSurfer to obtain the thalamic masks ------------ #
    cont_tech = config_dict["packages"]["freesurfer"]["con_tech"]
    cont_image = config_dict["packages"]["freesurfer"]["container"]
    
    
    fssubj_dir = os.path.join(deriv_dir, config_dict["outputs"]["freesurfer"])
    aseg_mgz = os.path.join(fssubj_dir, full_id, 'mri', 'aseg.auto.mgz')
    aseg_nii = os.path.join(fssubj_dir, full_id, 'tmp', 'aseg.nii')
    
    if not os.path.isfile(aseg_mgz):
        
        # Running autorecon1
        os.environ["SUBJECTS_DIR"] = fssubj_dir
        fssubj_dir = Path(fssubj_dir)
        
        # Creating the freesurfer directory
        try:
            fssubj_dir.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            print("The directory to store the spatial transformations does not have write access.")
            sys.exit()
        
        if not os.path.isfile(os.path.join(fssubj_dir, full_id, 'mri', 'brainmask.mgz')):
            cmd_bashargs = ['recon-all', '-i', t1, '-subjid', full_id, '-autorecon1']
            cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
            subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
        
        # Running the subcortical segmentation based on Fischl et al. 2002 at Neuron
        cmd_bashargs = ['recon-all', '-subjid', full_id, '-autorecon2-volonly']
        cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
        
        
    # Moving the resulting parcellation from conform space to native
    raw_vol = os.path.join(fssubj_dir, full_id, 'mri', 'rawavg.mgz')

    cmd_bashargs = ['mri_vol2vol', '--mov', aseg_mgz, '--targ', raw_vol,
                    '--regheader', '--o', aseg_nii, '--no-save-reg', '--interp', 'nearest']
    cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
    subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True)
    
    
    # ---------------- Creating Maximum probability Image ------------- #
    # Reading the thalamic parcellation
    spam_Ip          = nib.load(thalparc_spam)
    affine           = spam_Ip.affine
    spam_vol          = spam_Ip.get_fdata()

    spam_vol[spam_vol < prob_thresh] = 0
    spam_vol[spam_vol > 1]      = 1
    
    # Reading the aseg file
    aseg_Ip          = nib.load(aseg_nii)
    aseg_vol         = aseg_Ip.get_fdata()
    
    # 1. Left Hemisphere
    vol2rem_lh = info_dict["rh"]["index"]
    
    # Convert vol2rem_lh to numpy array
    vol2rem_lh = np.array(vol2rem_lh)-1
    
    # Creating the maxprob
    array_data = np.delete(spam_vol, vol2rem_lh, 3)
    ind              = np.where(np.sum(array_data, axis=3) == 0)
    maxprob_thl      = array_data.argmax(axis=3) + 1
    maxprob_thl[ind] = 0
    
    # Selecting the thalamus
    index        = np.where(aseg_vol != 10)
    maxprob_thl[index[0], index[1], index[2]] = 0
    
    # 2. Right Hemisphere
    vol2rem_rh = info_dict["lh"]["index"]
    
    # Convert vol2rem_lh to numpy array
    vol2rem_rh = np.array(vol2rem_rh)-1
    
    # Creating the maxprob
    array_data = np.delete(spam_vol, vol2rem_rh, 3)
    ind              = np.where(np.sum(array_data, axis=3) == 0)
    maxprob_thr      = array_data.argmax(axis=3) + 1
    maxprob_thr[ind] = 0
    
    # Selecting the thalamus
    index        = np.where(aseg_vol != 49)
    maxprob_thr[index[0], index[1], index[2]] = 0
    ind              = np.where(maxprob_thr != 0)
    maxprob_thr[ind] = maxprob_thr[ind] + 7
    
    # Saving the Nifti file
    imgcoll          = nib.Nifti1Image(maxprob_thr.astype('int16') + maxprob_thl.astype('int16'), affine)
    nib.save(imgcoll, thalparc_maxprob)
    
    # Preparting to create the colortable
    st_codes = info_dict["lh"]["index"] + info_dict["rh"]["index"]
    st_names = info_dict["lh"]["name"] + info_dict["rh"]["name"]
    st_colors = info_dict["lh"]["color"] + info_dict["rh"]["color"]
    
    # Saving the lut table
    cltparc.Parcellation.write_luttable(codes=st_codes, 
                        names=st_names, 
                        colors=st_colors,
                        out_file = thalparc_lutfile, force=True)
    
    # Saving the tsv table
    tsv_df = pd.DataFrame({"index": st_codes, "name": st_names, "color": st_colors})
    cltparc.Parcellation.write_tsvtable(tsv_df=tsv_df,
                        out_file = thalparc_tsvfile, force=True)
    # Remove the aseg file
    if os.path.isfile(aseg_nii):
        os.remove(aseg_nii)
    
    return thalparc_maxprob

# simple progress indicator callback function
def progress_indicator(future):
    """
    A simple progress indicator for the concurrent futures
    :param future: future object
    
    """
    global lock, n_subj, n_comp, pb, pb1
    # obtain the lock
    with lock:
        # update the counter
        n_comp += 1
        # report progress
        # print(f'{tasks_completed}/{n_subj} completed, {n_subj-tasks_completed} remain.')
        # pb.update(task_id=pb1, description= f'[red]Completed {n_comp}/{n_subj}', completed=n_subj)
        pb.update(task_id=pb1, description= f'[green]Finished ({n_comp}/{n_subj})', completed=n_comp) 

def main():
    # 0. Handle inputs
    parser = _build_args_parser()
    args = parser.parse_args()

    print(args)
    if args.verbose is not None:
        v = np.int(args.verbose[0])
    else:
        v = 0
        print('- Verbose set to 0\n')
    if v:
        print('\nInputs\n')
    #

    # Getting the path of the current running python file
    bids_dirercts   = args.bidsdir[0].split(sep=',')
    # deriv_dir    = args.derivdir[0]
    deriv_dirercts   = args.derivdir[0].split(sep=',')
    
    if args.t1s is not None:
        t1s2run_file = args.t1s[0]
    else:
        t1s2run_file = ''
    
    # Detecting the number of cores to be used
    ncores = os.cpu_count()
    nthreads = args.nthreads[0]
    if nthreads == 'auto':
        nthreads = ncores
        if nthreads > 4:
            nthreads = nthreads - 4
        else:
            nthreads = 1
    else:
        nthreads     = int(args.nthreads[0])
    
    
    # Avoiding the layout to decrease indexing time
    t1s = []

    for cont, bids_dir in enumerate(bids_dirercts):
        if len(bids_dirercts) == len(deriv_dirercts):
            deriv_dir = deriv_dirercts[cont]
        else:
            deriv_dir = deriv_dirercts[0]
            
        # Declaring global variables
        global layout, pb, pb1, n_subj, n_comp, lock
        
        # Selecting all the T1w images for each BIDS directory
        layout = BIDSLayout(bids_dir, validate=False)
        t1s = layout.get( extension=['nii.gz', 'nii'], suffix='T1w', return_type='filename')

        # Filtering the T1w images to be processed
        if os.path.isfile(t1s2run_file):
            t1s = cltmisc._select_ids_from_file(t1s, t1s2run_file)
        else:
            t12run = t1s2run_file.split(',')
            t1s = [s for s in t1s if any(xs in s for xs in t12run)]

        n_subj = len(t1s)
        
        with Progress() as pb:
        
            # create a lock for the counter
            lock = Lock()

            # Completed subjects
            n_comp = 0
            failed = []
            
            # print("Parcellation: % d"% (p+1), "of % d"% (n_parc))
            pb1 = pb.add_task(f'[green]Processing subject ', total=n_subj)
            if nthreads == 1:
                
                for i, t1 in enumerate(t1s):
                    # ent_dict = layout.parse_file_entities(t1)
                    
                    t1_name = os.path.basename(t1)
                    temp = t1_name.split("_")
                    full_id = '_'.join(temp[:-1])
                    
                    pb.update(task_id=pb1, description= f'[green]Processing: {full_id} ({i+1}/{n_subj})', completed=i+1) 
                    _compute_abased_thal_parc(t1, deriv_dir)
                
            else:
                start_time = time.perf_counter()
                
                # Adjusting the number of threads to the number of subjects
                if n_subj < nthreads:
                    nthreads = n_subj
                    
                # start the thread pool
                with ThreadPoolExecutor(nthreads) as executor:
                    # send in the tasks
                    futures = [executor.submit(_compute_abased_thal_parc, t1s[i],
                    deriv_dir) for i in range(n_subj)]
                    
                    # Just testing
                    # futures = [executor.submit(test, t1s[i]) for i in range(n_subj)]
                    
                    # register the progress indicator callback
                    for future in futures:
                        future.add_done_callback(progress_indicator)
                    # wait for all tasks to complete
            
            pb.update(task_id=pb1, description= f'[green]Finished ({n_subj}/{n_subj})', completed=n_subj) 

if __name__ == "__main__":
    main()