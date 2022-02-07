from os.path import join as opj
import os
import nipype.interfaces.utility as niu
from nipype.interfaces.fsl import (BET, ExtractROI, FAST, FLIRT,
                                   FNIRT, ApplyWarp, MCFLIRT,
                                   MotionOutliers, BinaryMaths, ImageMeants,
                                   DilateImage)
from nipype.interfaces.fsl.maths import MathsCommand
from nipype.interfaces import fsl
from nipype.interfaces.afni import (Despike, TShift, Volreg,
                                    TProject, MaskTool, Refit,
                                    Fourier, Detrend, Maskave,
                                    MaskTool, TStat, Localstat,
                                    Bandpass, Calc, BlurInMask)
from nipype.interfaces.ants import (Registration, ApplyTransforms)

from nipype.interfaces.io import DataSink, SelectFiles
from nipype.interfaces.utility import IdentityInterface, Function
from nipype import Workflow, Node
from nipype.interfaces.c3 import C3dAffineTool
from dsp_functions import *

"""
Set defaults for FSL outputs - should be nii.gz formatting for interoperability between AFNI and ANTS
"""

fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

experiment_dir = '/path/to/outputs'
output_directory = 'datasink'  # subdirectory name for datasink
working_directory = 'workingdir'  # subdirectory name for working directory (saving intermediate steps)

# list of subject identifiers, input your subject numbers.
subject_list = ['01', '02', '03']  # is your iterable list detailing each subject number

# list of session identifiers
ses_list = ['1', '2']  # is your iterable list identifying each session id number

# use this to specify the formatting of your input file path names - this is the BIDS standard format c. 2019
anatomical_filename = opj('sub-ID{subject_id}', 'ses-{ses}', 'anat', 'sub-ID{subject_id}_ses-{ses}*.nii.gz')
functional_filename = opj('sub-ID{subject_id}', 'ses-{ses}', 'func', 'sub-ID{subject_id}_ses-{ses}*.nii.gz')

input_dir = 'Path/To/BIDS/format/'


def preproc_workflow(input_dir, output_dir, anat_file, func_file, scan_size=477, bet_frac=0.37):
    """
    The preprocessing workflow used in the preparation of the psilocybin vs escitalopram rsFMRI scans.
    Workflows and notes are defined throughout. Inputs are designed to be general and masks/default MNI space is provided

    :param input_dir: The input file directory containing all scans in BIDS format
    :param output_dir: The output file directory
    :param anat_file: The format of the anatomical scan within the input directory
    :param func_file: The format of the functional scan within the input directory
    :param scan_size: The length of the scan by number of images, most 10 minutes scans are around 400-500 depending
    upon scanner defaults and parameters - confirm by looking at your data
    :param bet_frac: brain extraction fractional intensity threshold
    :return: the preprocessing workflow
    """
    preproc = Workflow(name='preproc')
    preproc.base_dir = output_dir

    # Infosource - a function free node to iterate over the list of subject names
    infosource = Node(IdentityInterface(fields=['subject_id', 'ses']), name="infosource")

    infosource.iterables = [('subject_id', subject_list), ('ses', ses_list)]

    # SelectFiles - to grab the data (alternative to DataGrabber)
    templates = {'anat': anat_file, 'func': func_file}  # define the template of each file input

    selectfiles = Node(SelectFiles(templates, base_directory=input_dir), name="selectfiles")

    # Datasink - creates output folder for important outputs
    datasink = Node(DataSink(base_directory=experiment_dir,
                             container=output_dir),
                    name="datasink")

    preproc.connect([(infosource, selectfiles, [('subject_id', 'subject_id'),
                                                ('ses', 'ses')])])

    ''' 
    This is your functional processing workflow, used to trim scans, despike the signal, slice-time correct, 
    and motion correct your data 
    '''

    fproc = Workflow(name='fproc')  # the functional processing workflow

    # ExtractROI - skip dummy scans at the beginning of the recording by removing the first three
    trim = Node(ExtractROI(t_min=3,
                           t_size=scan_size,
                           output_type='NIFTI_GZ'),
                name="trim")

    # 3dDespike - despike
    despike = Node(Despike(outputtype='NIFTI_GZ', args='-NEW'),
                   name="despike")
    fproc.connect([(trim, despike, [('roi_file', 'in_file')])])
    preproc.connect([(selectfiles, fproc, [('func', 'trim.in_file')])])

    # 3dTshift - slice time correction
    slicetime = Node(TShift(outputtype='NIFTI_GZ',
                            tpattern='alt+z2'),
                     name="slicetime")
    fproc.connect([(despike, slicetime, [('out_file', 'in_file')])])

    # 3dVolreg - correct motion and output 1d matrix
    moco = Node(Volreg(outputtype='NIFTI_GZ',
                       interp='Fourier',
                       zpad=4,
                       args='-twopass'),
                name="moco")
    fproc.connect([(slicetime, moco, [('out_file', 'in_file')])])

    moco_bpfdt = Node(MOCObpfdt(), name='moco_bpfdt')  # use the matlab function to correct the motion regressor
    fproc.connect([(moco, moco_bpfdt, [('oned_file', 'in_file')])])

    '''
    This is the co-registration workflow using FSL and ANTs
    '''

    coreg = Workflow(name='coreg')

    # BET - structural data brain extraction
    bet_anat = Node(BET(output_type='NIFTI_GZ',
                        frac=bet_frac,
                        robust=True), name="bet_anat")

    # FSL segmentation process to get WM map
    seg = Node(FAST(bias_iters=6,
                    img_type=1,
                    output_biascorrected=True,
                    output_type='NIFTI_GZ'), name="seg")
    coreg.connect([(bet_anat, seg, [('out_file', 'in_files')])])

    # functional to structural registration
    mean = Node(MCFLIRT(mean_vol=True,
                        output_type='NIFTI_GZ'), name="mean")

    # BBR using linear methods for initial transform fit
    func2struc = Node(FLIRT(cost='bbr',
                            dof=6,
                            output_type='NIFTI_GZ'), name='func2struc')
    coreg.connect([(seg, func2struc, [('restored_image', 'reference')])])
    coreg.connect([(mean, func2struc, [('mean_img', 'in_file')])])
    coreg.connect([(seg, func2struc, [(('tissue_class_files', pickindex, 2), 'wm_seg')])])

    # convert the FSL linear transform into a C3d format for AFNI
    f2s_c3d = Node(C3dAffineTool(itk_transform=True,
                                 fsl2ras=True), name='f2s_c3d')
    coreg.connect([(func2struc, f2s_c3d, [('out_matrix_file', 'transform_file')])])
    coreg.connect([(mean, f2s_c3d, [('mean_img', 'source_file')])])
    coreg.connect([(seg, f2s_c3d, [('restored_image', 'reference_file')])])

    # Functional to structural registration via ANTs non-linear registration
    reg = Node(Registration(fixed_image='default_images/MNI152_T1_2mm_brain.nii.gz',
                            transforms=['Affine', 'SyN'],
                            transform_parameters=[(0.1,), (0.1, 3.0, 0.0)],
                            number_of_iterations=[[1500, 1000, 1000], [100, 70, 50, 20]],
                            dimension=3,
                            write_composite_transform=True,
                            collapse_output_transforms=True,
                            metric=['MI'] + ['CC'],
                            metric_weight=[1] * 2,
                            radius_or_number_of_bins=[32] + [4],
                            convergence_threshold=[1.e-8, 1.e-9],
                            convergence_window_size=[20] + [10],
                            smoothing_sigmas=[[2, 1, 0], [4, 2, 1, 0]],
                            sigma_units=['vox'] * 2,
                            shrink_factors=[[4, 2, 1], [6, 4, 2, 1]],
                            use_histogram_matching=[False] + [True],
                            use_estimate_learning_rate_once=[True, True],
                            output_warped_image=True),
               name='reg')

    coreg.connect([(seg, reg, [('restored_image', 'moving_image')])])  # connect segmentation node to registration node

    merge1 = Node(niu.Merge(2), iterfield=['in2'], name='merge1')  # merge the linear and nonlinear transforms
    coreg.connect([(f2s_c3d, merge1, [('itk_transform', 'in2')])])
    coreg.connect([(reg, merge1, [('composite_transform', 'in1')])])

    # warp the functional images into MNI space using the transforms from FLIRT and SYN
    warp = Node(ApplyTransforms(reference_image='default_images/MNI152_T1_2mm_brain.nii.gz',
                                input_image_type=3), name='warp')
    coreg.connect([(moco, warp, [('out_file', 'input_image')])])
    coreg.connect([(merge1, warp, [('out', 'transforms')])])

    preproc.connect([(selectfiles, coreg, [('anat', 'bet_anat.in_file')])])
    preproc.connect([(fproc, coreg, [('moco.out_file', 'mean.in_file')])])

    '''
    Scrubbing workflow - find the motion outliers, bandpass filter, re-mean the data after bpf
    '''

    scrub = Workflow(name='scrub')

    # Generate the Scrubbing Regressor
    scrub_metrics = Node(MotionOutliers(dummy=4,
                                        out_file='FD_outliers.1D',
                                        metric='fd',
                                        threshold=0.4), name="scrub_metrics")

    # regress out timepoints
    scrub_frames = Node(Bandpass(highpass=0,
                                 lowpass=99999,
                                 outputtype='NIFTI_GZ'), name='scrub_frames')
    scrub.connect([(scrub_metrics, scrub_frames, [('out_file', 'orthogonalize_file')])])
    preproc.connect([(coreg, scrub, [('warp.output_image', 'scrub_frames.in_file')])])
    preproc.connect([(selectfiles, scrub, [('func', 'scrub_metrics.in_file')])])

    # mean image for remeaning after bandpass
    premean = Node(TStat(args='-mean',
                         outputtype='NIFTI_GZ'), name='premean')
    # remean the image
    remean2 = Node(Calc(expr='a+b',
                        outputtype='NIFTI_GZ'), name='remean2')
    scrub.connect([(scrub_frames, remean2, [('out_file', 'in_file_a')])])
    scrub.connect([(premean, remean2, [('out_file', 'in_file_b')])])
    preproc.connect([(coreg, scrub, [('warp.output_image', 'premean.in_file')])])

    '''
    Regressors for final cleaning steps
    '''

    regressors = Workflow(name='regressors')

    # Using registered structural image to create the masks for both WM and CSF
    regbet = Node(BET(robust=True,
                      frac=0.37,
                      output_type='NIFTI_GZ'), name='regbet')

    regseg = Node(FAST(img_type=1,
                       output_type='NIFTI_GZ',
                       no_pve=True,
                       no_bias=True,
                       segments=True), name='regseg')
    regressors.connect([(regbet, regseg, [('out_file', 'in_files')])])
    preproc.connect([(coreg, regressors, [('reg.warped_image', 'regbet.in_file')])])

    '''
    Create a cerebrospinal fluid (CSF) regressor 
    '''

    # subtract subcortical GM from the CSF mask
    subcortgm = Node(BinaryMaths(operation='sub',
                                 operand_file='default_images/subcortical_gm_mask_bin.nii.gz',
                                 output_type='NIFTI_GZ',
                                 args='-bin'), name='subcortgm')
    regressors.connect([(regseg, subcortgm, [(('tissue_class_files', pickindex, 0), 'in_file')])])

    # Fill the mask holes

    fillcsf = Node(MaskTool(fill_holes=True,
                            outputtype='NIFTI_GZ'), name='fillcsf')
    regressors.connect([(subcortgm, fillcsf, [('out_file', 'in_file')])])

    # Erode the mask

    erocsf = Node(MaskTool(outputtype='NIFTI_GZ',
                           dilate_inputs='-1'), name='erocsf')
    regressors.connect([(fillcsf, erocsf, [('out_file', 'in_file')])])

    # Take mean csf signal from functional image
    meancsf = Node(ImageMeants(output_type='NIFTI_GZ'), name='meancsf')
    regressors.connect([(erocsf, meancsf, [('out_file', 'mask')])])
    preproc.connect([(coreg, regressors, [('warp.output_image', 'meancsf.in_file')])])

    bpf_dt_csf = Node(CSFbpfdt(), name='bpf_dt_csf')
    regressors.connect([(meancsf, bpf_dt_csf, [('out_file', 'in_file')])])

    '''
    Creates a local white matter regressor
    '''

    # subtract subcortical gm
    subcortgm2 = Node(BinaryMaths(operation='sub',
                                  operand_file='default_images/subcortical_gm_mask_bin.nii.gz',
                                  output_type='NIFTI_GZ',
                                  args='-bin'), name='subcortgm2')
    regressors.connect([(regseg, subcortgm2, [(('tissue_class_files', pickindex, 2), 'in_file')])])

    # fill mask
    fillwm = Node(MaskTool(fill_holes=True,
                           outputtype='NIFTI_GZ'), name='fillwm')
    regressors.connect([(subcortgm2, fillwm, [('out_file', 'in_file')])])

    # erod mask
    erowm = Node(MaskTool(outputtype='NIFTI_GZ',
                          dilate_inputs='-1'), name='erowm')
    regressors.connect([(fillwm, erowm, [('out_file', 'in_file')])])

    # generate local wm
    localwm = Node(Localstat(neighborhood=('SPHERE', 25),
                             stat='mean',
                             nonmask=True,
                             outputtype='NIFTI_GZ'), name='localwm')
    regressors.connect([(erowm, localwm, [('out_file', 'mask_file')])])
    preproc.connect([(coreg, regressors, [('warp.output_image', 'localwm.in_file')])])

    # bandpass filter the local wm regressor
    localwm_bpf = Node(Fourier(highpass=0.01,
                               lowpass=0.08,
                               args='-retrend',
                               outputtype='NIFTI_GZ'), name='loacwm_bpf')
    regressors.connect([(localwm, localwm_bpf, [('out_file', 'in_file')])])

    # detrend the local wm regressor

    localwm_bpf_dt = Node(Detrend(args='-polort 2',
                                  outputtype='NIFTI_GZ'), name='localwm_bpf_dt')
    regressors.connect([(localwm_bpf, localwm_bpf_dt, [('out_file', 'in_file')])])

    '''
    Clean up your functional image with the regressors you have created above
    '''

    # create a mask for blurring filtering, and detrending

    clean = Workflow(name='clean')

    mask = Node(BET(mask=True,
                    functional=True), name='mask')

    mean_mask = Node(MCFLIRT(mean_vol=True,
                             output_type='NIFTI_GZ'), name="mean_mask")

    dilf = Node(DilateImage(operation='max',
                            output_type='NIFTI_GZ'), name='dilf')
    clean.connect([(mask, dilf, [('mask_file', 'in_file')])])
    preproc.connect([(scrub, clean, [('remean2.out_file', 'mask.in_file')])])

    fill = Node(MaskTool(in_file='default_images/MNI152_T1_2mm_brain.nii.gz',
                         fill_holes=True,
                         outputtype='NIFTI_GZ'), name='fill')

    axb = Node(Calc(expr='a*b',
                    outputtype='NIFTI_GZ'), name='axb')
    clean.connect([(dilf, axb, [('out_file', 'in_file_a')])])
    clean.connect([(fill, axb, [('out_file', 'in_file_b')])])

    bxc = Node(Calc(expr='ispositive(a)*b',
                    outputtype='NIFTI_GZ'), name='bxc')
    clean.connect([(mean_mask, bxc, [('mean_img', 'in_file_a')])])
    clean.connect([(axb, bxc, [('out_file', 'in_file_b')])])
    preproc.connect([(scrub, clean, [('remean2.out_file', 'mean_mask.in_file')])])

    #### BLUR, FOURIER BPF, and DETREND

    blurinmask = Node(BlurInMask(fwhm=6,
                                 outputtype='NIFTI_GZ'), name='blurinmask')
    clean.connect([(bxc, blurinmask, [('out_file', 'mask')])])
    preproc.connect([(scrub, clean, [('remean2.out_file', 'blurinmask.in_file')])])

    fourier = Node(Fourier(highpass=0.01,
                           lowpass=0.08,
                           retrend=True,
                           outputtype='NIFTI_GZ'), name='fourier')
    clean.connect([(blurinmask, fourier, [('out_file', 'in_file')])])

    tstat = Node(TStat(args='-mean',
                       outputtype='NIFTI_GZ'), name='tstat')
    clean.connect([(fourier, tstat, [('out_file', 'in_file')])])

    detrend = Node(Detrend(args='-polort 2',
                           outputtype='NIFTI_GZ'), name='detrend')
    clean.connect([(fourier, detrend, [('out_file', 'in_file')])])

    remean = Node(Calc(expr='a+b',
                       outputtype='NIFTI_GZ'), name='remean')
    clean.connect([(detrend, remean, [('out_file', 'in_file_a')])])
    clean.connect([(tstat, remean, [('out_file', 'in_file_b')])])

    concat = Node(ConcatModel(), name='concat')

    # Removes nuisance regressors via regression function
    clean_rs = Node(Bandpass(highpass=0,
                             lowpass=99999,
                             outputtype='NIFTI_GZ'), name='clean_rs')

    clean.connect([(concat, clean_rs, [('out_file', 'orthogonalize_file')])])

    remean1 = Node(Calc(expr='a+b',
                        outputtype='NIFTI_GZ'), name='remean1')
    clean.connect([(clean_rs, remean1, [('out_file', 'in_file_a')])])
    clean.connect([(tstat, remean1, [('out_file', 'in_file_b')])])

    preproc.connect([(regressors, clean, [('bpf_dt_csf.out_file', 'concat.in_file_a')])])
    preproc.connect([(fproc, clean, [('moco_bpfdt.out_file', 'concat.in_file_b')])])

    preproc.connect([(regressors, clean, [('localwm_bpf_dt.out_file', 'clean_rs.orthogonalize_dset')])])
    clean.connect([(remean, clean_rs, [('out_file', 'in_file')])])

    preproc.write_graph(graph2use='flat', format='png', simple_form=True)
    preproc.write_graph(graph2use='colored', format='png', simple_form=True)

    return preproc


