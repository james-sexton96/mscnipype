# resting state fMRI preprocessing with Nipype
This repository hosts the [Nipype](https://nipype.readthedocs.io/en/latest/) (Gorgolewski et al., 2011) scripted pipeline used to preprocess resting-state fMRI data for the second trial of psilocybin for Major Depressive Disorder from the [Imperial College London Centre for Psychedelic Research](https://www.imperial.ac.uk/psychedelic-research-centre/). This pipeline incorporates the following packages:

[FMRIB Software Library (FSL)](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) (Smith et al., 2004; Woolrich et al., 2009)

[AFNI](https://afni.nimh.nih.gov/) (Cox, 1996; Cox & Hyde, 1997; Gold et al., 1998)

[Advanced Normalization Tools (ANTS)](http://stnava.github.io/ANTs/) (Avants et al., 2014)

## Details

* rsfmri_preproc.py defines the primary pipeline for preprocessing this data

* dsp_functions.py defines MATLAB signal processing steps to pass sever

* the preprocessing_notebook.ipynb is used to import these and run the workflow

My other work and publications can be found here:
[Google Scholar](https://scholar.google.com/citations?user=3Z64U10AAAAJ&hl=en)

## Tutorial 
First install the necessary packages from [Nipype](https://nipype.readthedocs.io/en/latest/users/install.html). Docker is the easiest way to ensure you have the correct UNIX based tools
and Nipype version.

1) If your scans are stored in a BIDS formatted directory, this pipeline should work with no modifications (only specifying your paths).
However, if your data is not stored in BIDS format, please reformat the 'anatomical' and 'functional_filename' variables in the jupyter notebook 
to fit the formatting of your filenames. 
   
The following codeblock defines this path if it's the first subject's first scan: ' /sub-ID01/ses-01/anat/sub-ID01_ses-01.nii.gz '

```
anatomical_filename = opj('sub-ID{subject_id}', 'ses-{ses}', 'anat', 'sub-ID{subject_id}_ses-{ses}*.nii.gz')
```

2) The subject and ses lists define in the notebook can also be modified depending on your number of subjects and sessions. 
These numbers will be substituted into the filename {subject_id} and {ses} pathname locations. 
   
```
subject_list = ['01', '02', '03']  # I have three subjects I wish to iterate through
ses_list = ['1', '2']  # Each subject has two scans I wish to iterate through
```

3) Once you specify your inputs and outputs, you can pass these variables to the workflow function. This will write your 
graphs for each sub-workflow and the total workflow in both simple and detailed format
   
```
pipeline = preproc_workflow(input_dir=input_dir,
                             output_dir=experiment_dir,
                             anat_file=anatomical_filename,
                             func_file=functional_filename)
```

![WorkflowGraph](https://github.com/james-sexton96/rsfmri_preproc/blob/master/graphs/graph.png?raw=true)