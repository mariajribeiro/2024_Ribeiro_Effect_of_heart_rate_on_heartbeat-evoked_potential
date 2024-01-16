Code for

The heartbeat-evoked potential is modulated by heart rate on a task and mental effort dependent manner
Francesca Aprile, Marco Simões, Jorge Henriques, Paulo Carvalho, Miguel Castelo-Branco, Alejandra Sel, Maria J. Ribeiro

bioRxiv ; doi: https://doi.org/


Task and acquisition parameters are described in detail in
Age-related differences in event-related potentials and pupillary responses in cued reaction time tasks
Maria J. Ribeiro, Miguel Castelo-Branco. Neurobiol Aging. 2019 Jan;73:177-189. doi: 10.1016/j.neurobiolaging.2018.09.028


Analyses scripts run on MATLAB and depend on:

	eeglab functions - https://sccn.ucsd.edu/eeglab/index.php
		in topoplot function change line 275 to COLORARRAY  = { [0 0 0] [0.5 0 0] [0 0 0] }; to make circles black

	limo_eeg - eeglab plugin - https://github.com/LIMO-EEG-Toolbox/limo_meeg
	
	Pupillary waveform deblinking MATLAB code stublinks.m from https://sites.pitt.edu/~gsiegle/

	colour maps from http://www.fabiocrameri.ch/colourmaps.php
	https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

	mintersect - https://www.mathworks.com/matlabcentral/fileexchange/6144-mintersect-multiple-set-intersection

	jbfill - https://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves

	fooof toolbox – https://fooof-tools.github.io/fooof/index.html



Task scripts are available at https://github.com/CIBIT-UC/2021_Ribeiro_BrainVariability_Aging/tree/main/task_code_five_runs
and depend on:

	Psychophysics Toolbox Version 3 - http://psychtoolbox.org/


Pupil preprocessing scripts are available at https://github.com/CIBIT-UC/2021_Ribeiro_BrainVariability_Aging/tree/main/preprocessing_pupil
and depend on stublinks.m available here https://sites.pitt.edu/~gsiegle/


DATASET
Raw EEG, ECG, and pupil data available in OpenNeuro Repository - https://openneuro.org/datasets/ds003690/versions/1.0.0
Maria J. Ribeiro and Miguel Castelo-Branco (2021). EEG, ECG and pupil data from young and older adults: rest and auditory cued reaction time tasks. OpenNeuro. 
[Dataset] doi: 10.18112/openneuro.ds003690.v1.0.0


For questions, feel free to contact by email mjribeiro@fmed.uc.pt; twitter @Ribeiro_neurosc