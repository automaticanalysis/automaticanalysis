% Automatic analysis
% User master script example (aa version 5.*.*)
%
% Tibor Auer, MRC-CBSU
% 01-02-2016

%% INITIALISE
clear

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_fmri_fromnifti.xml');
aap = aas_configforSPM12(aap);


% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% typical value localsingle
aap.options.autoidentifyfieldmaps=0;  							% typical value 1
aap.options.NIFTI4D = 1;										% typical value 0
aap.options.email='All.Knowing@mrc-cbu.cam.ac.uk';
aap.tasksettings.aamod_slicetiming.autodetectSO = 1;
aap.tasksettings.aamod_slicetiming.refslice = 17;
aap.tasksettings.aamod_firstlevel_model.xBF.UNITS = 'secs';     % OPTIONS: 'scans'|'secs' for onsets and durations, typical value 'secs'
aap.tasksettings.aamod_firstlevel_model.includemovementpars = 1;% Include/exclude Moco params in/from DM, typical value 1

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything'; 
aap.directory_conventions.analysisid = 'Nature_Paper'; 

% Add data
% Assumed data storage
% 	All data in /imaging/xy00/aa_nifti/0_data
% 	Subjects have specific folders with format vol_xxxx (e.g. vol_0001)
% 	In each subject specific folder there are four files:
%		xxxx_fMRI.nii	- fMRI in 4D NIFTI
%		xxxx_T1.nii		- structural in 3D NIFTI
%		xxxx_fMRI.nii	- fMRI in 4D NIFTI
%		xxxx_ev1.txt	- event file in a tab-delimited 2-column-format text file, one row for each trial (e.g. 18	2)
%			first column	- onset times in seconds (e.g. 18 above)
%			second column	- duration times in seconds (e.g. 2 above)
aap.directory_conventions.rawdatadir = '/imaging/xy00/aa_nifti/0_data';
aap.directory_conventions.subjectoutputformat = 'vol_%s';
aap.directory_conventions.subject_directory_format = 3;
aap.acq_details.numdummies = 0;
d = dir(fullfile(aap.directory_conventions.rawdatadir,'vol*'));
for s = 1:numel(d)
    subjstr = d(s).name;
    subj = d(s).name(5:end);
    aap = aas_addsubject(aap,subj,...
        'structural',fullfile(subjstr,[subj '_T1.nii']),...
        'functional', ...
            struct('fname',fullfile(subjstr,[subj '_fMRI.nii']), ...
                    'hdr',fullfile(subjstr,[subj '_fMRI.json'])));
    % '_fMRI.json' is a file that stores various acquisition parameters, which are usually extracted from the DICOM.
    % In this example, the same parameters are specified for all measurements.
    % Once the files are created, you can comment out these lines.
    hdr.RepetitionTime =    2;      % TR in s
    hdr.EchoTime =          0.027;  % TE in s
    hdr.SliceTiming =       1:32;   % slice order (you can also specify slice acquisition times, which is more accurate)
    makedicomhdr(hdr,fullfile(subjstr,[subj '_fMRI.json']));
 
    aap = aas_addsession(aap,'run1');
    ev = importdata(fullfile(aap.directory_conventions.rawdatadir,subjstr,[subj '_ev1.txt']));
    aap = aas_addevent(aap,'aamod_firstlevel_model',subjstr,'run1','task',ev(:,1)',ev(:,2)');
    ev = importdata(fullfile(aap.directory_conventions.rawdatadir,subjstr,[subj '_ev2.txt']));
    aap = aas_addevent(aap,'aamod_firstlevel_model',subjstr,'run1','rest',ev(:,1)',ev(:,2)');
end
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:run1',[1 -1],'task-rest','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));