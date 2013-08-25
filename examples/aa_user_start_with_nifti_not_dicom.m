% Automatic analysis
% User master script
% Rhodri Cusack BMI LondON Nov 2011

% This is an example of a branched analysis specified from within a user script
aap=aarecipe('aap_tasklist_neonate.xml');

aap=aas_localconfig(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=4.0; % will only work on aa version 4.0 or above

% Location of raw DICOM data
aap.directory_conventions.rawdatadir='/home/rcusack/Neonatal_fMRI/rawdata';

% Where to put the analyzed data 
aap.acq_details.root = '/home/rcusack/Neonatal_fMRI';
aap.directory_conventions.analysisid='analaysis_v4_10images_2slicesperchunk';


% The subjects
aap=aas_addsubject(aap,'NNU982',{'NNU982-fMRIbal.nii'});
%aap=aas_addsubject(aap,'NNU1056',{'NNU1056_fMRIbal.nii'});
%aap=aas_addsubject(aap,'NNU1070',{'NNU1070-fMRIbal.nii'});
%aap=aas_addsubject(aap,'NNU922',{'NNU922_fMRIbal.nii'});
%aap=aas_addsubject(aap,'NNU961',{'NNU961-fMRI_bal.nii'});

aap=aas_addsession(aap,'motor');

% Model and contrasts
%aap=aas_addevent(aap,'aamod_firstlevel_model','*','passive','anysound',-2+ [ 4.0100   10.1300   16.1500   22.3000   28.4300   34.4700   40.5300   46.6600 52.7000   58.9100   64.9400   70.9700   77.0700   83.1300   89.3000   95.3400 101.5100  107.5500  113.7600  120.0100  126.2100  132.2400  138.3600  144.4700],3);
%aap=aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:attentionmemory_block3',[0 1 0 -1],'count_slow_vs_remember_slow');


% DO PROCESSING

aa_doprocessing(aap);


