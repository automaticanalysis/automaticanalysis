% Jonathan Peelle - August 2009

% To make sure SPM8 is used for DARTEL scripts, run this AA setup script
% from a Matlab session that has SPM8 in the path (i.e., start a session
% with 'spm 8'). This seemed to do the trick for me; even though an SPM5
% job is still launched, spm8 is the first installation in the path.
%
% Most of the DARTEL AA modules will log the version of the DARTEL command
% they are using, so it's a good idea to see if it's finding the one you
% want it to find in your path. -JP

clear all

% load the sample (subjects in S, id numbers S(i).id)
load /imaging/local/structurals/nca/scripts/jp_sample_59perdecile_2009-8-2_18-30-10.mat


% Add this if not in your path
%addpath /imaging/local/spm/aa_svn_new/devel/
addpath /imaging/jp01/software/aa/devel
aa_ver3_devel


%% Analysis recipe
% You may want to make a copy of the dartelvbm8.xml file in your local
% directory and comment it so that it runs through segmentation but not
% template creation; this gives you a chance to make sure all segmentations
% look appropriate etc. before continuing. Alternatively, of course, you
% can run it all the way through, and run it again if something needs to be
% fixed.
aap = aarecipe('aap_parameters_defaults.xml', 'aap_tasklist_dartelvbm8.xml');

%% Define study-specific parameters

aap.acq_details.root = '/imaging/jp01/examples/SPM8_aa_DARTEL'; % the study directory
aap.options.autoidentifyfieldmaps = 0;
aap.options.autoidentifystructural_chooselast = 1;
aap.options.copystructuraltocentralstore = 0;
aap.options.deletestructuralaftercopyingtocentralstore = 0;
aap.options.aa_minver = 3;

% This is the subdirectory of structurals where the rc* and u_* files live
aap.directory_conventions.dartelsubjdirname='dartel01';

% make the timeout longer (needed probably just for the template stage)
aap.timeouts.busy=24*60*4; % units are minutes, this is 4 days

% Normal segmentation gives us GM/GM; with segment 8 we should be able to
% use all tissue classes (although not entirely sure this matters/helps?)
aap.tasksettings.aamod_dartel_createtemplate.numtissueclasses = 6;

aap.tasksettings.aamod_segment8.samp = 1;              % sampling distance in mm
aap.tasksettings.aamod_dartel_normmnistruct.fwhm = 12; % smoothing for output images


%% NB. This keeps the AA output the 'old' way; otherwise output is saved in
%% a separate directory for each stage. Ask me if questions. -JP
aap.directory_conventions.outputformat = [];


%% Add subjects
for i=1:2
    aap = aas_addsubject(aap, sprintf('%s_*/*', S(i).id));
end



%% Do processing
aap = aa_doprocessing(aap);

