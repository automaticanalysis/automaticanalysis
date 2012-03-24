function [aap resp]=aamod_emeg_icaproj(varargin)
% project out blink (and pulse? and heog?) components
% Adapted from Jason's _example_pathway.m
% Jason recommends running this on epoched data after artefact rejection,
% but can also run on continuous data.
% 
% Danny Mitchell 04/04/08
%

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% add paths for eeglab
%addpath /imaging/jt03/spm_eeglab/eeglab_current/eeglab2007November16_beta/;
addpath /imaging/local/spm_eeglab/
addpath /imaging/local/eeglab

%%
for f=1:length(files);

%% project out blink (and pulse?) (and HEOG?) components...
    % perhaps shouldn't try to project out pulse?
    % See users.fmrib.ox.ac.uk/~rami/fmribplugin
    load(files{f})
    
    S.D=files{f};
    S.newfname=['p' D.fname];
    if ~exist(fullfile(D.path,S.newfname),'file') || settings.Overwrite
        if settings.ProjectOutPulse
            S.ic2rem=[D.ica.class.blink; D.ica.class.heog; D.ica.class.pulse];
            fprintf('\nProjecting out eye movement and pulse components...please wait...\n')
        else
            S.ic2rem=[D.ica.class.blink; D.ica.class.heog];
            fprintf('\nProjecting out eye movement components...please wait...\n')            
        end
        clear D varargin
        spm_eeglab_icaproj(S);
    end

    % NOTE, you can create a dataset consisting of a single IC
    % component of interest projected to the sensors using:
    % ic2proj=7  % for example
    % clear S
    % S.D=oldfname;
    % S.ic2proj=ic2proj;
    % S.newfname=['proj_' sprintf('%d_',ic2proj) D.fname];
    % [proj D] = spm_eeglab_icaproj(S);
end

return