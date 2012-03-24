% Automatic analysis cloud client
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Feb 2010

% Script for the special circumstance where different subjects have
% different numbers of sessions. It runs each subject as a separate
% analysis. This precludes "study level", across subject analyses - these
% would have to be done as a separate aa analysis stage.

aa_ver4_devel

aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_typical_fmri.xml');

% Clear all jobs for this user from the queue 
% Don't do this if you have other analyses running!
% Can only clear queue from website unless on AWS at present
%aas_clear_jobq(aap,'rhodri','aws');

sss=load('/imaging/rhodri/camneuro2/asn_test/subjectsstudiesandseries.mat');

global cc

import camneuro.aacc.*;
if (~exist('cc','var') || isempty(cc))
    cc=Aacc();
end;
if (~cc.isSessionAlive())
    cc.login();
end;
% Add subjects sequentially with separate aa_doprocessing calls, as they
% have different numbers of sessions
for subjind=1:1 %length(allsubj)
    
    % DEFINE STUDY SPECIFIC PARAMETERS
    aap.options.aa_minver=4.0; % will only work on aa version 4.0 or above

    % This is an example of a basic analysis 
    aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_typical_fmri.xml');
    
    aap=aas_localconfig(aap);
    
    % Directory for analysed data
    aap.directory_conventions.analysisid='analysis_asn2';
   
    aap.tasksettings.aamod_slicetiming.sliceorder=[1:32];
    
    % The name of your scanner
    aap.directory_conventions.rawdatadir='Machine_MRC-CBU_MRC35119_TrioTim';
    
    % Add subject & sessions
    aap=aas_addsubject(aap,sprintf('Patient_%s_.*/.*/',sss.allsubj{subjind}),sss.seriesnum{subjind});
    for sessind=1:length(sss.seriesnum{subjind})
        aap=aas_addsession(aap,sprintf('session_%d',sessind));
    end;
    
    % SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS
    
    % DO PROCESSING
    aa_doprocessing_camneuro(aap)
    
end;