function aatest_ds000114_fmri(deleteprevious,wheretoprocess)
% This script runs a basic modeling pipeline (preprocessing through secondlevel)
% using ds000114 motor and line bisection tasks. The data should be dowloaded 
% prior to running the script

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'),deleteprevious);

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.options.autoidentifystructural_choosefirst = 1;
aap.options.autoidentifystructural_chooselast = 0;

aap.options.NIFTI4D = 1;
aap.options.wheretoprocess = wheretoprocess;
aap.acq_details.numdummies = 4;	
aap.acq_details.numdummies = 1;
aap.acq_details.input.correctEVfordummies = 1;

aap.tasksettings.aamod_segment8.writenormimg = 0;
aap.tasksettings.aamod_segment8.samp = 2;
aap.tasksettings.aamod_slicetiming.autodetectSO = 1;
aap.tasksettings.aamod_slicetiming.refslice = 16;
aap.tasksettings.aamod_smooth.FWHM = 5;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

aap.acq_details.input.combinemultiple = true;

% five subjects is about minimum to run a second level model

aap = aas_processBIDS(aap,[],{'finger_foot_lips','line_bisection'},{'sub-01','sub-02','sub-03','sub-04'});

% -------------------------------------------------------------------------
% modeling - contrast specification
% -------------------------------------------------------------------------

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[1 0 0],'All:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[0 1 0],'All:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[0 0 1],'All:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sessions:+1xfinger_foot_lips_test|-1xfinger_foot_lips_retest',[1 0 0],'T-RT:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sessions:+1xfinger_foot_lips_test|-1xfinger_foot_lips_retest',[0 1 0],'T-RT:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sessions:+1xfinger_foot_lips_test|-1xfinger_foot_lips_retest',[0 0 1],'T-RT:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_test',[1 0 0],'Loc_T:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_test',[0 1 0],'Loc_T:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_test',[0 0 1],'Loc_T:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_retest',[1 0 0],'Loc_RT:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_retest',[0 1 0],'Loc_T:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_retest',[0 0 1],'Loc_T:Lips','T');

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sameforallsessions',[1 0 0 -1 0],'All:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sameforallsessions',[0 0 -1 0 1],'All:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sessions:+1xline_bisection_test|-1xline_bisection_retest',[1 0 0 -1 0],'T-RT:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sessions:+1xline_bisection_test|-1xline_bisection_retest',[0 0 -1 0 1],'T-RT:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_test',[1 0 0 -1 0],'LB_T:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_test',[0 0 -1 0 1],'LB_T:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_retest',[1 0 0 -1 0],'LB_RT:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_retest',[0 0 -1 0 1],'LB_RT:Control:Resp-NoResp','T');

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);


