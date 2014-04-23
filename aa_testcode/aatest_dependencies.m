% aa testing module
% Tests dependecy generation, goes with figure (that should be uploaded to
% wiki)
% Rhodri Cusack, BMI Western, Nov 2013
%
% output should look like this:
% *** possiblestreamlocations ***
% Source is subject, target is session[ 3 4 ]: subject[ 3 ] session[ 3 4 ] searchlight_package[ 3 4 1 ] searchlight_package[ 3 4 2 ] 
% Source is session, target is subject[ 3 ]: session[ 3 1 ] searchlight_package[ 3 1 1 ] searchlight_package[ 3 1 2 ] session[ 3 2 ] searchlight_package[ 3 2 1 ] searchlight_package[ 3 2 2 ] session[ 3 3 ] searchlight_package[ 3 3 1 ] searchlight_package[ 3 3 2 ] diffusion_session[ 3 1 ] diffusion_session_bedpostx[ 3 1 ] 
% Source is hyperalignment_searchlight_package, target is subject[ 3 ]: hyperalignment_searchlight_package[ 1 ] hyperalignment_subject[ 1 1 ] hyperalignment_subject[ 1 2 ] hyperalignment_subject[ 1 3 ] hyperalignment_subject[ 1 4 ] hyperalignment_searchlight_package[ 2 ] hyperalignment_subject[ 2 1 ] hyperalignment_subject[ 2 2 ] hyperalignment_subject[ 2 3 ] hyperalignment_subject[ 2 4 ] 
% Source is diffusion_bedpostx, target is diffusion_session[ 3 1 ]: diffusion_session_bedpostx[ 3 1 ] 
% 
% *** doneflaglocations ***
% Source is subject, target is session[ 3 4 ]: subject[ 3 ] 
% Source is session, target is subject[ 3 ]: session[ 3 1 ] session[ 3 2 ] session[ 3 3 ] 
% Source is hyperalignment_searchlight_package, target is subject[ 3 ]: hyperalignment_searchlight_package[ 1 ] hyperalignment_searchlight_package[ 2 ] 
% Source is diffusion_bedpostx, target is diffusion_session[ 3 1 ]: diffusion_session_bedpostx[ 3 1 ] 


% A recipe to set up aap structure
aap=aarecipe('aap_tasklist_typical_fmri.xml');

% Four simulated subjects
for subjind=1:4
    aap=aas_addsubject(aap,sprintf('subj%d',subjind));
end;

% Three simulated sessions
for sessind=1:3
    aap=aas_addsession(aap,sprintf('sess%d',sessind));
end;

% Hyperalignment searchlight packages
aap.options.searchlight.Npackage=2;

% Now testing. Can't test 'doneflaglocations_thatexist' as need aap
% structure from full running analysis (which includes aap.internal)
% But, dependencies are the same as for 'doneflaglocations'
wtr_list={'possiblestreamlocations','doneflaglocations'};
for wtrind=1:length(wtr_list)
    fprintf('\n*** %s ***\n', wtr_list{wtrind});
    clear global aacache;
    [deps commonind desc]=aas_getdependencies_bydomain(aap,'subject','session',[3 4],wtr_list{wtrind},1);
    fprintf('Source is subject, target is session[ 3 4 ]: %s\n',desc);
    
        [deps commonind desc]=aas_getdependencies_bydomain(aap,'subject','session',[2 2],wtr_list{wtrind},1);
    fprintf('Source is subject, target is session[ 2 2 ]: %s\n',desc);
    
%     clear global aacache;
    [deps commonind desc]=aas_getdependencies_bydomain(aap,'session','subject',[3],wtr_list{wtrind},1);
    fprintf('Source is session, target is subject[ 3 ]: %s\n',desc);
    
%     clear global aacache;
    [deps commonind desc]=aas_getdependencies_bydomain(aap,'hyperalignment_searchlight_package','subject',[3],wtr_list{wtrind},1);
    fprintf('Source is hyperalignment_searchlight_package, target is subject[ 3 ]: %s\n',desc);
    
%     clear global aacache;
    [deps commonind desc]=aas_getdependencies_bydomain(aap,'diffusion_session_bedpostx','diffusion_session',[3 1],wtr_list{wtrind},1);
    fprintf('Source is diffusion_bedpostx, target is diffusion_session[ 3 1 ]: %s\n',desc);
end;