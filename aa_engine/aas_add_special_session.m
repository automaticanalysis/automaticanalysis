function [aap]=aas_add_special_session(aap,name)
% Add a special session (e.g. ASL, MPM, MTI), which are not explicitly defined, to the analysis. These acquisitions are typically conducted once per subject, hence they are not really multiple sessions (of the same acquisition).
% However, they are independent from each other and therfore can be analyzed in parallel. The function can be called multiple times to add more kinds or more of the same kind; however, the 'session' names must be unique in any 
% case. Example useage can be seen in aa_user_ASL, aa_user_MPM and aa_user_MTI.
%
% FORMAT function aap = aas_add_special_session(aap, name)
% This will add a new special session as specified in the name. Also, the input/output streams of the corresponding processing modules (aamod_get_dicom_specialseries and aamod_convert_specialseries) will be created to match the
% actual series. E.g.: aap=aas_add_special_session(aap,'ASL'); results in streams called dicom_ASL, ASL and ASL_dicom_header.
%
% aap           - aap structure with parameters and tasklist
% name          - name of the session for future reference

% Blank template for a session entry
thissess.name=name;

% And put into acq_details, replacing a single blank entry if it exists
firstsession = numel(aap.acq_details.special_sessions)==1 && isempty(aap.acq_details.special_sessions.name);
if firstsession
    aap.acq_details.special_sessions=thissess;
else
    doAdd = true;
    for iSess = 1:numel(aap.acq_details.special_sessions)
        if strcmp(aap.acq_details.special_sessions(iSess).name,thissess.name)
            doAdd = false;
        end
    end
    if doAdd, aap.acq_details.special_sessions(end+1) = thissess; end
end;

%% Adjust streams
stages = arrayfun(@(x) aas_getstagetag(aap,x), 1:numel(aap.tasklist.main.module),'UniformOutput',false);
for c = cell_index(stages,'aamod_get_dicom_specialseries')'
    if firstsession, aap = aas_renamestream(aap,stages{c},'dicom_specialseries',['dicom_' name],'output');
    else aap = aas_renamestream(aap,stages{c},'append',['dicom_' name],'output'); end
    aas_log(aap,false,['INFO: ' stages{c} ' output stream: ''dicom_' name '''']);
end
for c = cell_index(stages,'aamod_convert_specialseries')'
    if firstsession, aap = aas_renamestream(aap,stages{c},'dicom_specialseries',['dicom_' name],'input');
    else aap = aas_renamestream(aap,stages{c},'append',['dicom_' name],'input'); end
    aas_log(aap,false,['INFO: ' stages{c} ' input stream: ''dicom_' name '''']);
    if firstsession
        aap = aas_renamestream(aap,stages{c},'specialseries',name,'output');
        aap = aas_renamestream(aap,stages{c},'specialseries_dicom_header',[name '_dicom_header'],'output');
    else
        aap = aas_renamestream(aap,stages{c},'append',name,'output');
        aap = aas_renamestream(aap,stages{c},'append',[name '_dicom_header'],'output');
    end
    aas_log(aap,false,['INFO: ' stages{c} ' output stream: ''' name '''']);
    aas_log(aap,false,['INFO: ' stages{c} ' output stream: ''' name '_dicom_header''']);
end
