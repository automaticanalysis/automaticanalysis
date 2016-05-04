% This module finds all of the DICOM files associated with the given 
% specialsession using aas_listdicomfiles and the name of the special
% session, and copies them into the special session directory of this
% module, either across the local filesystem or from s3. Suffixes (i.e. 
% substring after "_") will be treated as subsessions (similarly to the two
% subsessions of the fieldmap) and copied into subfolders. E.g.:
%
%   aap.tasklist.currenttask.outputstreams.stream = {'MT'}
%   aap.acq_details.special_sessions.name = {'MT_baseline' 'MT_MT'}
%   aap.acq_details.subjects.specialseries = [14 15]
%
%   Series 14 and 15 will go to subjpath/MT/baseline and subjpath/MT/MT,
%   respectively.
%
% It then creates the output stream based on the name of the outputstream
% containing all the subsessions in the order as given in the 
% special_sessions.name.

function [aap resp]=aamod_get_dicom_specialseries(aap,task,subj)
global aaworker

resp='';

switch task
    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,subj);
        streamname = strrep(aas_getstreams(aap,'output'),'dicom_',''); streamname = streamname{1};
        sesspath=fullfile(subjpath,streamname);
        
        % Obtain sereies index from the special_session names
        sessind = cell_index({aap.acq_details.special_sessions.name},streamname);
        if ~sessind, aas_log(aap,true,sprintf('Special session not found: %s\n',streamname)); end
        for f = 1:numel(sessind)
            sessfolds{f} = list_index(aap.acq_details.special_sessions(sessind(f)).name,1,2);
        end
        
        % Go through each subsessions
        out=[];
        for seriesind=1:numel(sessind)
            [d, mriser] = aas_get_series(aap,'special',subj,seriesind);
            [aap, dicom_files_src]=aas_listdicomfiles(aap,[subj d],mriser);
            
            % Now copy files to this module's directory
            foldpath = fullfile(sesspath, sessfolds{seriesind});
            aas_makedir(aap,foldpath);
            outstream={};
            switch(aap.directory_conventions.remotefilesystem)
                case 'none'
                    for ind=1:numel(dicom_files_src)
                        copyfile(deblank(dicom_files_src{ind}),foldpath);
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        outstream{ind}=fullfile(foldpath,[nme ext]);
                    end;
                case 's3'
                    s3fles={};
                    for ind=1:length(dicom_files_src)
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        s3fles=[s3fles [nme ext]];
                        outstream{ind}=fullfile(foldpath,s3fles{ind});
                    end;
                    s3_copyfrom_filelist(aap,foldpath,s3fles,aaworker.bucketfordicom,pth);
            end;
            out=[out outstream];
        end;
        
        %% To Edit
        % DICOM dictionary
        dict = load(aas_getsetting(aap,'DICOMdictionary'));
        if isempty(getenv('DCMDICTPATH'))
            setenv('DCMDICTPATH',fullfile(aap.directory_conventions.DCMTKdir,'share','dcmtk','dicom.dic'));
        end
        
        % Fields to edit
        toEditsetting = aas_getsetting(aap,'toEdit');
        toEdit = toEditsetting(strcmp({toEditsetting.subject},aas_getsubjname(aap,subj)));
        
        % do it
        if ~isempty(toEdit)
            for f = {toEdit.DICOMfield}
                group = dict.group(strcmp({dict.values.name}',f{1}.FieldName));
                element = dict.element(strcmp({dict.values.name}',f{1}.FieldName));
                
                for imnum = 1:numel(out)
                    aas_shell(sprintf('dcmodify -m "(%04x,%04x)=%s" %s',group,element,f{1}.Value,out{imnum}));
                end
            end
        end
        
        %% Output
        aap=aas_desc_outputs(aap,subj,['dicom_' streamname],out);
end;
