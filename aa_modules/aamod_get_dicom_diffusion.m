% This module finds all of the DICOM files associated with the diffusion
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_structural(aap,task,subjind,diffsessind)

function [aap resp]=aamod_get_dicom_diffusion(aap,task,subj,sess)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        %% Get
        
        [d, mriser] = aas_get_series(aap,'diffusion',subj,sess);
        
        for ind=1:length(mriser)
            [aap fns]=aas_listdicomfiles(aap,[subj d],mriser(ind));

            if length(mriser)==1
                dsesspth= aas_getpath_bydomain(aap,'diffusion_session',[subj sess]);
            else
                dsesspth= aas_getpath_bydomain(aap,'diffusion_session_phaseencode_direction',[subj sess ind]);
            end;
            
            aap=aas_makedir(aap,dsesspth);
            
            outfns={};
            for fnind=1:length(fns)
                copyfile(fns{fnind},dsesspth);
                [pth nme ext]=fileparts(fns{fnind});
                outfns{end+1}=[nme ext];
            end;
        
            %% To Edit
            % DICOM dictionary
            dict = load(aas_getsetting(aap,'DICOMdictionary'));
            if isempty(getenv('DCMDICTPATH'))
                setenv('DCMDICTPATH',fullfile(aap.directory_conventions.DCMTKdir,'share','dcmtk','dicom.dic'));
            end

            % Fields to edit
            toEditsetting = aas_getsetting(aap,'toEdit');
            toEditsubj = toEditsetting(strcmp({toEditsetting.subject},aas_getsubjname(aap,subj)));
            toEdit = [];
            for s = 1:numel(toEditsubj)
                sessnames = regexp(toEditsubj(s).session,':','split');
                if any(strcmp(sessnames,aap.acq_details.sessions(sess).name)),
                    toEdit = toEditsubj(s);
                    break;
                end
            end

            % do it
            if ~isempty(toEdit)
                for f = {toEdit.DICOMfield}
                    group = dict.group(strcmp({dict.values.name}',f{1}.FieldName));
                    element = dict.element(strcmp({dict.values.name}',f{1}.FieldName));

                    for imnum = 1:numel(outfns)
                        aas_shell(sprintf('%s/dcmodify -m "(%04x,%04x)=%s" %s',fullfile(aap.directory_conventions.DCMTKdir,'bin'),group,element,f{1}.Value,out{imnum}));
                    end
                end
            end

            if length(mriser)==1
                %% Output
                aap=aas_desc_outputs(aap,'diffusion_session',[subj sess],'dicom_diffusion',outfns);
            else
                %% Output
                aap=aas_desc_outputs(aap,'diffusion_session_phaseencode_direction',[subj sess ind],'dicom_diffusion',outfns);
            end;
        end;
end
end

