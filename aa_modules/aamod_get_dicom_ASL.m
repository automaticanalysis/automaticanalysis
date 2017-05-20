% In case of multiple series specified for the same sessions (e.g. for ASL,
% PWI, CBF), subfolders raw, PWI and CBF will be created and images will be
% copied accordingly (automatically based on DICOM header):
%
%   aap.acq_details.special_sessions(1).name = 'ASL';
%   aap.acq_details.subjects(1).specialseries{1}{1} = [14 16 17]
%   % assuming that 14, 16 and 17 are raw ASL, PWI and CBF, respectively
%
%   Series 14, 16 and 17 will go to subjpath/ASL/raw, subjpath/ASL/PWI and
%   subjpath/ASL/CBF, respectively.
%
% It then creates three output streams 'dicom_ASL_raw', 'dicom_ASL_PWI' and
% 'dicom_ASL_CBF' containing all the subsessions.

function [aap resp]=aamod_get_dicom_ASL(aap,task,subj,sess)
global aaworker

resp='';

switch task
    case 'report'
    case 'doit'
        % Go through each subsessions
        [d, mriser] = aas_get_series(aap,'special',subj,sess);
        for seriesind=1:numel(mriser)
            [aap, dicom_files_src]=aas_listdicomfiles(aap,[subj d],mriser(seriesind));
            
            hdr = spm_dicom_headers(dicom_files_src{1}); hdr = hdr{1};
            if ~isempty(strfind(hdr.ImageType,'ORIGINAL')), subsess = 'raw'; end
            if ~isempty(strfind(hdr.ImageType,'DERIVED')) && isempty(strfind(hdr.ImageType,'RELCBF')), subsess = 'PWI'; end
            if ~isempty(strfind(hdr.ImageType,'DERIVED')) && ~isempty(strfind(hdr.ImageType,'RELCBF')), subsess = 'CBF'; end
            
            % Now copy files to this module's directory
            foldpath = fullfile(aas_getsesspath(aap,subj,sess), subsess);
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
            
            %% To Edit
            % DICOM dictionary
            dict = load(aas_getsetting(aap,'DICOMdictionary'));
            if isempty(getenv('DCMDICTPATH'))
                setenv('DCMDICTPATH',fullfile(aap.directory_conventions.DCMTKdir,'share','dcmtk','dicom.dic'));
            end
            % Add the dcmodify path if it does not exist (necessary to edit dicom headers)
            if ~isempty(strfind(getenv('PATH'), fullfile(aap.directory_conventions.DCMTKdir,'bin')))
                setenv('PATH', [getenv('PATH') ':' fullfile(aap.directory_conventions.DCMTKdir,'bin')]);
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
                    
                    for imnum = 1:numel(outstream)
                        aas_shell(sprintf('dcmodify -m "(%04x,%04x)=%s" %s',group,element,f{1}.Value,outstream{imnum}));
                    end
                end
            end
            
            %% Output
            outputs = aas_getstreams(aap,'output');
            aap=aas_desc_outputs(aap,'special_session',[subj,sess],outputs{cell_index(outputs,subsess)},outstream);
        end
end;
