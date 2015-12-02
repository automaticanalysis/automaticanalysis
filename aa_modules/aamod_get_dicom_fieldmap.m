% This module finds all of the DICOM files associated with the fieldmaps
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_fieldmap(aap,task,subj)

function [aap resp]=aamod_get_dicom_fieldmap(aap,task,subj)
global aaworker

resp='';

switch task
    case 'description'
        resp=sprintf('Getting fieldmaps DICOM files');
        
    case 'summary'
        resp=sprintf('Getting fieldmaps DICOM files\n');
        
    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,subj);
        fieldpath=fullfile(subjpath,aap.directory_conventions.fieldmapsdirname);
        
        fieldfolds = {'rawmag' 'rawphase'};
        
        % Manually specified value for fieldmaps series number over-rides automatically scanned value
        if ~isempty(aap.acq_details.subjects(subj).fieldmaps)
            fieldseries=aap.acq_details.subjects(subj).fieldmaps;
        else
            % Load up automatically scanned value, validate
            ais=load(fullfile(subjpath,'autoidentifyseries_saved.mat'));
            % Backward compatibility
            if ~isfield(ais,'series_fieldmap'), ais.series_fieldmap = ais.series_newfieldmap; end
                
            switch numel(ais.series_fieldmap)
                case 2
                    fieldseries=ais.series_fieldmap;
                case 1
                    aas_log(aap,true,'ERROR: Was expecting two fieldmaps, but autoidentify series found only one!');
                case 0
                    aas_log(aap,true,'ERROR: No fieldmaps found.');
                otherwise
                    aas_log(aap,false,sprintf('WARNING: Was expecting only two fieldmaps, but autoidentify series found %d. Will proceed with last two, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',numel(ais.series_fieldmap)));
                    fieldseries=ais.series_fieldmap(end-1:end);
            end
        end;
        
        dicom_files_src = {};
        % Go through each fieldmap
        out=[];
        for seriesind=1:length(fieldseries)
            [aap dicom_files_src]=aas_listdicomfiles(aap,subj,fieldseries(seriesind));
            
            % Now copy files to this module's directory
            aas_makedir(aap,fullfile(fieldpath, fieldfolds{seriesind}));
            outstream={};
            switch(aap.directory_conventions.remotefilesystem)
                case 'none'
                    for ind=1:length(dicom_files_src)
                        copyfile(deblank(dicom_files_src{ind}),fullfile(fieldpath, fieldfolds{seriesind}));
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        outstream{ind}=fullfile(fullfile(fieldpath, fieldfolds{seriesind}),[nme ext]);
                    end;
                case 's3'
                    s3fles={};
                    for ind=1:length(dicom_files_src)
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        s3fles=[s3fles [nme ext]];
                        outstream{ind}=fullfile(fullfile(fieldpath, fieldfolds{seriesind}),s3fles{ind});
                    end;
                    s3_copyfrom_filelist(aap,fullfile(fieldpath, fieldfolds{seriesind}),s3fles,aaworker.bucketfordicom,pth);
            end;
            out=[out outstream];
        end;
        
        aap=aas_desc_outputs(aap,subj,'dicom_fieldmap',out);
end;
