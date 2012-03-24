% This module finds all of the DICOM files associated with the fieldmaps
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_fieldmap(aap,task,i)

function [aap resp]=aamod_get_dicom_fieldmap(aap,task,i)
global aaworker

resp='';

switch task
    case 'description'
        resp=sprintf('Getting fieldmaps DICOM files');
        
    case 'summary'
        resp=sprintf('Getting fieldmaps DICOM files\n');
        
    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        fieldpath=fullfile(subjpath,aap.directory_conventions.fieldmapsdirname);
        
        fieldfolds = {'rawmag' 'rawphase'};
        
        % Manually specified value for fieldmaps series number over-rides automatically scanned value
        if (~isempty(aap.acq_details.subjects(i).fieldmaps))
            fieldseries=aap.acq_details.subjects(i).fieldmaps;
        else
            % Load up automatically scanned value, validate
            aisfn=fullfile(subjpath,'autoidentifyseries_saved.mat');
            ais=load(aisfn);
            if (length(ais.series_newfieldmap)>2)
                aas_log(aap,false,sprintf('Was expecting only two fieldmaps, but autoidentify series found %d. Will proceed with last t2o, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_newfieldmap)));
                fieldseries=ais.series_newfieldmap(end-1:end);
            elseif (isempty(ais.series_newfieldmap))
                aas_log(aap,true,'No fieldmaps found.');
            else
                fieldseries=ais.series_newfieldmap;
            end;
        end;
        
        dicom_files_src = {};
        % Go through each fieldmap
        out=[];
        for seriesind=1:length(fieldseries)
            [aap dicom_files_src]=aas_listdicomfiles(aap,i,fieldseries(seriesind));
            
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
        
        aap=aas_desc_outputs(aap,i,'dicom_fieldmap',out);
end;
