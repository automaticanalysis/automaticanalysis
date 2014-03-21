% This module finds all of the DICOM files associated with the structural
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_structural(aap,task,i)

function [aap resp]=aamod_get_dicom_t2(aap,task,i)

global aaworker

resp='';

switch task
        
    case 'description'
        resp=sprintf('Getting T2 DICOM files');
        
    case 'summary'
        resp=sprintf('Getting T2 DICOM files\n');
        
    case 'report'
        
    case 'doit'
        
        subjpath=aas_getsubjpath(aap,i);
        structpath=fullfile(subjpath,aap.directory_conventions.structdirname);
        
        % Manually specified value for structural series number over-rides automatically scanned value
        %if (~isempty(aap.acq_details.subjects(i).t2))
        %    structseries=aap.acq_details.subjects(i).t2;
        %else
            
        % Load up automatically scanned value, validate
        aisfn=fullfile(subjpath,'autoidentifyseries_saved.mat');
        ais=load(aisfn);
        
        if (length(ais.series_t2)>1)
           
           if aap.options.autoidentifyt2_multiple, 
            aas_log(aap,true,sprintf('Was expecting only one T2, but autoidentify series found %d. You might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_t2))); 
            structseries=ais.series_t2; 
           elseif aap.options.autoidentifyt2_chooselast
            aas_log(aap,false,sprintf('Was expecting only one T2, but autoidentify series found %d. Will proceed with last, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_t2)));
            structseries=ais.series_t2(end);
           elseif aap.options.autoidentifyt2_choosefirst
            aas_log(aap,false,sprintf('Was expecting only one T2 but autoidentify series found %d. Will proceed with first, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_t2)));
            structseries = ais.series_t2(1);
           else
            aas_log(aap,true,sprintf('Was expecting only one T2, but autoidentify series found %d. You might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_t2)));
           end
           
        elseif (isempty(ais.series_t2))
           aas_log(aap,true,'No T2 found.');
        else
           structseries=ais.series_t2;
        end
               
        % Make structural directory
        aas_makedir(aap,structpath);
            
        % In case there is more than one T2
        out = [];
        for seriesind=1:length(structseries)
            [aap dicom_files_src]=aas_listdicomfiles(aap,i,structseries(seriesind));
            
            % Now copy files to this module's directory
            outstream=[];
            switch(aap.directory_conventions.remotefilesystem)
                case 'none'
                    for ind=1:length(dicom_files_src)
                        copyfile(deblank(dicom_files_src{ind}),structpath);
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        outstream{ind}=fullfile(structpath,[nme ext]);
                    end
                case 's3'
                    s3fles={};
                    for ind=1:length(dicom_files_src)
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        s3fles=[s3fles [nme ext]];
                        outstream{ind}=fullfile(structpath,s3fles{ind});
                    end
                    s3_copyfrom_filelist(aap,structpath,s3fles,aaworker.bucketfordicom,pth);
            end
            out=[out outstream];
        end
        
        aap=aas_desc_outputs(aap,i,'dicom_t2',out);
        
case 'checkrequirements'
        
otherwise
        
    aas_log(aap,1,sprintf('Unknown task %s',task));
    
end
    
end