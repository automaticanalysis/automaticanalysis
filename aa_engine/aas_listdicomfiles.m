function [aap,resp]=aas_listdicomfiles(aap,i,seriesnum)
global aaworker
resp={};

dicomdirsearchpth=fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).mriname,sprintf(aap.directory_conventions.seriesoutputformat,seriesnum));

switch(aap.directory_conventions.remotefilesystem)
    case 'none'
        subjpath=aas_getsubjpath(aap,i);
        aisfn=fullfile(subjpath,'autoidentifyseries_saved.mat');
        ais=load(aisfn);
        resp=ais.alldicomfiles{seriesnum};
    case 's3'
        
        % Just get everything with that prefix
        [aap s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,[dicomdirsearchpth '_']);
        if (~isfield(s3resp,'Contents'))
            aas_log(aap,false,sprintf('Can''t find dicom files for series %d in path %s',seriesnum,fullfile(dicomdirsearchpth,aap.directory_conventions.dicomfilter)));
        else
            
            for i=1:length(s3resp.Contents)
                resp=[resp s3resp.Contents(i).Key];
            end
        end
end