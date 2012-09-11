% This module finds all of the DICOM files associated with the structural
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_structural(aap,task,i)

function [aap resp]=aamod_get_dicom_t1(aap,task,i)
global aaworker

resp='';

switch task
    case 'description'
        resp=sprintf('Getting T1 structural DICOM files');

    case 'summary'
        resp=sprintf('Getting T1 structural DICOM files\n');

    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        structpath=fullfile(subjpath,aap.directory_conventions.structdirname);

        % Manually specified value for structural series number over-rides automatically scanned value
        if (~isempty(aap.acq_details.subjects(i).structural))
            structseries=aap.acq_details.subjects(i).structural;
        else
            % Load up automatically scanned value, validate
            aisfn=fullfile(subjpath,'autoidentifyseries_saved.mat');
            ais=load(aisfn);
            if (length(ais.series_spgr)>1)
                if aap.options.autoidentifystructural_multiple
                    aas_log(aap,false,sprintf('Was expecting multiple T1 structurals and found %d.',length(ais.series_spgr)));
                    structseries=ais.series_spgr;
                elseif aap.options.autoidentifystructural_chooselast
                    aas_log(aap,false,sprintf('Was expecting only one T1 structural, but autoidentify series found %d. Will proceed with last, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_spgr)));
                    structseries=ais.series_spgr(end);
                elseif aap.options.autoidentifystructural_choosefirst
                    aas_log(aap,false,sprintf('Was expecting only one structural, but autoidentify series found %d. Will proceed with first, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_spgr)));
                    structseries = ais.series_spgr(1);
                else
                    aas_log(aap,true,sprintf('Was expecting only one T1 structural, but autoidentify series found %d. You might want to try using the ignoreseries field in aas_addsubject in your user script.',length(ais.series_spgr)));
                end
            elseif (isempty(ais.series_spgr))
                aas_log(aap,true,'No T1 structural found.');
            else
                structseries=ais.series_spgr;
            end
        end

        % Make structural directory
        aas_makedir(aap,structpath);

        % In case there is more than one structural, e.g. MP2RAGE [AVG]
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

        aap=aas_desc_outputs(aap,i,'dicom_t1',out);
end
