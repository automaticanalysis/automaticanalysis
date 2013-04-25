% This module finds all of the DICOM files associated with this EPI session
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_epi(aap,task,i,j)

function [aap resp]=aamod_get_dicom_epi(aap,task,i,j)
global aaworker
resp='';

switch task
    case 'description'
        resp=sprintf('Getting EPI DICOM files');
        
    case 'summary'
        resp=sprintf('Getting EPI DICOM files\n');
        
    case 'report'
    case 'doit'
        % Multi-echo enabled or not?
        if (iscell(aap.acq_details.subjects(i).seriesnumbers))
            seriesnum=aap.acq_details.subjects(i).seriesnumbers{j};
        else
            seriesnum=aap.acq_details.subjects(i).seriesnumbers(j);
        end
        % Go through each echo
        out=[];
        for seriesind=1:length(seriesnum)
            [aap dicom_files_src]=aas_listdicomfiles(aap,i,seriesnum(seriesind));

            % trim number of files to a multiple of "multipleof"
            multipleof=aap.tasklist.currenttask.settings.multipleof;
            nfiles=floor(length(dicom_files_src)/multipleof)*multipleof;
            
            % and trim according to ignoreafter
            nfiles=min(nfiles,aap.tasklist.currenttask.settings.ignoreafter);
            dicom_files_src=dicom_files_src(1:nfiles);
            
            %    Put this in to accelerate testing - but remove before use!!!
            %        fprintf('Truncating to first 10 DICOMs\n');
            %         dicom_files_src=dicom_files_src(1:10);
            
            % Now copy files to this module's directory
            if (length(seriesnum)==1)
                echopath=aas_getsesspath(aap,i,j);
            else
                echopath=fullfile(aas_getsesspath(aap,i,j),sprintf('echo%d',seriesind));
            end
            aas_makedir(aap,echopath);
            outstream={};
            switch(aap.directory_conventions.remotefilesystem)
                case 'none'
                    for ind=1:length(dicom_files_src)
                        copyfile(deblank(dicom_files_src{ind}),echopath);
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        outstream{ind}=fullfile(echopath,[nme ext]);
                    end
                case 's3'
                    s3fles={};
                    for ind=1:length(dicom_files_src)
                        [pth nme ext]=fileparts(dicom_files_src{ind});
                        s3fles=[s3fles [nme ext]];
                        outstream{ind}=fullfile(echopath,s3fles{ind});
                    end
                    s3_copyfrom_filelist(aap,echopath,s3fles,aaworker.bucketfordicom,pth);
                    
                    % Check Drupal registration of series
                    series_numbers=aap.acq_details.subjects(i).seriesnumbers(j);
                    if (~iscell(series_numbers))
                        series_numbers={series_numbers};
                    end
                    for seriesind=1:length(series_numbers)
                        dicomdirsearchpth=fullfile(aap.acq_details.subjects(i).mriname,sprintf(aap.directory_conventions.seriesoutputformat,series_numbers{seriesind}));
                        field_attr=[];
                        field_attr.seriesid.value=dicomdirsearchpth;
                        field_attr.indataset.nid=aap.acq_details.subjects(i).drupalnid;
                        
                        attr.parent_node=aap.acq_details.subjects(i).drupalnid;
                        
                        [pth nme ext]=fileparts(aas_getsubjpath(aap,i));
                        field_attr.seriesalias.value=fullfile(aap.acq_details.sessions(j).name);
                        [aap waserror drupalnid]=drupal_checkexists(aap,'series',sprintf(aap.directory_conventions.seriesoutputformat,series_numbers{seriesind}),field_attr);
                        aap.acq_details.subjects(i).seriesnumbers_drupalnid{j}{seriesind}=drupalnid;
                    end
            end
            out=[out outstream];
        end
        aap=aas_desc_outputs(aap,i,j,'dicom_epi',out);
end