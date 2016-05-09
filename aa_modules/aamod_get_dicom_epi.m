% This module finds all of the DICOM files associated with this EPI session
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_epi(aap,task,subj,sess)

function [aap resp]=aamod_get_dicom_epi(aap,task,subj,sess)
global aaworker
resp='';

switch task
    case 'description'
        resp=sprintf('Getting EPI DICOM files');
        
    case 'summary'
        resp=sprintf('Getting EPI DICOM files\n');
        
    case 'report'
    case 'doit'
        %% Locate series
        [d, seriesnum] = aas_get_series(aap,'functional',subj,sess);
        
        %% Go through each echo
        out=[];
        for seriesind=1:length(seriesnum)
            [aap dicom_files_src]=aas_listdicomfiles(aap,[subj d],seriesnum(seriesind));

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
                echopath=aas_getsesspath(aap,subj,sess);
            else
                echopath=fullfile(aas_getsesspath(aap,subj,sess),sprintf('echo%d',seriesind));
            end
            aas_makedir(aap,echopath);
            outstream={};
            switch(aap.directory_conventions.remotefilesystem)
                case 'none'
                    for ind=1:length(dicom_files_src)
                        [success,msg,msgid] = copyfile(deblank(...
                            dicom_files_src{ind}),echopath);
                        % allow copyfile failures due to permissions issues
                        % (e.g. if copying from networked system)
                        assert(success || strfind(msg,'chflags'),...
                            'copyfile failed!')
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
                    series_numbers=aap.acq_details.subjects(subj).seriesnumbers(sess);
                    if (~iscell(series_numbers))
                        series_numbers={series_numbers};
                    end
                    for seriesind=1:length(series_numbers)
                        dicomdirsearchpth=fullfile(aap.acq_details.subjects(subj).mriname{1},sprintf(aap.directory_conventions.seriesoutputformat,series_numbers{seriesind}));
                        field_attr=[];
                        field_attr.seriesid.value=dicomdirsearchpth;
                        field_attr.indataset.nid=aap.acq_details.subjects(subj).drupalnid;
                        
                        attr.parent_node=aap.acq_details.subjects(subj).drupalnid;
                        
                        [pth nme ext]=fileparts(aas_getsubjpath(aap,subj));
                        field_attr.seriesalias.value=fullfile(aap.acq_details.sessions(sess).name);
                        [aap waserror drupalnid]=drupal_checkexists(aap,'series',sprintf(aap.directory_conventions.seriesoutputformat,series_numbers{seriesind}),field_attr);
                        aap.acq_details.subjects(subj).seriesnumbers_drupalnid{sess}{seriesind}=drupalnid;
                    end
            end
            out=[out outstream];
        end
        
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
                
                for imnum = 1:numel(out)
                    aas_shell(sprintf('dcmodify -m "(%04x,%04x)=%s" %s',group,element,f{1}.Value,out{imnum}));
                end
            end
        end
        
        %% Output        
        aap=aas_desc_outputs(aap,subj,sess,'dicom_epi',out);
end
