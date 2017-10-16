% This module finds all of the DICOM files associated with the fieldmaps
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_fieldmap(aap,task,subj)

function [aap resp]=aamod_get_dicom_fieldmap(aap,task,subj,sess)
global aaworker

resp='';

switch task
    case 'report'
    case 'doit'
        sesspath=aas_getsesspath(aap,subj,sess);
        fieldpath=fullfile(sesspath,aap.directory_conventions.fieldmapsdirname);
        
        %% Locate series
        % Manually specified value for fieldmaps series number over-rides automatically scanned value
        if ~isempty(cell2mat(aap.acq_details.subjects(subj).fieldmaps))
            fieldseries0=aap.acq_details.subjects(subj).fieldmaps;
        else
            % Load up automatically scanned value
            ais=load(aas_getfiles_bystream(aap,subj,'autoidentifyseries'));
            fieldseries0=ais.series_fieldmap;
        end
        % locate epi
        [d, mriser] = aas_get_series(aap,'functional',subj,sess);

        % locate the first pair of fieldmaps after epi
        if numel(mriser) > 1, mriser = mriser(end); end % last echo
        fieldseries0 = fieldseries0{d};
        fieldseries = fieldseries0(fieldseries0>mriser);
        if numel(fieldseries)>aap.options.autoidentifyfieldmaps_number
            aas_log(aap,false,sprintf('INFO:autoidentifyseries found %d fieldmaps after EPI serie %d',numel(fieldseries),mriser));
            aas_log(aap,false,sprintf('INFO:Will proceed with the first %d, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',aap.options.autoidentifyfieldmaps_number));
            fieldseries=fieldseries(1:aap.options.autoidentifyfieldmaps_number);
        elseif numel(fieldseries)<aap.options.autoidentifyfieldmaps_number
            if numel(fieldseries0)>=aap.options.autoidentifyfieldmaps_number
                aas_log(aap,false,sprintf('INFO:autoidentifyseries found no fieldmaps after EPI serie %d',mriser));
                aas_log(aap,false,sprintf('INFO:Will proceed with the last %d acquired before the EPI.',aap.options.autoidentifyfieldmaps_number));
                fieldseries=fieldseries0(end-(aap.options.autoidentifyfieldmaps_number-1):end);
            else
                aas_log(aap,true,sprintf('ERROR: Was expecting %d fieldmaps after EPI serie %d, but autoidentifyseries found only %d',aap.options.autoidentifyfieldmaps_number,mriser,numel(fieldseries)));
            end
        end
        
        %% Obtain
        % Go through each fieldmap
        out=[];
        for seriesind=1:length(fieldseries)
            [aap, dicom_files_src]=aas_listdicomfiles(aap,[subj d],fieldseries(seriesind));
            newpath = fullfile(fieldpath, sprintf('serie%02d',seriesind));
            
            % Now copy files to this module's directory
            aas_makedir(aap,newpath);
            switch(aap.directory_conventions.remotefilesystem)
                case 'none'
                    for ind=1:length(dicom_files_src)
                        copyfile(deblank(dicom_files_src{ind}),newpath);
                    end;
                case 's3'
                    s3fles=spm_file(dicom_files_src,'filename');
                    s3_copyfrom_filelist(aap,fullfile(fieldpath, sprintf('serie%02d',seriesind)),s3fles,aaworker.bucketfordicom,pth);
            end;
            out=[out spm_file(dicom_files_src,'path',newpath)];
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
                    aas_shell(sprintf('%s/dcmodify -m "(%04x,%04x)=%s" %s',fullfile(aap.directory_conventions.DCMTKdir,'bin'),group,element,f{1}.Value,out{imnum}));
                end
            end
        end
        
        %% Output
        aap=aas_desc_outputs(aap,subj,sess,'dicom_fieldmap',out);
end