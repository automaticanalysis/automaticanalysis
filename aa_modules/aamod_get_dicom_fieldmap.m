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
        [d, mriser] = aas_get_series(aap,strtok(aas_getsesstype(aap),'_'),subj,sess);

        % locate the first pair of fieldmaps after epi
        msg = {{'after' 'first'} ...
            {'before' 'last'}};
        
        if numel(mriser) > 1, mriser = mriser(end); end % last echo
        fieldseries0 = fieldseries0{d};

        fieldseries{1} = fieldseries0(fieldseries0>mriser); % first choice: fieldmap after EPI
        fieldseries{2} = fieldseries0(fieldseries0<mriser); % second choice: fieldmap before EPI
        selection{1} = 1:abs(aap.options.autoidentifyfieldmaps_number);
        selection{2} = numel(fieldseries{2})-(abs(aap.options.autoidentifyfieldmaps_number)-1):numel(fieldseries{2});
        
        if aap.options.autoidentifyfieldmaps_number < 0 % first choice fieldmap before EPI
            fieldseries = fieldseries([2 1]);
            msg = msg([2 1]);
            selection = selection([2 1]);
        end
        
        if numel(fieldseries{1})==abs(aap.options.autoidentifyfieldmaps_number)
            fieldseries=fieldseries{1}(selection{1});
        elseif numel(fieldseries{1})>abs(aap.options.autoidentifyfieldmaps_number)
            aas_log(aap,false,sprintf('INFO:autoidentifyseries found %d fieldmaps %s EPI serie %d',numel(fieldseries{1}),msg{1}{1}, mriser));
            aas_log(aap,false,sprintf('INFO:Will proceed with the %s %d, but you might want to try using the ignoreseries field in aas_addsubject in your user script.',msg{1}{2},abs(aap.options.autoidentifyfieldmaps_number)));
            fieldseries=fieldseries{1}(selection{1});
        elseif numel(fieldseries{1})<abs(aap.options.autoidentifyfieldmaps_number)
            if numel(fieldseries{2})>=abs(aap.options.autoidentifyfieldmaps_number)
                aas_log(aap,false,sprintf('INFO:autoidentifyseries found no(t enough) fieldmaps %s EPI serie %d',msg{1}{1},mriser));
                aas_log(aap,false,sprintf('INFO:Will proceed with the %s %d acquired %s the EPI.',msg{2}{2},abs(aap.options.autoidentifyfieldmaps_number),msg{2}{1}));
                fieldseries=fieldseries{2}(end-(aap.options.autoidentifyfieldmaps_number-1):end);
            else
                aas_log(aap,false,sprintf('ERROR: Was expecting %d fieldmaps %s or %s EPI serie %d, but autoidentifyseries found only %d and %d',...
                    abs(aap.options.autoidentifyfieldmaps_number),msg{1}{1},msg{2}{1},mriser,numel(fieldseries{1}),numel(fieldseries{2})));
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
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj,sess],'dicom_fieldmap',out);
end