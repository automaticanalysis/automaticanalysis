% AA initialisation module - Identify dicom headers from Tim Trio
% Rhodri Cusack MRC CBU Cambridge Dec 2005

function  [aap,resp]=aamod_autoidentifyseries_siemens(aap,task,i)

resp='';

switch task
    case 'report'
        respc={};
        if (aap.options.autoidentifyfieldmaps)
            respc{end+1}='fieldmaps';
        end
        if (aap.options.autoidentifystructural)
            respc{end+1}='structural';
        end
        if (aap.options.autoidentifyt2)
            respc{end+1}='t2';
        end
        if (aap.options.autoidentifytmaps)
            respc{end+1}='realtime t-maps';
        end
        
        if isempty(respc)
            resp='No automatic identification done.';
        else
            resp='Automatically identified ';
            for r = 1:numel(respc)-1
                resp = [resp respc{r} ', '];
            end
            resp = [resp(1:end-2) ' and ' respc{end} '.\n'];
        end
        
    case 'doit'
        global aaworker

        aisfn=fullfile(aas_getsubjpath(aap,i),'autoidentifyseries_saved.mat');
        
        %% Get a listing of all of the series for this subject
        alldicomfiles = cell(1,numel(aap.acq_details.subjects(i).mriname));
        rawdata_allseries = cell(1,numel(aap.acq_details.subjects(i).mriname));
        for d = 1:numel(aap.acq_details.subjects(i).mriname)
            rawdata_subj=aas_findvol(aap,[i d]);
            switch aap.directory_conventions.remotefilesystem
                case 'none'
                    % Looks for DICOMs in here and in all subdirectories
                    allpths=textscan(genpath(rawdata_subj),'%s','delimiter',':'); allpths = allpths{1};
                    alldicomfiles{d} = {};
                    serieslist=[];
                    acqlist=[];
                    for thispth = allpths'
                        fn=dir(fullfile(thispth{1},aap.directory_conventions.dicomfilter));
                        fn = fn(~[fn.isdir]);
                        if isempty(fn), continue; end
                        aas_log(aap,false,sprintf('Examining folder %s',thispth{1}));
                        fullfn = spm_file(char({fn.name}'),'path',thispth{1});
                        H = spm_dicom_headers(fullfn);
                        for fnind=1:numel(H)
                            if (isfield(H{fnind},'SeriesNumber') && isfield(H{fnind},'AcquisitionNumber'))
                                serieslist=[serieslist H{fnind}.SeriesNumber];
                                acqlist=[acqlist H{fnind}.AcquisitionNumber];
                                if (H{fnind}.SeriesNumber>length(alldicomfiles{d}))
                                    alldicomfiles{d}{H{fnind}.SeriesNumber}=[];
                                end
                                alldicomfiles{d}{H{fnind}.SeriesNumber}{end+1}=deblank(fullfn(fnind,:));
                            end
                        end
                    end
                    
                    %  Check that at least one dicom found, or paths are
                    %  probably wrong
                    if isempty(serieslist)
                        aas_log(aap,true,sprintf('No DICOM files found in directory\n %s\nusing dicom filter\n %s, quitting.',rawdata_subj, aap.directory_conventions.dicomfilter));
                    end;
                    
                    % Go through all of the series we've found
                    rawdata_allseries{d}=unique(serieslist);
                    
                    for j=1:length(rawdata_allseries{d});
                        hdr=spm_dicom_headers(alldicomfiles{d}{rawdata_allseries{d}(j)}{1});
                        aas_log(aap,false,sprintf('Series %d (%s) with %d dicom files',rawdata_allseries{d}(j), hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field), length(alldicomfiles{d}{rawdata_allseries{d}(j)})));
                    end
                case 's3'
                    % Use delimiter to get series names as CommonPrefixes
                    [aap, s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,rawdata_subj,[],'/');
                    rawdata_allseries{d}={s3resp.CommonPrefixes.Prefix};
            end
        end
        
        %% Prepeare outputs
        series_spgr=cell(1,numel(aap.acq_details.subjects(i).mriname));
        series_fieldmap=cell(1,numel(aap.acq_details.subjects(i).mriname));
        series_tmaps=cell(1,numel(aap.acq_details.subjects(i).mriname));
        series_t2=cell(1,numel(aap.acq_details.subjects(i).mriname));
        
        tmpdir=aas_gettempfilename();
        
        %% Go through each series, and examine type
        for d = 1:numel(aap.acq_details.subjects(i).mriname)
            
            aas_log(aap,false,sprintf('Examining mridata %s',aap.acq_details.subjects(i).mriname{d}));
            for j=1:length(rawdata_allseries{d})
                % Get the path to a single dicom file from series "j", downloading from S3 first if necessary
                switch aap.directory_conventions.remotefilesystem
                    case 'none'
                        seriesnum=rawdata_allseries{d}(j);
                        dicomfilepath=alldicomfiles{d}{seriesnum}{1};
                    case 's3'
                        seriespth=rawdata_allseries{d}{j};
                        while (seriespth(end)==filesep)
                            seriespth=seriespth(1:(end-1));
                        end
                        [junk, nme, ext]=fileparts(seriespth);
                        dicomseriesname=[nme ext];
                        searchpath=fullfile(rawdata_subj,dicomseriesname);
                        aas_log(aap,false,sprintf('Checking here on s3 %s',searchpath));
                        [aap, s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,searchpath,[],[],1);
                        
                        [keypth, nme, ext]=fileparts(s3resp.Contents(1).Key);
                        dicomfilename=[nme ext];
                        s3_copyfrom_filelist(aap,tmpdir,dicomfilename,aaworker.bucketfordicom,keypth);
                        dicomfilepath=fullfile(tmpdir,dicomfilename);
                        % Option to specify series numbers in terms of numbering of
                        % scanner, or ordering of files
                        
                        if (aap.directory_conventions.rawseries_usefileorder)
                            seriesnum=j;
                        else
                            [aap, seriesnum]=aas_getseriesnumber(aap,dicomseriesname);
                        end
                end
                
                % For this series, find type from a single DICOM file
                if ~isempty(dicomfilepath)
                    hdr=spm_dicom_headers(dicomfilepath);

                    % Decide whether to ignore this series [djm 20/3/06]
                    % absolute number rather than index (jc)
                    if ~any(aap.acq_details.subjects(i).ignoreseries{d} == ...
                            rawdata_allseries{d}(j))
                        
                        if (aap.options.autoidentifyfieldmaps)
                            if (strfind(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),aap.directory_conventions.protocol_fieldmap))
                                if numel(series_fieldmap{d})>aap.options.autoidentifyfieldmaps_number
                                    aas_log(aap,true,['ERROR: autoidentifyseries failed - more than a pair of Siemens fieldmap acquisition were found:' sprintf(' %d',series_fieldmap{d})]);
                                end
                                series_fieldmap{d}=[series_fieldmap{d} seriesnum];
                            end
                        end
                        
                        if (aap.options.autoidentifystructural)
                            if (strfind(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),aap.directory_conventions.protocol_structural))
                                series_spgr{d}=[series_spgr{d} seriesnum];
                            end
                        end
                        
                        if (aap.options.autoidentifyt2)
                            if (strfind(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),aap.directory_conventions.protocol_t2))
                                series_t2{d}=[series_t2{d} seriesnum];
                            end
                        end
                        
                        if (aap.options.autoidentifytmaps)
                            % Use directory name rather than protocol to
                            % recognise t maps
                            if (strfind(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),'EvaSeries_tTest'))
                                series_tmaps{d}=[series_tmaps{d} seriesnum];
                            end
                        end
                    end
                end
                
            end
        end
        
        %% Save output
        aas_makedir(aap,fileparts(aisfn));
        if exist('alldicomfiles','var')
            save(aisfn,'series_fieldmap','series_spgr','series_tmaps','series_t2','alldicomfiles','rawdata_allseries');
        else
            save(aisfn,'series_fieldmap','series_spgr','series_tmaps','series_t2');
        end
        aap=aas_desc_outputs(aap,i,'autoidentifyseries',aisfn);

        %% Check and Select
        comment='\n';
        
        for d = 1:numel(aap.acq_details.subjects(i).mriname)
            comment=[comment sprintf('Checking mridata %s\n',aap.acq_details.subjects(i).mriname{d})];
            
            if aap.options.autoidentifyfieldmaps
                aap.acq_details.subjects(i).fieldmaps{d}=[];
                if numel(series_fieldmap{d})==aap.options.autoidentifyfieldmaps_number
                    comment=[comment '  fieldmap series'];
                    % Generalisation of fieldmap number...
                    for n = 1:aap.options.autoidentifyfieldmaps_number
                        comment=[comment sprintf(' %d',series_fieldmap{d}(n))];
                    end
                    comment=[comment '\n'];
                    aap.acq_details.subjects(i).fieldmaps{d}=series_fieldmap{d};
                elseif isempty(aap.acq_details.subjects(i).seriesnumbers{d}) || isempty(cell2mat(aap.acq_details.subjects(i).seriesnumbers{d}))
                    comment=[comment '  no fieldmap required\n'];
                else
                    aas_log(aap,true,'ERROR: autoidentifyseries failed - one of the fieldmap acquisitions was not found.');
                end
            end
            
            if (aap.options.autoidentifytmaps)
                aap.acq_details.subjects(i).tmaps{d}=series_tmaps{d};
                comment=[comment ['  T-map series' sprintf(' %d',series_tmaps{d})]];
                comment=[comment '\n'];
            end
            
            if aap.options.autoidentifystructural && (numel(series_spgr{d})>1) && ~isempty(aap.acq_details.subjects(i).structural{d})
                series_spgr{d} = intersect(series_spgr{d},aap.acq_details.subjects(i).structural{d});
            end
        end
        
        if (aap.options.autoidentifystructural)
            switch numel(cell2mat(series_spgr))
                case 0
                    aas_log(aap,1,'ERROR: autoidentifyseries is failed - structural series not found');
                case 1
                    ser_sel = [];
                    for d = 1:numel(aap.acq_details.subjects(i).mriname)
                        if ~isempty(series_spgr{d})
                            ser_sel = d;
                            break;
                        end
                    end
                otherwise
                    ser_sel = [];
                    for d = 1:numel(aap.acq_details.subjects(i).mriname)
                        if ~isempty(series_spgr{d})
                            ser = series_spgr{d};
                            series_spgr = cell(1,numel(aap.acq_details.subjects(i).mriname));
                            if aap.options.autoidentifystructural_choosefirst
                                ser_sel = d;
                                series_spgr{d} = ser(1);
                                break;
                            elseif aap.options.autoidentifystructural_chooselast
                                ser_sel = d;
                                series_spgr{d} = ser(end);
                            else
                                ser_sel(end+1) = d;
                                series_spgr{d} = ser;
                            end
                        end
                    end
            end
            aap.acq_details.subjects(i).structural=series_spgr;
            for s = ser_sel
                comment=[comment ['  structural series' sprintf(' %d',series_spgr{s})]];
                comment=[comment '\n'];
            end
        end
        
        if (aap.options.autoidentifyt2)
            switch numel(cell2mat(series_t2'))
                case 0
                    aas_log(aap,1,'ERROR: autoidentifyseries is failed - t2 series not found');
                case 1
                    ser_sel = [];
                    for d = 1:numel(aap.acq_details.subjects(i).mriname)
                        if ~isempty(series_t2{d})
                            ser_sel = d;
                            break;
                        end
                    end
                otherwise
                    ser_sel = [];
                    for d = 1:numel(aap.acq_details.subjects(i).mriname)
                        if ~isempty(series_t2{d})
                            ser = series_t2{d};
                            series_t2 = cell(1,numel(aap.acq_details.subjects(i).mriname));
                            if aap.options.autoidentifyt2_choosefirst
                                ser_sel = d;
                                series_t2{d} = ser(1);
                                break;
                            elseif aap.options.autoidentifyt2_chooselast
                                ser_sel = d;
                                series_t2{d} = ser(end);
                            else
                                ser_sel(end+1) = d;
                                series_t2{d} = ser;
                            end
                        end
                    end
            end
            for s = ser_sel
                comment=[comment ['  T2 series' sprintf(' %d',series_t2{s})]];
                comment=[comment '\n'];
            end
        end        
        if ~isempty(comment), aas_log(aap,false,comment); end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
