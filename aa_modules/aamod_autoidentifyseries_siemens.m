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
        
        % Get a listing of all of the series for this subject
        rawdata_subj=aas_findvol(aap,i);
        switch (aap.directory_conventions.remotefilesystem)
            case 'none'
                % Looks for DICOMs in here and in all subdirectories
                allpths=textscan(genpath(rawdata_subj),'%s','delimiter',':'); allpths = allpths{1};
                serieslist=[];
                acqlist=[];
                alldicomfiles={};
                for thispth = allpths'
                    fprintf('Examining folder %s\n',thispth);
                    fn=dir(fullfile(thispth,aap.directory_conventions.dicomfilter));
                    for fnind=1:length(fn)
                        if ~fn(fnind).isdir
                            fullfn=fullfile(thispth,fn(fnind).name);
                            H=aas_dicom_headers_light(fullfn);
                            if (isfield(H{1},'SeriesNumber') && isfield(H{1},'AcquisitionNumber'))
                                serieslist=[serieslist H{1}.SeriesNumber];
                                acqlist=[acqlist H{1}.AcquisitionNumber];
                                if (H{1}.SeriesNumber>length(alldicomfiles))
                                    alldicomfiles{H{1}.SeriesNumber}=[];
                                end
                                alldicomfiles{H{1}.SeriesNumber}{end+1}=fullfn;
                            end
                        end
                    end
                    
                end
                
                %  Check that at least one dicom found, or paths are
                %  probably wrong
                if isempty(serieslist)
                    aas_log(aap,true,sprintf('No DICOM files found in directory\n %s\nusing dicom filter\n %s, quitting.',rawdata_subj, aap.directory_conventions.dicomfilter));
                end;
                
                % Go through all of the series we've found
                rawdata_allseries=unique(serieslist);
                
                for j=1:length(rawdata_allseries);
                    hdr=spm_dicom_headers(alldicomfiles{rawdata_allseries(j)}{1});
                    aas_log(aap,false,sprintf('Series %d (%s) with %d dicom files',rawdata_allseries(j), hdr{1}.ProtocolName, length(alldicomfiles{rawdata_allseries(j)}))); 
                end
            case 's3'
                % Use delimiter to get series names as CommonPrefixes
                [aap s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,rawdata_subj,[],'/');
                rawdata_allseries={s3resp.CommonPrefixes.Prefix};
        end
        
        % Check if there is any preference for this subject
        % structural
        if ~isempty(aap.acq_details.subjects(i).structural)
            structural_choose = aap.acq_details.subjects(i).structural;
        else
            structural_choose = 0;
        end  
        
        series_spgr=[];
        series_fieldmap=[];
        series_tmaps=[];
        series_t2=[];
        
        filenumber=0;
        
        tmpdir=aas_gettempfilename();
        
        protocolnames=[];
        
        % Go through each series, and examine type
        for j=1:length(rawdata_allseries)
            % Get the path to a single dicom file from series "j", downloading from S3 first if necessary
            switch(aap.directory_conventions.remotefilesystem)
                case 's3'
                    seriespth=rawdata_allseries{j};
                    while (seriespth(end)==filesep)
                        seriespth=seriespth(1:(end-1));
                    end
                    [pth nme ext]=fileparts(seriespth);
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
                case 'none'
                    seriesnum=rawdata_allseries(j);
                    dicomfilepath=alldicomfiles{seriesnum}{1};
            end
            
            
            
            
            % For this series, find type from a single DICOM file
            if (~isempty(dicomfilepath))
                hdr=spm_dicom_headers(dicomfilepath);
                
                % Decide whether to ignore this series [djm 20/3/06]
                % absolute number rather than index (jc)
                if ~any(aap.acq_details.subjects(i).ignoreseries == ...
                        rawdata_allseries(j))
                    
                    if (aap.options.autoidentifyfieldmaps)
                        if (findstr(hdr{1}.ProtocolName,aap.directory_conventions.protocol_fieldmap))
                            if numel(series_fieldmap)>aap.options.autoidentifyfieldmaps_number
                                aas_log(aap,true,['ERROR: autoidentifyseries failed - more than a pair of Siemens fieldmap acquisition were found:' sprintf(' %d',series_fieldmap)]);
                            end
                            series_fieldmap=[series_fieldmap seriesnum];
                        end
                    end
                    
                    if (aap.options.autoidentifystructural)
                        if (findstr(hdr{1}.ProtocolName,aap.directory_conventions.protocol_structural))
                            series_spgr=[series_spgr seriesnum];
                        end
                    end
                    
                    if (aap.options.autoidentifyt2)
                        if (findstr(hdr{1}.ProtocolName,aap.directory_conventions.protocol_t2))
                            series_t2=[series_t2 seriesnum];
                        end
                    end
                    
                    if (aap.options.autoidentifytmaps)
                        % Use directory name rather than protocol to
                        % recognise t maps
                        if (findstr(hdr{1}.ProtocolName,'EvaSeries_tTest'))
                            series_tmaps=[series_tmaps seriesnum];
                        end
                    end
%                    fprintf('Protocol %s\n',hdr{1}.ProtocolName);
                end
            end
            
        end
        
        % Save file
        aas_makedir(aap,fileparts(aisfn));
        if exist('alldicomfiles','var')
            save(aisfn,'series_fieldmap','series_spgr','series_tmaps','series_t2','alldicomfiles','rawdata_allseries');
        else
            save(aisfn,'series_fieldmap','series_spgr','series_tmaps','series_t2');
        end
        
        % Make comment
        comment='';
        if aap.options.autoidentifyfieldmaps
            aap.acq_details.subjects(i).fieldmaps=[];
            if numel(series_fieldmap)==aap.options.autoidentifyfieldmaps_number
                comment=[comment '\n  ' sprintf('fieldmap series %d',series_fieldmap(1))];
                % Generalisation of fieldmap number...
                for n = 2:aap.options.autoidentifyfieldmaps_number
                    comment=[comment ' ' sprintf(' and %d',series_fieldmap(n))];
                end
                aap.acq_details.subjects(i).fieldmaps=series_fieldmap;
            else
                aas_log(aap,true,'ERROR: autoidentifyseries failed - one of the fieldmap acquisitions was not found.');
            end
        end        
        
        if (aap.options.autoidentifystructural)
            if length(series_spgr)>1
                if any(structural_choose)
                    series_spgr = intersect(series_spgr,structural_choose);
                    if isempty(series_spgr)
                        aas_log(aap,1,sprintf('ERROR: autoidentifyseries is failed -%s.',sprintf(' structural serie %d not found',structural_choose)));
                    end
                elseif aap.options.autoidentifystructural_chooselast
                    series_spgr = series_spgr(length(series_spgr));
                elseif aap.options.autoidentifystructural_choosefirst
                    series_spgr = series_spgr(1);
                end
            end
            aap.acq_details.subjects(i).structural=series_spgr;
            comment=[comment sprintf('\n  Structural series %d ',series_spgr)];
        end
        
        if (aap.options.autoidentifyt2)
            aap.acq_details.subjects(i).t2=series_t2;
            comment=[comment ['\n  T2 series ' sprintf('%d\t',series_t2)]];
        end
        
        if (aap.options.autoidentifytmaps)
            aap.acq_details.subjects(i).tmaps=series_tmaps;
            comment=[comment ['\n  T-map series ' sprintf('%d\t',series_tmaps)]];
        end
        if ~isempty(comment), aas_log(aap,false,comment); end
        
        aap=aas_desc_outputs(aap,i,'autoidentifyseries',aisfn);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
