% This module finds all of the DICOM files associated with the structural
% using aas_listdicomfiles, and copies them into the session directory of
% this module, either across the local filesystem or from s3. It then
% creates the output stream.
% function aap=aamod_get_dicom_structural(aap,task,i)

function [aap, resp]=aamod_get_dicom_structural(aap,task,i)
global aaworker

resp='';

switch task
    case 'description'
        resp=sprintf('Getting structural DICOM files');
        
    case 'summary'
        resp=sprintf('Getting structural DICOM files\n');
        
    case 'report'
    case 'doit'
        outstreamname = aas_getstreams(aap,'output'); outstreamname = outstreamname{1};
        modality = strrep(outstreamname,'dicom_','');
        
        subjpath=aas_getsubjpath(aap,i);
        structpath=fullfile(subjpath,aap.directory_conventions.structdirname);
        
        % Make structural directory
        aas_makedir(aap,structpath);
        
        % Load up automatically scanned value, validate
        aisfn=fullfile(subjpath,'autoidentifyseries_saved.mat');
        ais=load(aisfn);
        switch modality
            case 'structural'
                series = ais.series_spgr;
                % Manually specified value for structural series number over-rides automatically scanned value
                structural_choose = '';
                for d = 1:numel(aap.acq_details.subjects(i).mriname)
                    if ~isempty(aap.acq_details.subjects(i).structural{d})
                        structural_choose = sprintf('%s %s - serie(s)%s',structural_choose,...
                            aap.acq_details.subjects(i).mriname{d},...
                            sprintf(' %d',aap.acq_details.subjects(i).structural{d}));
                        series{d}=intersect(series{d},aap.acq_details.subjects(i).structural{d});
                    end
                end
            case 't2'
                series = ais.series_t2;
        end
        
        switch numel(cell2mat(series))
            case 0
                aas_log(aap,true,sprintf('ERROR: %s not found%s.',modality,structural_choose));
            case 1
                structseries = series;
            otherwise
                if aap.options.(['autoidentify' modality '_multiple']) || aap.options.(['autoidentify' modality '_average'])
                    aas_log(aap,false,sprintf('INFO: Was expecting multiple %s and found %d.',modality,numel(cell2mat(series))));
                    structseries=series;
                else
                    ser_sel = [];
                    for d = 1:numel(aap.acq_details.subjects(i).mriname)
                        if ~isempty(series{d})
                            ser = series{d};
                            structseries = cell(1,numel(aap.acq_details.subjects(i).mriname));
                            if aap.options.(['autoidentify' modality '_choosefirst'])
                                ser_sel = d;
                                structseries{d} = ser(1);
                                break;
                            elseif aap.options.(['autoidentify' modality '_chooselast'])
                                ser_sel = d;
                                structseries{d} = ser(end);
                            else
                                ser_sel(end+1) = d;
                                structseries = series;
                            end
                        end
                    end
                end
                
                % Inform
                if aap.options.(['autoidentify' modality '_choosefirst']) || aap.options.(['autoidentify' modality '_chooselast'])
                    aas_log(aap,false,sprintf(['INFO: Was expecting only one %s, but autoidentify series found %d.\n' ...
                        '    Will proceed with %s - serie %d\n' ...
                        '    but you might want to try using the ignoreseries field in aas_addsubject in your user script.'],...
                        modality,...
                        numel(cell2mat(series)),...
                        aap.acq_details.subjects(i).mriname{ser_sel},...
                        structseries{ser_sel}));
                else
                    aas_log(aap,true,sprintf(['ERROR: Was expecting only one %s, but autoidentify series found %d.\n' ...
                        '    but you might want to try using the ignoreseries field in aas_addsubject in your user script.'],...
                        modality,...
                        numel(cell2mat(series))));
                end
        end
        
        % In case there is more than one structural, e.g. MP2RAGE [AVG]
        out = [];
        for d = 1:numel(aap.acq_details.subjects(i).mriname)
            for seriesind=1:length(structseries{d})
                [aap, dicom_files_src]=aas_listdicomfiles(aap,[i d],structseries{d}(seriesind));
                
                % Now copy files to this module's directory
                outstream=cell(1,numel(dicom_files_src));
                switch(aap.directory_conventions.remotefilesystem)
                    case 'none'
                        for ind=1:numel(dicom_files_src)
                            copyfile(deblank(dicom_files_src{ind}),structpath);
                            [pth, nme, ext]=fileparts(dicom_files_src{ind});
                            outstream{ind}=fullfile(structpath,[nme ext]);
                        end
                    case 's3'
                        s3fles=cell(1,numel(dicom_files_src));
                        for ind=1:numel(dicom_files_src)
                            [pth, nme, ext]=fileparts(dicom_files_src{ind});
                            s3fles{ind}=[nme ext];
                            outstream{ind}=fullfile(structpath,s3fles{ind});
                        end
                        s3_copyfrom_filelist(aap,structpath,s3fles,aaworker.bucketfordicom,pth);
                end
                out=[out outstream];
            end
        end
        
        aap=aas_desc_outputs(aap,i,outstreamname,out);
end
