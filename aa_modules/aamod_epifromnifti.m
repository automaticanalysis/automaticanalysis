% AA module - EPI from NIFTI

function [aap,resp]=aamod_epifromnifti(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        if (not(iscell(aap.acq_details.subjects(subj).seriesnumbers)))
            aas_log(aap,true,'Was exepcting list of filenames in cell array instead of series numbers, check aas_addsubject command in user script');
        end;
        rejectimages=aap.tasklist.currenttask.settings.rejectimages;
        if ischar(rejectimages)
            rejectimages= str2num(rejectimages);
        end;
        for sess=1:length(aap.acq_details.subjects(subj).seriesnumbers)
            finalepis={};
            if (iscell(aap.acq_details.subjects(subj).seriesnumbers{sess}))
                % lots of files 3D
                % [AVG] we cannot cocatenate the single root with
                % multiple 3D image files, so this expects already the
                % full location of the 3D images instead...
                DICOMHEADERS{1}.RepetitionTime = [];
                DICOMHEADERS{1}.EchoTime = [];
				imageFns = aap.acq_details.subjects(subj).seriesnumbers{sess};
                
                sesspth=aas_getsesspath(aap,subj,sess);
                
                if ~exist(sesspth, 'dir')
                    mkdir(sesspth);
                end
                
                for f = 1:length(imageFns)
                    % [AVG] Expects a cell array of images at the moment
                    [root, fn, ext] = fileparts(imageFns{f});
                    % [AVG] Copy file to module location
                    unix(['cp ' imageFns{f} ' ' fullfile(sesspth, [fn ext])]);
                    % [AVG] Add file to what will be described as output...
                    finalepis = [finalepis fullfile(sesspth, [fn ext])];
                end
                V0 = spm_vol(finalepis);
            else
                %Only one file, assume 4D
                DICOMHEADERS{1}.RepetitionTime = [];
                DICOMHEADERS{1}.EchoTime = [];
                niftifile = '';
                niftisearchpth=aas_findvol(aap,'');
                if ~isempty(niftisearchpth)
                    niftifile = fullfile(niftisearchpth,aap.acq_details.subjects(subj).seriesnumbers{sess});
                end
                V = spm_vol(niftifile);
                V0 = V;
                sesspth=aas_getsesspath(aap,subj,sess);
                aas_makedir(aap,sesspth);
                [pth fle ext]=fileparts(aap.acq_details.subjects(subj).seriesnumbers{sess});
                for fileind=1:length(V)
					%Search for header information using space-separated format
                    if exist(fullfile(pth,[fle '.dcmhdr']),'file')
                        fid = fopen(fullfile(pth,[fle '.dcmhdr']),'r');
                        flines = {}; while ~feof(fid), flines{end+1} = fgetl(fid); end
                        fclose(fid);
                        inds = strfind(flines,'Repetition Time');
                        for i = 1:numel(inds)
                            if ~isempty(inds{i}), break; end
                        end
                        str = textscan(flines{i}(inds{i}:end),'%s');
                        DICOMHEADERS{fileind}.RepetitionTime = str2double(str{1}{end});
                        inds = strfind(flines,'Echo Time');
                        for i = 1:numel(inds)
                            if ~isempty(inds{i}), break; end
                        end
                        str = textscan(flines{i}(inds{i}:end),'%s');
                        DICOMHEADERS{fileind}.EchoTime = str2double(str{1}{end});
                    end
                    Y=spm_read_vols(V(fileind));
                    fn=fullfile(sesspth,[sprintf('%s_%04d',fle,fileind) ext]);
                    
                    V(fileind).fname=fn;
                    V(fileind).n=[1 1];
                    spm_write_vol(V(fileind),Y);
                    if (not(any(rejectimages==fileind)))
                        finalepis=[finalepis fn];
                    end;
                end;
            end;
            % Write out the files

            % Now move dummy scans to dummy_scans directory
            dummylist=[];
            if aap.acq_details.numdummies
                dummypath=fullfile(sesspth,'dummy_scans');
                aap=aas_makedir(aap,dummypath);
                for d=1:aap.acq_details.numdummies
                    cmd=['mv ' finalepis{d} ' ' dummypath];
                    [pth nme ext]=fileparts(finalepis{d});
                    dummylist=strvcat(dummylist,fullfile('dummy_scans',[nme ext]));
                    [s w]=aas_shell(cmd);
                    if (s)
                        aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',convertedfns{d},dummypath));
                    end
                end
            else
                d = 0;
            end
            finalepis = {finalepis{d+1:end}};
            % 4D conversion [TA]
            if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D
                finalepis = finalepis{1};
                ind = find(finalepis=='-');
                if isempty(ind), ind = find(finalepis=='.'); end; ind(2) = ind(end);
                finalepis = [finalepis(1:ind(2)-1) '.nii'];
                spm_file_merge(V0(aap.acq_details.numdummies+1:end),finalepis,0);
            end
            % And describe outputs
            aap=aas_desc_outputs(aap,subj,sess,'epi',finalepis);
            aap = aas_desc_outputs(aap,subj,sess,'dummyscans',dummylist);
            dcmhdrfn = fullfile(sesspth,'dicom_headers.mat');
            save(dcmhdrfn,'DICOMHEADERS');
            aap = aas_desc_outputs(aap,subj,sess,'epi_dicom_header',dcmhdrfn);
        end;
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
