% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_specialseries(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        % From header of this module
        if ~isempty(aas_getsetting(aap,'numdummies'))
            ndummies=aas_getsetting(aap,'numdummies');
        else % backward compatibility
            ndummies = aap.acq_details.numdummies;
        end
        
        for inpstream = aas_getstreams(aap,'input')            
            streamname = strrep(inpstream{1},'dicom_','');
            
            [aap, convertedfns, dcmhdr] = aas_convertseries_fromstream(aap, subj, sess, ['dicom_' streamname]);
            
            outstream = {};
            % Restructure outputs!
            if iscell(convertedfns{1}) % subdirs
                for c = 1:length(convertedfns)
                    outstream = [outstream; convertedfns{c}];
                    dcmhdr{c} = dcmhdr{c}{1}; 
                end
            else
                outstream = convertedfns;
            end
            
            % select images in the same folder
            dirs = unique(spm_file(outstream,'path'))';
            for d = dirs
                ind4D = cellfun(@(x) ~isempty(x), strfind(outstream,d{1}));
                if sum(ind4D) > 1
                    img4D = outstream(ind4D);
                    outstream(ind4D) = [];
                    
                    % dummies
                    if ndummies
                        dummypath=fullfile(d{1},'dummy_scans');
                        aap=aas_makedir(aap,dummypath);
                        for dummyind=1:ndummies
                            cmd=['mv ' img4D{dummyind} ' ' dummypath];
                            [s w]=aas_shell(cmd);
                            if (s)
                                aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',niifiles{dummyind},dummypath));
                            end
                        end
                        img4D(1:ndummies) = [];
                    end
                    
                    if aap.options.NIFTI4D
                        fname4D = getuniformfields(img4D);
                        spm_file_merge(char(img4D),fname4D,0);
                        img4D = cellstr(fname4D);
                    end
                    outstream = [outstream; img4D];
                end
            end
            
            aap = aas_desc_outputs(aap, 'special_session', [subj, sess], streamname, outstream);
            dcmhdrfn = fullfile(aas_getsesspath(aap,subj,sess),[streamname '_dicom_headers.mat']);
            save(dcmhdrfn,'dcmhdr');
            aap = aas_desc_outputs(aap, 'special_session', [subj, sess], [streamname '_dicom_header'], dcmhdrfn);
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function outfname = getuniformfields(fnames,delimiter)
if nargin < 2, delimiter = '-'; end
fields = cell(0);
for f = 1:numel(fnames)
    l = textscan(fnames{f},'%s','Delimiter',delimiter);
    fields = vertcat(fields,l{1}');
end
outfname = '';
for f = 1:size(fields,2)
    if numel(unique(fields(:,f))) == 1, outfname = [outfname char(unique(fields(:,f))) delimiter]; end
end
outfname(end) = [];
end
