% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_specialseries(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        % From header of this module
        if ~isempty(aas_getsetting(aap,'numdummies',sess))
            ndummies=aas_getsetting(aap,'numdummies',sess);
        else % backward compatibility
            ndummies = aap.acq_details.numdummies;
        end
        if ~isempty(aas_getsetting(aap,'NIFTI4D',sess))
            NIFTI4D=aas_getsetting(aap,'NIFTI4D',sess);
        else % backward compatibility
            NIFTI4D = aap.options.NIFTI4D;
        end
        
        
        for inpstream = aas_getstreams(aap,'input')
            if ~aas_stream_has_contents(aap,'special_session',[subj sess],inpstream{1}), continue; end
            
            streamname = strrep(inpstream{1},'dicom_','');
            
            [aap, convertedfns, dcmhdr] = aas_convertseries_fromstream(aap, subj, sess, ['dicom_' streamname]);
            
            outstream = {};
            % Restructure outputs!
            if iscell(convertedfns{1}) % subdirs
                for c = 1:length(convertedfns)
                    outstream = [outstream; convertedfns{c}];
                    dcmhdr{c} = header_unique(dcmhdr{c}); 
                end
            else
                outstream = convertedfns;
            end
            
            % select images in the same folder
            dirs = unique(spm_file(outstream,'path'))';
            for d = 1:numel(dirs)
                ind4D = cellfun(@(x) ~isempty(x), strfind(outstream,dirs{d}));
                if sum(ind4D) > 1
                    img4D = outstream(ind4D);
                    
                    % dummies
                    if ndummies
                        dummypath=fullfile(dirs{d},'dummy_scans');
                        aap=aas_makedir(aap,dummypath);
                        for dummyind=1:ndummies
                            cmd=['mv ' img4D{dummyind} ' ' dummypath];
                            [s w]=aas_shell(cmd);
                            if (s)
                                aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',img4D{dummyind},dummypath));
                            end
                        end
                        img4D(1:ndummies) = [];
                    end
                    
                    if NIFTI4D
                        fname4D = [getuniformfields(img4D) '.nii'];
                        spm_file_merge(char(img4D),fname4D,0,dcmhdr{d}.volumeTR);
                        img4D = cellstr(fname4D);
                        % replace first instance with the 4D and remove the rest
                        outstream(find(ind4D,1,'first')) = img4D;
                        ind4D(find(ind4D,1,'first')) = false;
                        outstream(ind4D) = [];
                    end
                end
            end
            
            aap = aas_desc_outputs(aap, 'special_session', [subj, sess], streamname, outstream);
            dcmhdrfn = fullfile(aas_getsesspath(aap,subj,sess),[streamname '_dicom_headers.mat']);
            save(dcmhdrfn,'dcmhdr');
            aap = aas_desc_outputs(aap, 'special_session', [subj, sess], [streamname '_dicom_header'], dcmhdrfn);
        end
        if ~exist('streamname','var'), aas_log(aap,true,sprintf('No input found for special session %s',aas_getsessname(aap,sess))); end
    case 'checkrequirements'
%          % get input
%         [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
%         stageindex = sscanf(index,'_%05d');
%         in = aap.tasksettings.(stagename)(stageindex).inputstreams.stream; if ~iscell(in), in = {in}; end
% 
%         if ~any(strcmp(in,['dicom_' aas_getsessname(aap,sess)]))
%             aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',['dicom_' aas_getsessname(aap,sess)],'input');
%             aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' input stream: ''' ['dicom_' aas_getsessname(aap,sess)] '''']);
%         end
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function hdrout = header_unique(hdr)
% based on spm_dicom_convert.m 6639 2015-12-11 09:50:49Z volkmar

% Single DICOM
if numel(hdr) == 1
    hdrout = hdr{1}; return;
end

% Initialise structure
hdrout = hdr{1}; hdrout(1) = [];
i = 1;

while i < numel(hdr)
    hdrout(end+1) = hdr{i};
    new = false;
    for i = i:numel(hdr)
        if isfield(hdr{i},'ImageType') && isfield(hdrout(end), 'ImageType')
            new = new || ~strcmp(hdr{i}.ImageType, hdrout(end).ImageType);
        end
        if isfield(hdr{i},'SequenceName') && isfield(hdrout(end), 'SequenceName')
            new = new || ~strcmp(hdr{i}.SequenceName,hdrout(end).SequenceName);
        end
        if isfield(hdr{i},'SeriesInstanceUID') && isfield(hdrout(end), 'SeriesInstanceUID')
            new = new || ~strcmp(hdr{i}.SeriesInstanceUID,hdrout(end).SeriesInstanceUID);
        end
        if isfield(hdr{i},'EchoNumbers')  && isfield(hdrout(end), 'EchoNumbers')
            new = new || ~(hdr{i}.EchoNumbers == hdrout(end).EchoNumbers);
        end
        if isfield(hdr{i},'GE_ImageType')  && isfield(hdrout(end), 'GE_ImageType')
            new = new || ~(hdr{i}.GE_ImageType == hdrout(end).GE_ImageType);
        end
        if new, break; end
    end    
end
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