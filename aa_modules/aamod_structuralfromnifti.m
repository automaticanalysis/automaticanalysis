% AA module - structural from NIFTI

function [aap,resp]=aamod_structuralfromnifti(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        allseries = horzcat(aap.acq_details.subjects(subj).structural{:}); % uniform input assumed
        if isstruct(allseries{1})
            allseries = cell2mat(allseries);
        elseif iscellstr(allseries)
            allseries = cellfun(@(x) struct('fname',x,'hdr',''), allseries);
        else
            aas_log(aap,true,'ERROR: aap.acq_details.subjects.structural has a wrong format')
            help aas_addsubject
        end
        allstreams = aas_getstreams(aap,'output'); allstreams(cell_index(allstreams,'dicom_header')) = [];
        sfxs = textscan(aas_getsetting(aap,'sfxformodality'),'%s','delimiter',':'); sfxs = sfxs{1}';
        
        for m = 1:numel(sfxs)
            stream = allstreams{m};
            
            % Select
            series = allseries(cellfun(@(x) ~isempty(x), strfind({allseries.fname},sfxs{m})));
            switch numel(series)
                case 0
                    aas_log(aap,false,['WARNING: no ' stream ' image found']);
                    continue;
                case 1
                otherwise
                    if aap.options.(['autoidentify' stream '_choosefirst']), series = series(1); end
                    if aap.options.(['autoidentify' stream '_chooselast']), series = series(end); end
                    if (numel(series) > 1) && ...
                            ~aap.options.(['autoidentify' stream '_multiple']) && ~aap.options.(['autoidentify' stream '_average'])
                        aas_log(aap,true,sprintf('ERROR: %d %s image found but one expected',numel(series),stream));
                    end
            end
            
            sesspth=fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.structdirname);
            aas_makedir(aap,sesspth);
            Ys = 0;
            for s = 1:numel(series)
                niftifile = series(s).fname;
                hdrfile = series(s).hdr;
                
                %% Image
                if ~exist(niftifile,'file')
                    niftisearchpth=aas_findvol(aap,subj);
                    if ~isempty(niftisearchpth)
                        niftifile = fullfile(niftisearchpth,niftifile);
                        if ~exist(niftifile,'file'), aas_log(aap,true,['ERROR: image ' niftifile ' not found']); end
                    end
                end                
                comp = false;
                if strcmp(spm_file(niftifile,'Ext'),'gz')
                    comp = true;
                    gunzip(niftifile,fullfile(sesspth,'temp'));
                    niftifile = spm_file(niftifile(1:end-3),'path',fullfile(sesspth,'temp'));
                end
                
                fn(s,:)=spm_file(niftifile,'path',sesspth,'suffix','_0001');
                V(s)=spm_vol(niftifile);
                Y=spm_read_vols(V(s));  
                V(s).fname=deblank(fn(s,:));
                V(s).n=[1 1];
                if comp, rmdir(fullfile(sesspth,'temp'),'s'); end                
                
                if aap.options.(['autoidentify' stream '_average'])
                    Ys = Ys + Y/numel(series);
                else                    
                    spm_write_vol(V(s),Y);
                end
                
                %% header
                if ~isempty(hdrfile)
                    if isstruct(hdrfile)
                        dcmhdr = hdrfile;
                    else
                        switch spm_file(hdrfile,'Ext')
                            case 'mat'
                                dat = load(hdrfile);
                                dcmhdr = dat.dcmhdr;
                            case 'json'
                                hdrfile = loadjson(hdrfile);
                                % convert timings to ms (DICOM default)
                                for f = fieldnames(hdrfile)'
                                    if strfind(f{1},'Time'), dcmhdr{s}.(f{1}) = hdrfile.(f{1})*1000; end
                                end
                                if isfield(dcmhdr{s},'RepetitionTime'), dcmhdr{s}.volumeTR = dcmhdr{s}.RepetitionTime/1000; end
                                if isfield(dcmhdr{s},'EchoTime'), dcmhdr{s}.volumeTE = dcmhdr{s}.EchoTime/1000; end                    
                        end
                    end
                else
                    aas_log(aap,false,'WARNING: No header provided!');
                end
            end
            if aap.options.(['autoidentify' stream '_average'])
                fn = deblank(fn(1,:));
                V = V(1); 
                V.fname = fn;
                spm_write_vol(V,Ys);
                dcmhdr = dcmhdr(1);
            end
            aap=aas_desc_outputs(aap,subj,stream,fn);
            if exist('dcmhdr','var')
                dcmhdrfn=fullfile(sesspth,[stream '_dicom_header.mat']);
                save(dcmhdrfn,'dcmhdr');
                aap=aas_desc_outputs(aap,subj,[stream '_dicom_header'],dcmhdrfn);
            end
            % remove variable to avoid writing wrong header to subsequent sfxs
            clear dcmhdr
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
