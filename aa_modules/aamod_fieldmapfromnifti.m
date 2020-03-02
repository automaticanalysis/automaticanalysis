% AA module - fieldmap from NIFTI

function [aap,resp]=aamod_fieldmapfromnifti(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        sesspth=fullfile(aas_getsesspath(aap,subj,sess),aap.directory_conventions.fieldmapsdirname);
        aas_makedir(aap,sesspth);
               
        %% locate
        if ~iscell(aap.acq_details.subjects(subj).fieldmaps)
            aas_log(aap,true,'Was exepcting list of filenames in cell array');
        end;
        niftistruct = [];
        % try session-specific --> visit-specific --> rest
        d = aas_get_series(aap,strtok(aas_getsesstype(aap),'_'),subj,sess);
        fieldmaps = [aap.acq_details.subjects(subj).fieldmaps(d:end) aap.acq_details.subjects(subj).fieldmaps(1:d-1)];
        for n = horzcat(fieldmaps{:})
            if any(strcmp(n{1}.session,aap.acq_details.([aas_getsesstype(aap) 's'])(sess).name)) || strcmp(n{1}.session,'*'), niftistruct = n{1}; end
        end
        if isempty(niftistruct), aas_log(aap,true,sprintf('ERROR: No fieldmap found for session %s',aap.acq_details.sessions(sess).name)); end

        niftifile = niftistruct.fname; % checks first only
        hdrfile = niftistruct.hdr;
        if ~exist(niftifile{1},'file')
            niftisearchpth=aas_findvol(aap,'');
            if ~isempty(niftisearchpth)
                niftifile = spm_file(niftifile,'path',niftisearchpth);
                if ischar(hdrfile), hdrfile = spm_file(hdrfile,'path',niftisearchpth); end
            end
        end
        
        %% images
        for f = 1:numel(niftifile)
            comp = strcmp(spm_file(niftifile{f},'Ext'),'gz');
            if comp
                gunzip(niftifile{f},fullfile(sesspth,'temp'));
                niftifile{f} = spm_file(niftifile{f}(1:end-3),'path',fullfile(sesspth,'temp'));
            end
            
            aas_makedir(aap,fullfile(sesspth,sprintf('serie%02d',f)));
            V=spm_vol(niftifile{f});
            for iV = V'
                Y=spm_read_vols(iV);
                iV.fname = spm_file(niftifile{f},'path',fullfile(sesspth,sprintf('serie%02d',f)));
                spm_write_vol(iV,Y);
            end
            fn{f} = iV.fname;
            
            if comp, rmdir(fullfile(sesspth,'temp'),'s'); end
        end
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],'fieldmap',char(fn));
        
        %% header
        if ~isempty(hdrfile)
            if iscell(hdrfile)
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
                            if strfind(f{1},'Time'), dcmhdr{1}.(f{1}) = hdrfile.(f{1})*1000; end
                        end
                        if isfield(dcmhdr{1},'RepetitionTime'), dcmhdr{1}.volumeTR = dcmhdr{1}.RepetitionTime/1000; end
                        if isfield(dcmhdr{1},'EchoTime1') && isfield(dcmhdr{1},'EchoTime2'), dcmhdr{1}.volumeTE = [dcmhdr{1}.EchoTime1 dcmhdr{1}.EchoTime2]/1000; end
                end
            end
        end
        
        %% Output
        dcmhdrfn=fullfile(aas_getsesspath(aap,subj,sess),'fieldmap_dicom_header.mat');
        save(dcmhdrfn,'dcmhdr');
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],'fieldmap_dicom_header',dcmhdrfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end