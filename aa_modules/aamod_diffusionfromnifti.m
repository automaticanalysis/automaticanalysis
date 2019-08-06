% AA module - diffusion from NIFTI

function [aap,resp]=aamod_diffusionfromnifti(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        %% Init
        domain = aap.tasklist.currenttask.domain;
        indices = [subj sess];
        sesspth = aas_getpath_bydomain(aap,domain,indices);
        numdummies = aap.acq_details.numdummies;
        if ~isempty(aap.tasklist.currenttask.settings.numdummies)
            numdummies=aap.tasklist.currenttask.settings.numdummies;
        end
        
        %% Select
        series = horzcat(aap.acq_details.subjects(subj).diffusion_seriesnumbers{:});
        if ~iscell(series) || ~isstruct(series{1})
            aas_log(aap,true,'ERROR: Was expecting list of struct(s) of fname+bval+bvec in cell array');
        end
        series = series{sess};
        
        %% Process - assume 4D
        finalepis={}; 
        % Files
        niftifile = series.fname;
        headerfile = series.hdr;
        bvalfile = series.bval;
        bvecfile = series.bvec;
        if ~exist(niftifile,'file') % assume path realtive to (first) rawdatadir
            niftisearchpth=aas_findvol(aap,'');
            if ~isempty(niftisearchpth)
                niftifile = fullfile(niftisearchpth,niftifile);
                headerfile = fullfile(niftisearchpth,headerfile);
                bvalfile = fullfile(niftisearchpth,bvalfile);
                bvecfile = fullfile(niftisearchpth,bvecfile);
            end
        end
        comp = false;
        if strcmp(spm_file(niftifile,'Ext'),'gz')
            comp = true;
            gunzip(niftifile,fullfile(sesspth,'temp'));
            niftifile = spm_file(niftifile(1:end-3),'path',fullfile(sesspth,'temp'));
        end
        
        % BVals/BVecs
        bvals = dlmread(bvalfile); if size(bvals,2) > size(bvals,1), bvals = bvals'; end
        bvecs = dlmread(bvecfile); if size(bvecs,2) > size(bvecs,1), bvecs = bvecs'; end
        
        bvals_fn=fullfile(sesspth,'bvals');
        fid=fopen(bvals_fn,'w');
        fprintf(fid,'%d ',bvals);
        fprintf(fid,'\n');
        fclose(fid);
        
        bvecs_fn=fullfile(sesspth,'bvecs');
        fid=fopen(bvecs_fn,'w');
        for ln=1:3
            fprintf(fid,'%.14f ',bvecs(:,ln));
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        % Header
        DICOMHEADERS{1} = struct;
        if ~isempty(headerfile)
            if strcmp(spm_file(headerfile,'Ext'),'mat')
                dat = load(headerfile);
                DICOMHEADERS = dat.DICOMHEADERS;
            else
                switch spm_file(headerfile,'Ext')
                    case 'dcmhdr'
                        DICOMHEADERS{1} = header_dcmhdr(headerfile);
                    case 'json'
                        DICOMHEADERS{1} = header_json(headerfile); % BIDS
                end
                DICOMHEADERS{1}.volumeTR = DICOMHEADERS{1}.RepetitionTime/1000;
                DICOMHEADERS{1}.volumeTE = DICOMHEADERS{1}.EchoTime/1000;
                DICOMHEADERS{1}.slicetimes = DICOMHEADERS{1}.SliceTiming/1000;
                [junk, DICOMHEADERS{1}.sliceorder] = sort(DICOMHEADERS{1}.slicetimes);
                DICOMHEADERS{1}.echospacing = DICOMHEADERS{1}.EchoSpacing/1000;
            end
        else
            aas_log(aap,false,'WARNING: No header provided!');
        end 
        
        % Image
        V = spm_vol(niftifile);
        aas_makedir(aap,sesspth);
        [junk, fle, ext]=fileparts(niftifile);
        for fileind=1:numel(V)
            Y=spm_read_vols(V(fileind));
            fn=fullfile(sesspth,[sprintf('%s_%04d',fle,fileind) ext]);
            V(fileind).fname=fn;
            V(fileind).n=[1 1];
            spm_write_vol(V(fileind),Y);
            finalepis=[finalepis fn];
        end
        
        % Write out the files        
        % Now move dummy scans to dummy_scans directory
        dummylist=[];
        if numdummies
            dummypath=fullfile(sesspth,'dummy_scans');
            aap=aas_makedir(aap,dummypath);
            for d=1:numdummies
                cmd=['mv ' finalepis{d} ' ' dummypath];
                dummylist=strvcat(dummylist,spm_file(finalepis{d},'path','dummy_scans'));
                s=aas_shell(cmd);
                if (s)
                    aas_log(aap,1,sprintf('ERROR:Problem moving dummy scan\nERROR:    %s to\nERROR:    %s\n',finalepis{d},dummypath));
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
            sfx = '';
            if isempty(ind)
                ind = find(finalepis=='.');
                sfx = '_4D';
            end
            ind(2) = ind(end);
            finalepis = [finalepis(1:ind(2)-1) sfx '.nii'];
            if iscell(V), V = cell2mat(V); end
            spm_file_merge(char({V(numdummies+1:end).fname}),finalepis,0);
        end
        
        %% Describe outputs
        if comp, rmdir(fullfile(sesspth,'temp'),'s'); end
        aap=aas_desc_outputs(aap,domain,indices,'dummyscans',dummylist);
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_data',finalepis);
        dcmhdrfn = fullfile(sesspth,'dicom_headers.mat');
        save(dcmhdrfn,'DICOMHEADERS');
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_dicom_header',dcmhdrfn);
        aap=aas_desc_outputs(aap,domain,indices,'bvals',bvals_fn);
        aap=aas_desc_outputs(aap,domain,indices,'bvecs',bvecs_fn);  
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end