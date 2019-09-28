% Get diffusion dicom and convert nii.gz
% Converting Data to Analyze format and extracting the Gradient Directions
% extract the gradient direction and b-values as a text file

function [aap resp]=aamod_convert_diffusion(aap,task,subjind,diffsessind,phaseencodeind)
resp='';

switch task
    case 'report'
    case 'doit'
        
        if exist('phaseencodeind','var')
            indices = [subjind diffsessind phaseencodeind];
            domain = 'diffusion_session_phaseencode_direction';
        else
            indices = [subjind diffsessind];
            domain = 'diffusion_session';
        end;
        domainpath=aas_getpath_bydomain(aap,domain,indices);
        
        % Get DICOM filenames from stream and bvals and bvecs from the
        % header
        [aap niifiles DICOMHEADERS subdirs]=aas_convertseries_fromstream(aap,domain,indices,'dicom_diffusion');
        
        if (numel(niifiles) > 1)
            for fileind=1:numel(niifiles)
                V(fileind)=spm_vol(niifiles{fileind});
            end
        else
            temp=spm_vol(niifiles{1});
            for fileind=1:numel(temp)
                V(fileind)=temp(fileind);
            end
        end            
        
        % Now move dummy scans to dummy_scans directory
        
        % From header of this module
        if isfield(aap.tasklist.currenttask.settings,'numdummies') &&...
                ~isempty(aap.tasklist.currenttask.settings.numdummies)
            ndummies=aap.tasklist.currenttask.settings.numdummies;
        else % backward compatibility
            ndummies = aap.acq_details.numdummies;
        end
        
        dummylist=[];
        if ndummies
            dummypath=fullfile(domainpath,'dummy_scans');
            aap=aas_makedir(aap,dummypath);
            for d=1:ndummies
                cmd=['mv ' niifiles{d} ' ' dummypath];
                [pth nme ext]=fileparts(niifiles{d});
                dummylist=strvcat(dummylist,fullfile('dummy_scans',[nme ext]));
                [s w]=aas_shell(cmd);
                if (s)
                    aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',niifiles{d},dummypath));
                end
            end
        else
            d = 0;
        end
        niifiles = niifiles(d+1:end);
        DICOMHEADERS = DICOMHEADERS(d+1:end);
        
        % 4D conversion [TA]

        if numel(niifiles) > 1 && isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D
            niifiles = niifiles{1};
            ind = find(niifiles=='-');
            if numel(ind) > 1, ind = ind(2); end
            niifiles = [niifiles(1:ind(1)-1) '.nii'];
            spm_file_merge(char({V(ndummies+1:end).fname}),niifiles,0,DICOMHEADERS{1}.volumeTR);
        end
        
        % bvals and bvecs based on Guy William's algorithm as implemented by Matthew Brett
        % Init
        y_flipper = diag([1 -1 1]);
        
        % Get voxel to dicom rotation matrix
        if isfield(DICOMHEADERS{1},'ImageOrientationPatient')
            orient           = reshape(DICOMHEADERS{1}.ImageOrientationPatient,[3 2]);
        elseif isfield(DICOMHEADERS{1},'PerFrameFunctionalGroupsSequence') % Philips (stack)
            orient           = reshape(DICOMHEADERS{1}.PerFrameFunctionalGroupsSequence{1}.PlaneOrientationSequence{1}.ImageOrientationPatient,[3 2]);
        end
        orient(:,3)      = null(orient');
        if det(orient)<0, orient(:,3) = -orient(:,3); end;
        vox_to_dicom = orient;
        vox_to_dicom = vox_to_dicom * y_flipper;
        
        n_hdrs = numel(DICOMHEADERS);
        bvals = zeros(n_hdrs, 1);
        bvecs = zeros(n_hdrs, 3);

        if isfield(DICOMHEADERS{1},'CSAImageHeaderInfo') % Siemens
            for h = 1:n_hdrs
                % Read B_matrix info
                bm = aas_get_numaris4_numval(DICOMHEADERS{h}.CSAImageHeaderInfo,'B_matrix')';
                % If no B_matrix, this is 0 B value
                if isempty(bm)
                    continue
                end
                B_matrix = [bm(1:3); bm(2) bm(4) bm(5); bm(3) bm(5) bm(6)];
                % find max eigenvalue, eigenvector from B_matrix
                [vecs, vals] = eig(B_matrix);
                vals = max(vals);
                [bvals(h), i] = max(vals);
                dbvec = vecs(:,i);
                % For convenience, turn vectors to point towards positive X
                if dbvec(1) < 0
                    dbvec = dbvec * -1;
                end
                bvecs(h,:) = [vox_to_dicom\dbvec]';
            end
        elseif all(isfield(DICOMHEADERS{1},{'Private_0019_10bb' 'Private_0019_10bc' 'Private_0019_10bd'})) % GE
            isb = cellfun(@(x) isfield(x,'DiffusionBValue'), DICOMHEADERS);
            if any(isb)
                bvals(1:numel(DICOMHEADERS)) = DICOMHEADERS{find(isb,1,'first')}.DiffusionBValue;
                bvals(~isb) = 0;
            elseif isfield(DICOMHEADERS{1},'Private_0043_1039') % Signa Excite 12.0 or later
                bvals = cellfun(@(x) x.Private_0043_1039(1), DICOMHEADERS);
            else
                aas_log(aap,true,'No field for b-value found!')
            end
            
            bvecs = cell2mat(cellfun(@(x) [x.Private_0019_10bb; x.Private_0019_10bc; x.Private_0019_10bd], DICOMHEADERS,'UniformOutput', false))';
            
            if strcmp(deblank(DICOMHEADERS{1}.InPlanePhaseEncodingDirection),'COL')
                % do nothing
            elseif strcmp(deblank(DICOMHEADERS{1}.InPlanePhaseEncodingDirection),'ROW')
                bvecs(:,[1 2]) = bvecs(:,[2 1]); % swap row and col
            end
            bvecs = cell2mat(arrayfun(@(x) vox_to_dicom\bvecs(x,:)', 1:size(bvecs,1),'UniformOutput',0))'; % rotate to scanner space
        elseif isfield(DICOMHEADERS{1},'PerFrameFunctionalGroupsSequence') % Philips (stack)
            H = cell2mat(DICOMHEADERS{1}.PerFrameFunctionalGroupsSequence);
            pos = cell2mat([H.FrameContentSequence]);
            pos = reshape([pos.DimensionIndexValues],4,[])';
            H = H(pos(:,2) == 1); % first slices
            
            bvals = arrayfun(@(x) x.MRDiffusionSequence{1}.DiffusionBValue, H);
            bvecs(bvals~=0,:) = cell2mat(arrayfun(@(x) x.MRDiffusionSequence{1}.DiffusionGradientDirectionSequence{1}.DiffusionGradientOrientation', H(bvals~=0),'UniformOutput',false))';
        end        
        
        % Output final data
        sesspth=aas_getpath_bydomain(aap,domain,indices);
        
        % Write bvals
        bvals_fn=fullfile(sesspth,'bvals');
        fid=fopen(bvals_fn,'w');
        fprintf(fid,'%d ',bvals);
        fprintf(fid,'\n');
        fclose(fid);
        
        % Write bvecs
        bvecs_fn=fullfile(sesspth,'bvecs');
        fid=fopen(bvecs_fn,'w');
        for ln=1:3
            fprintf(fid,'%.14f ',bvecs(:,ln));
            fprintf(fid,'\n');
        end;
        fclose(fid);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,domain,indices,'dummyscans',dummylist);
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_data',niifiles);
        dcmhdrfn = fullfile(domainpath,'dicom_headers.mat');
        save(dcmhdrfn,'DICOMHEADERS');
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_dicom_header',dcmhdrfn);
        aap=aas_desc_outputs(aap,domain,indices,'bvals',bvals_fn);
        aap=aas_desc_outputs(aap,domain,indices,'bvecs',bvecs_fn);    
       
end
end

