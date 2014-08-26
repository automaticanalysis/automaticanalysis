% Get diffusion dicom and convert nii.gz
% Converting Data to Analyze format and extracting the Gradient Directions
% extract the gradient direction and b-values as a text file

function [aap resp]=aamod_convert_diffusion(aap,task,subjind,diffsessind)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        
        domain = 'diffusion_session';
        indices = [subjind diffsessind];
        domainpath=aas_getpath_bydomain(aap,domain,indices);
        
        % Get DICOM filenames from stream and bvals and bvecs from the
        % header
        [aap niifiles dicomheader subdirs]=aas_convertseries_fromstream(aap,domain,indices,'dicom_diffusion');
        
        % Now move dummy scans to dummy_scans directory
        
        dummylist=[];
        if aap.acq_details.numdummies
            dummypath=fullfile(domainpath,'dummy_scans');
            aap=aas_makedir(aap,dummypath);
            for ndummies=1:aap.acq_details.numdummies
                cmd=['mv ' niifiles{ndummies} ' ' dummypath];
                [pth nme ext]=fileparts(niifiles{ndummies});
                dummylist=strvcat(dummylist,fullfile('dummy_scans',[nme ext]));
                [s w]=aas_shell(cmd);
                if (s)
                    aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',niifiles{ndummies},dummypath));
                end
            end
        else
            ndummies = 0;
        end
        niifiles = {niifiles{ndummies+1:end}};
        % 4D conversion [TA]
        for fileind=1:numel(niifiles)
            V(fileind)=spm_vol(niifiles{fileind});
        end
        if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D
            niifiles = niifiles{1};
            ind = find(niifiles=='-');
            niifiles = [niifiles(1:ind(2)-1) '.nii'];
			spm_file_merge(char({V.fname}),niifiles,0);
        end
        
        bvecs=zeros(length(dicomheader),3);
        for headerind=1:length(dicomheader)
                bvals(headerind)=aas_get_numaris4_numval(dicomheader{headerind}.CSAImageHeaderInfo,'B_value');
                if bvals(headerind)~=0
                    bvecs(headerind,:)=aas_get_numaris4_numval(dicomheader{headerind}.CSAImageHeaderInfo,'DiffusionGradientDirection');
                end;
        end;
        
        % Output final data
        sesspth=aas_getpath_bydomain(aap,domain,indices);
        
        % Write bvals
        bvals_fn=fullfile(sesspth,'bvals');
        fid=fopen(bvals_fn,'w');
        fprintf(fid,'%d ',bvals(ndummies+1:end));
        fprintf(fid,'\n');
        fclose(fid);
        
        % Write bvecs
        bvecs_fn=fullfile(sesspth,'bvecs');
        fid=fopen(bvecs_fn,'w');
        for ln=1:3
            fprintf(fid,'%.14f ',bvecs(ndummies+1:end,ln));
            fprintf(fid,'\n');
        end;
        fclose(fid);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,domain,indices,'dummyscans',dummylist);
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_data',niifiles);
        aap=aas_desc_outputs(aap,domain,indices,'bvals',bvals_fn);
        aap=aas_desc_outputs(aap,domain,indices,'bvecs',bvecs_fn);    
       
end
end

