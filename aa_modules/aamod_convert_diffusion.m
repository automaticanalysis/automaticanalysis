% Get diffusion dicom and convert nii.gz
% Converting Data to Analyze format and extracting the Gradient Directions
% extract the gradient direction and b-values as a text file

function [aap resp]=aamod_convert_diffusion(aap,task,subjind,diffsessind)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        
        % Get DICOM filenames from stream and bvals and bvecs from the
        % header
        [aap niifiles dicomheader subdirs]=aas_convertseries_fromstream(aap,'diffusion_session',[subjind diffsessind],'dicom_diffusion');
        bvecs=zeros(length(dicomheader),3);
        for headerind=1:length(dicomheader)
                bvals(headerind)=aas_get_numaris4_numval(dicomheader{headerind}.CSAImageHeaderInfo,'B_value');
                if bvals(headerind)~=0
                    bvecs(headerind,:)=aas_get_numaris4_numval(dicomheader{headerind}.CSAImageHeaderInfo,'DiffusionGradientDirection');
                end;
        end;
        
        % Output final data
        sesspth=aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]);
        
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
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'diffusion_data',niifiles);
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'bvals',bvals_fn);
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'bvecs',bvecs_fn);    
       
end
end

