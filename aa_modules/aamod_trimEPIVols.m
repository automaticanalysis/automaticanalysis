function [aap,resp] = aamod_trimEPIVols(aap, task, subj, sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        sessDir = aas_getsesspath(aap, subj, sess);
        
        fileNames = aas_getimages_bystream(aap, subj, sess, 'epi');
        numImgs = size(fileNames, 1);
        volumeI = zeros(numImgs, 1);      % indices of volumes to drop
        
        V = spm_vol(fileNames(1, :));
        sliceI = zeros(V.dim(3), 1);      % indices of slices to drop
        
        trimSettings = aap.tasklist.currenttask.settings.data;
        
        for d = 1 : length(trimSettings)
            
            if ischar(trimSettings(d).subject)
                if strcmp(trimSettings(d).subject, '*');
                    subjI = subj;
                else
                    subjI = find(strcmp({aap.acq_details.subjects}, trimSettings(d).subject));
                end
            else
                subjI = trimSettings(d).subject;
            end
            
            if ischar(trimSettings(d).session)
                if strcmp(trimSettings(d).session, '*');
                    sessI = sess;
                end
            else
                sessI = trimSettings(d).session;
            end
            
            if (subj == subjI) && (sessI == sess)
               
                dropImgs = trimSettings(d).dropImgs;
                if ischar(dropImgs)
                    dropImgs = str2num(dropImgs);
                end
           
                volumeI(dropImgs) = 1;
                
                trimSlices = trimSettings(d).trimSlices;
                if ischar(trimSlices)
                    trimSlices = str2num(trimSlices);
                end
                
                sliceI(trimSlices) = 1;          
            end
        end
        
        outImgs = [];
        
        fprintf('%60s', '');
        
        for f = 1 : numImgs
            
            fprintf('%s%-60s', repmat(sprintf('\b'),1,60), sprintf('Trimming Image %d/%d', f, size(fileNames,1)));
            
            if ~(volumeI(f))
                
                V = spm_vol(fileNames(f, :));
                [path file ext] = fileparts(V.fname);
                keepSlicesInd = find(~sliceI);
                Y = spm_read_vols(V);
                Y = Y(:, :, keepSlicesInd);
                V.dim(3) = V.dim(3)-length(find(sliceI));
                V.fname = fullfile(sessDir, ['t' file ext]);
                spm_write_vol(V, Y);
                outImgs = strvcat(outImgs, V.fname);
            end
            
        end
        
        fprintf('\nDropped %d images.\n', sum(volumeI));
        
        aap = aas_desc_outputs(aap, subj, sess,'epi', outImgs);
             
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end