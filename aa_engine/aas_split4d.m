%aas_split4d(inFile) reads in a 4d file and writes out separate 3d
%files. Useful when, e.g., working with EPI data.
%
% outFiles = aas_split4d returns full paths of the output files.
%
% aas_split4d(inFile, true) deletes the 4d file when done. The default is to
% not delte anything.
%
% If the requested file is actually a 3d file, it is returned unchanged and
% nothing is deleted, regardless of what was requested.
function outFiles = aas_split4d(inFile, deleteafter)

if nargin < 2 || isempty(deleteafter)
    deleteafter = true;
end

useFSL = true;

Vin = spm_vol(inFile);
[pth, nm, ext] = fileparts(inFile);


% If only one volume, don't split, but return this same file. Otherwise,
% split, and delete 4d (if requested).
if length(Vin)==1
    outFiles = inFile;
else
    
    outFiles = [];
    
    if useFSL        
        cmd = sprintf('fslsplit %s %s_ -t', inFile, fullfile(pth, nm));
        aas_runfslcommand(aap, cmd);
        outFiles = spm_select('fplist', pth, sprintf('%s_[0-9]*.nii', nm));
    else
        % choose how many places to write out volumes (important if more than 4)
        nZeros = max(4, length(num2str(length(Vin))));
        formatString = sprintf('%%s_%%0%dd%%s', nZeros);
        
        for volInd = 1:length(Vin)
            Y = spm_read_vols(Vin(volInd));
            thisFile = fullfile(pth, sprintf(formatString, nm, volInd, ext));
            Vin(volInd).fname = thisFile;
            Vin(volInd).n = [1 1];
            spm_write_vol(Vin(volInd), Y);
            
            outFiles = strvcat(outFiles, thisFile);
        end
        
        
    end
    
    if deleteafter
        delete(inFile);
    end
end