% function [V Y XYZ nimg]=aas_spm_vol(allfn,dontdelete)
% Use like spm_vol, except that:
%  (1) it will read a single 3D, list of 3D or a single 4D file, returning an array of headers in V, and a 4d
%   array in Y with final dimension "time". nimg contains number of time
%   points
%  (2) it will read .nii.gz files
%
% Inputs:
%  allfn - cell array or char array of filenames
%  dontdelete - boolean to specify don't delete compressed file (useful if
%   you will reuse it manually), default false (or delete the file)
% Outputs:
%  V - header as returned by spm_vol (array of structures if multiple
%   timepoints)
%  Y - data, 3d or 4d
%  XYZ - coordinates, if you provide a third output argument
%  nimg - number of timepoints
%
% Useful for FSL integration.
% Rhodri Cusack, BMI Western, Canada, June 2013-16

function [V Y XYZ nimg]=aas_spm_vol(allfn,dontdelete)
if ~exist('dontdelete','var')
    dontdelete=false;
end;

% Cell or char array of filenames, get length
if iscell(allfn)
    ninputfn=length(allfn);
else
    ninputfn=size(allfn,1);
end;

% May be multiple input files
for imgind=1:ninputfn
    % Cell or char array of filenames
    if iscell(allfn)
        fn=allfn{imgind};
    else
        fn=allfn(imgind,:);
    end;
    [pth nme ext]=aas_fileparts(fn);
    
    % Uncompress .gz if needed
    if strcmp(ext,'.nii.gz')
        wascompressed=true;
        % Uncompress and read
        [temp_pth temp_nme temp_ext]=fileparts(tempname);
        tf=fullfile(pth,[temp_nme ext]);
        copyfile(deblank(fn),tf);
        [s w]=aas_shell(['gunzip ' tf]);
        niifn=tf(1:end-3);
    else
        wascompressed=false;
        niifn=fn;
    end;
    
    if ninputfn>1
        % List of 3D files
        V(imgind)=spm_vol(niifn);
        if nargout>=3
            [Y(:,:,:,imgind) XYZ]=spm_read_vols(V(imgind));
        else
            Y(:,:,:,imgind)=spm_read_vols(V(imgind));
        end;
    else
        % Single 3D or 4D file
        V=spm_vol(niifn);
        if nargout>=3
            [Y XYZ]=spm_read_vols(V);
        else
            Y=spm_read_vols(V);
        end;
    end;
    
    % Tidy up if needed
    if wascompressed && ~dontdelete
        delete(tf(1:end-3));
    end;
    
end

% Get number of timepoints
V=V(:);
nimg=size(Y,4);
end
