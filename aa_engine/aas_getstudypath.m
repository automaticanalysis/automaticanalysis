% aa - get study path
% examples:
%  pth=aas_getstudypath(aap);
%  pth=aas_getstudypath(aap,'s3');   % path on s3
%  pth=aas_getstudypath(aap,2);      % path for module 2
%  pth=aas_getstudypath(aap,'s3',2);


function [root]=aas_getstudypath(aap,varargin)
global aaworker


v=varargin;

% First, see if filesystem explicitly specified
if (~isempty(v) && ischar(v{1}))
    remotefilesystem=v{1};
    v=v(2:end);
else
    remotefilesystem='none';
end;

% Numeric argument on the end corresponds to a specified input source
if (~isempty(v))
    k=v{end};
    
    % Get analysis id and analysis id suffix
    if isfield(aap, 'internal')
        analysisid=aap.internal.aap_initial.directory_conventions.analysisid;
    else
        analysisid=aap.directory_conventions.analysisid;
    end
    try
        analysisid_suffix=aap.tasklist.main.module(k).extraparameters.aap.directory_conventions.analysisid_suffix;
    catch
        analysisid_suffix=aap.internal.aap_initial.directory_conventions.analysisid_suffix;
    end;
    
    % Get the basic root directory for the current filesystem
    switch (remotefilesystem)
        case 'none'
            try
                if isfield(aap, 'internal')
                    root=aap.internal.aap_initial.acq_details.root;
                else
                    root = aap.acq_details.root;
                end
            catch
                root=fullfile('/cn_data',aaworker.username,aap.acq_details.root);
            end;
        otherwise
            try
                root=aap.internal.aap_initial.acq_details.(remotefilesystem).root;
            catch
                root=aap.acq_details.(remotefilesystem).root;
            end;
    end;
    
    % Add suffixes
    root=fullfile(root,[analysisid analysisid_suffix]);
    % Add split by module layer
    
    outputformat=aap.directory_conventions.outputformat;
    
    switch(outputformat)
        case 'splitbymodule'
            suffix=aas_getstagetag(aap,k);
        otherwise
            suffix='';
    end;
    root=fullfile(root,suffix);
else
    % otherwise, just use the root we've been given - will already have
    % been adorned approriately for the current module
    switch (remotefilesystem)
        case 'none'
            root=aap.acq_details.root;
        otherwise
            root=aap.acq_details.(remotefilesystem).root;
    end;
end;

