% aa - get study path
% examples:
%  pth=aas_getstudypath(aap);
%  pth=aas_getstudypath(aap,'s3');   % path on s3
%  pth=aas_getstudypath(aap,2);      % path for module 2
%  pth=aas_getstudypath(aap,'s3',2);


function [root]=aas_getstudypath(aap,varargin)

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
    
    % Get analysis id and analysis id suffix, ensuring any
    % stage-specific overloading is applied
    try
        extra_aap=aap.tasklist.main.module(k).extraparameters.aap;
    catch
    end;
    try
        analysisid=extra_aap.directory_conventions.analysisid;
    catch
        analysisid=aap.internal.aap_initial.directory_conventions.analysisid;
    end;
    try
        analysisid_suffix=extra_aap.directory_conventions.analysisid_suffix;
    catch
        analysisid_suffix=aap.internal.aap_initial.directory_conventions.analysisid_suffix;
    end;
    
    % Check on the filesystem
    switch (remotefilesystem)
        case 'none'
            try
                root=extra_aap.acq_details.root;
            catch
                root=aap.internal.aap_initial.acq_details.root;
            end;
        otherwise
            
            root=aap.acq_details.(remotefilesystem).root;
    end;
    
    % Add suffixes
    root=fullfile(root,[analysisid analysisid_suffix]);
    
    % Add split by module layer
    switch(aap.directory_conventions.outputformat)
        case 'splitbymodule'
            suffix=aas_getstagetag(aap,k);
        otherwise
            suffix='';
    end;
else
    % otherwise, just use the root we've been given
    root=aap.acq_details.root;
    suffix='';
    
end;



% % Is there a character argument - taken to be a type of filesystem
% %  (e.g., s3) that will have its own root
% if (~isempty(v) && ~isempty(v{1}))
%     if (ischar(varargin{1}))
%         if (~strcmp(v{1},'none'))
%             root=aap.acq_details.(v{1}).root;
%         end;
%         v=v(2:end);
%     end;
% end;




root=fullfile(root,suffix);

