% aacloudclient launch script
%  Rhodri Cusack - Cambridge Neuroimaging 2010

function aacloudclient()

% Check Java
if ~usejava('jvm')
    error([mfilename ' requires Java to run.']);
end;

if ~usejava('swing')
    error([mfilename ' requires swing Java component to run (email Rhodri if this is a problem).']);
end;

% See whether it is already running
loadclass=true;
try
    import camneuro.aacc.*;
    tmp=Aacc();
    if isa(tmp,'camneuro.aacc.Aacc')
        loadclass=false;
    end;
catch
end;

% Only load if not already available
if (loadclass)
    [pth nme ext]=fileparts(mfilename('fullpath'));
    
    javaaddpath(fullfile(pth,'aacc.jar'));
    
    libpth=fullfile(pth,'lib');
    
    fn=dir(fullfile(libpth,'*.jar'));
    
    for ind=1:length(fn)
        javaaddpath(fullfile(libpth,fn(ind).name));
    end;
    
    try
        tmp=camneuro.aacc.Aacc();
    catch
        fprintf('There was a problem setting up the cloud client.\n');
    end;
end;
