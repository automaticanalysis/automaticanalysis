% aacloudclient launch script
%  Rhodri Cusack - Cambridge Neuroimaging 2010

function aacloudclient()
global cc

% Check Java
if ~usejava('jvm')
    error([mfilename ' requires Java to run.']);
end;

if ~usejava('swing')
    error([mfilename ' requires swing Java component to run (email Rhodri if this is a problem).']);
end;

% See whether it is already running
if ~isa(cc,'camneuro.aacc.Aacc')
    import camneuro.aacc.*;
    [pth nme ext]=fileparts(mfilename('fullpath'));
    
    javaaddpath(fullfile(pth,'aacc.jar'));
    
    libpth=fullfile(pth,'lib');
    
    fn=dir(fullfile(libpth,'*.jar'));
    
    for ind=1:length(fn)
        javaaddpath(fullfile(libpth,fn(ind).name));
    end;
    cc=camneuro.aacc.Aacc();
end;
