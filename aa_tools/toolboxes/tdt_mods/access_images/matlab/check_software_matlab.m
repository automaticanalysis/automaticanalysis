function checked = check_software_matlab(software)

if length(software)>=6 && strcmpi(software(1:6),'matlab')
    try
        matlab_ver = ver;
    catch %#ok<CTCH>
        error('MATLAB is not appropriately configured or the path is corrupt!\nPlease run restoredefaultpath or change cfg.software to another software')
    end
    assert(isstruct(matlab_ver) && isfield(matlab_ver,'Name') && any(strcmpi({matlab_ver.Name},'MATLAB')), 'invalid version information');
    matlab_ver = matlab_ver(strcmpi({matlab_ver.Name},'MATLAB'));
    if length(software)>6
        softwarever = software(7:end);
        matlabver = regexp(matlab_ver.Release,'[rR]{1}[0-9]{4}[abAB]{1}','match');
        assert(strcmpi(softwarever,matlabver),sprintf('cfg.software = %s, but MATLAB version = %s\nPlease decide what you want to use.',softwarever,matlabver{1}))
    end
    % seems to be fine, so return this
    checked = true;
else 
    % software does not start with MATLAB, should not happen here
    error('cfg.software does not start with MATLAB. This should not happen here')
end