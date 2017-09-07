% Automatic analysis - check versions of MATLAB and aa
% specify min and max version numbers for aa

function aas_requiresversion(aap)

MATLAB_minver = {'7.14' 'R2012a'}; % (hardcoded)
global aa

%% MATLAB
currver = sscanf(version,'%d.%d%*s')'*[100 1]'; % only main and secondary version
if currver>=sscanf(MATLAB_minver{1},'%d.%d%*s')'*[100 1]'
    aas_log(aap,0,sprintf('INFO: Current MATLAB version %s suitable for this aa %s %s\n',version, aa.Version, aa.Date));
else
    aas_log(aap,1,sprintf('MATLAB version %s earlier than required %s (%s) for aa %s %s\n',version,MATLAB_minver{1},MATLAB_minver{2},aa.Version,aa.Date));
end

%% aa
currver = sscanf(aa.Version,'%d.%d.%d%*s')'*[100^2 100 1]';
minver  = sscanf(aap.options.aa_minver,'%d.%d.%d')'*[100^2 100 1]';
maxver  = sscanf(aap.options.aa_maxver,'%d.%d.%d')'*[100^2 100 1]';

if (currver>=minver && currver<=maxver) 
    aas_log(aap,0,sprintf('INFO: Current aa version %s %s suitable for this user script\n',aa.Version,aa.Date));
else
    aas_log(aap,1,sprintf('aa version %s outside of range %s to %s suitable for user script\n',aa.Version,aap.options.aa_minver,aap.options.aa_maxver));
end

end
