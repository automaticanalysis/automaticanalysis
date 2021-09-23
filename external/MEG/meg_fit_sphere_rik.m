function spherefit = meg_fit_sphere_rik(aap,infname,outfname)

if nargin < 2
    outfname = infname;
end

% Assumes FIFACCESS tool hpipoints and Neuromag utility
% fit_sphere_to_points are on path

%% Read digitized points from FIF file
[pth fstem ext] = fileparts(outfname);
    
hptxtfname = fullfile(pth,[fstem '_headpoints.txt']); % specifying head points text file
if exist(hptxtfname,'file')~=2
    [co ki] = hpipoints(infname); % HPI points  
    if size(co,2)~=3; co=co'; end % In case linux command above transposes on different machines (only fail if 3 headpoints too!)
    headpoints = co(ki>1,:);  % don't include the fiducial points
    headpoints = headpoints(find(~(headpoints(:,2)>0 & headpoints(:,3)<0)),:); % Remove nose points:
    save(hptxtfname,'-ASCII','headpoints');
end

%% Fitting sphere to headshape points
sphtxtfname = fullfile(pth,[fstem '_sphere_fit.txt']);
if exist(sphtxtfname,'file')~=2
    cmd_fit = [fullfile(aap.directory_conventions.neuromagdir,'bin','util','fit_sphere_to_points') ' ' hptxtfname];
    [status, out] = unix(cmd_fit);
    spherefit=regexp(out,'[0-9]{1}\.[0-9]{6}','match');
    if status ~= 0 || isempty(spherefit)
        error('Spherefit failed! - %s',out)
    end    
    spherefit = str2double(spherefit)*1000;
    save(sphtxtfname,'-ASCII','spherefit');
else
    spherefit = textread(sphtxtfname,'%f');
end

return