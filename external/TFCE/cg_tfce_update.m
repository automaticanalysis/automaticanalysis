function cg_tfce_update(update)
% check for new updates
%
% FORMAT cg_tfce_update(update)
% update - allow installation of update
% 
% This function will connect itself to the SBM server, compare the
% version number of the updates with the one of the VBM8 installation 
% currently in the MATLAB path and will display the outcome.
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_tfce_update.m 50 2011-06-21 15:01:54Z gaser $

rev = '$Rev: 50 $';

if nargin == 0
  update = 0;
end

r = 0;

% get current release number
A = ver;
for i=1:length(A)
  if strcmp(A(i).Name,'TFCE Toolbox')
    r = str2double(A(i).Version);
  end
end

url = 'http://dbm.neuro.uni-jena.de/tfce/';

% get new release number
if usejava('jvm')
  [s,sts] = urlread(url);
  if ~sts
    fprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\n. You can download your update at %s\n',url,url); 
    return
  end
else
  fprintf('Please enable Java (JVM) to use update function.\n. You can download your update at %s\n',url); 
  return
end

n = regexp(s,'tfce_r(\d.*?)\.zip','tokens');
if isempty(n)
  fprintf('There are no new releases available yet.\n');
  return;
else
  % get largest release number
  rnew = [];
  for i=1:length(n)
    rnew = [rnew str2double(n{i})];
  end
  rnew = max(rnew);
end

if rnew > r
  fprintf('A new version of TFCE is available on: %s\n',url);
  fprintf('Your version: %d - New version: %d\n',r,rnew);
  if ~update
    fprintf('In order to update use Toolbox|tfce|Check for updates\n',r,rnew);
  end

  if update
    d = fullfile(spm('Dir'),'toolbox'); 
    overwrite = spm_input('Update',1,'m','Do not update|Download zip-file only|Overwrite old TFCE installation',[-1 0 1],3);
    switch overwrite
    case 1
      try
        % list mex-files and delete these files to prevent that old
        % compiled files are used
        mexfiles = dir(fullfile(d,'tfce','*.mex*'));
        for i=1:length(mexfiles)
          name = fullfile(d,'tfce',mexfiles(i).name);
          spm_unlink(name);
        end
        fprintf('Download TFCE\n');
        s = unzip([url sprintf('tfce_r%d.zip',rnew)], d);
        fprintf('%d files have been updated.\nSPM should be restarted.\n',numel(s));
        restart = spm_input('Restart SPM',1,'m','no|yes',[0 1],2);
        if restart
          rehash
          toolbox_path_cache
          eval(['spm fmri']);
        end
      catch
        fprintf('Update failed: check file permissions. Download zip-file only.\n');
        web([url sprintf('tfce_r%d.zip',rnew)],'-browser');
        fprintf('Unzip file to %s\n',d);
      end
    case 0
      web([url sprintf('tfce_r%d.zip',rnew)],'-browser');
      fprintf('Unzip file to %s\n',d);
    end
  end
elseif update
  fprintf('You already have the newest version %d.\n',r);
end