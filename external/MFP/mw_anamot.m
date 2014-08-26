function mw_anamot
%
% Little companion script that reads out result files from
% mw_mfp (so run this first) and summarizes motion in the form of
% total and scan-to-scan displacement at average cortical distance.
%
% To run, type 'mw_anamot' in the matlab window and select
% the directory that contains the mw_motion.mat files.
%
% If you specify a top-level directory, all sub-directories
% will be searched and analyzed; results are written to a 
% mw_anamot_results.txt file in this directory (so make sure
% you have write permissions). Additionally, a short summary
% for each session is printed on the screen.
%
% For further information, please refer to and cite 
%
%   Wilke M: An alternative approach towards assessing and
%   accounting for individual motion in fMRI timeseries;
%   NeuroImage 2012, 59: 2062-2072
% 
% available at http://dx.doi.org/10.1016/j.neuroimage.2011.10.043
%
% It requires the following functions to run properly:
%
%   subdir (included, courtesy of the Mathworks File Exchange)
%   
% Idea and implementation by Marko Wilke, see file for version information.
%
% This software comes with absolutely *no* guarantees whatsoever,
% explicit or implied or otherwise, so use it at your own risk!
%
% V1: version at first inclusion in the mw_mfp-package


% ==========================================================================================================
%                                          Check, set options
% ==========================================================================================================


% get an nice working environment
  clc;
  ori = pwd;


% search where?
  mydata = spm_select(1,'dir','Please select data directory to search',[],pwd);
  cd(mydata);


% search this folder
  myparams = subdir(['mw_motion.mat']);
  if ~isempty(myparams)

	  disp(['... found ' num2str(size(myparams,1)) ' sessions, please wait...']);
  else
	  error(['... Sorry, no motion parameters (mw_motion.mat) were found, aborting!']);
  end;


% settings
  thrs = [1 2 3 4 5];    % assess these values [in mm] as indicators for total discplacement
  thrs_s = thrs./2;      % assess these values [in mm] as indicators for scan-to-scan discplacement
  dvs = 3;               % if automatic detection did not work, assume voxel size to be this value


% take a note
  outfile = [mydata filesep 'mw_anamot_results.txt'];
  fid = fopen(outfile,'At+');
  fprintf(fid, ['Source path' '\t' ['TD:#>' num2str(thrs(1,1)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,2)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,3)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,4)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,5)) 'mm'] '\t'  ['STS:#>' num2str(thrs_s(1,1)) 'mm'] '\t'   ['STS:#>' num2str(thrs_s(1,2)) 'mm'] '\t'   ['STS:#>' num2str(thrs_s(1,3)) 'mm'] '\t'   ['STS:#>' num2str(thrs_s(1,4)) 'mm'] '\t'   ['STS:#>' num2str(thrs_s(1,5)) 'mm'] '\n']);


% ==========================================================================================================
%                                          Get going...
% ==========================================================================================================

% loop over parameters
  for i = 1:size(myparams,1)


	% get, load data
	  curr = myparams(i).name;
	  load(curr);


	% get fMRI data in order to get voxel sizes
	  [p nm e v] = spm_fileparts(curr);
	  try,

		  currp = spm_select('list', p, 'rp_.*.txt');
		  currp = spm_select('List', p, ['^' currp(4:end-4) '.(img|nii)']);
		  V = spm_vol([p filesep currp]);
		  vs = (abs(V.mat(1,1))+abs(V.mat(2,2))+abs(V.mat(3,3))) / 3;

	  catch,

		  vs = dvs;

	  end;


	% check if data is available
	  try,  sum(td);   catch,  td  = NaN;  end;
	  try,  sum(sts);  catch,  sts = NaN;  end;


	% simplified output to screen
	  if max(td) > vs

		disp(['   ... session ' num2str(i) ' from ' p ': npoints > vs = ' num2str(sum(td>vs)) ', max: ' num2str(max(td)) ]);

	  elseif max(td) < vs

		disp(['   ... session ' num2str(i) ' from ' p ': no points exceed voxel size!']);

	  elseif isnan(td)

		disp(['   ... session ' num2str(i) ' from ' p ': data not available!']);

	  end;


	% more comprehensive output to file
	  p = strrep(p,'\','/');
	  fprintf(fid, [p '\t' num2str(sum(td>thrs(1,1))) '\t' num2str(sum(td>thrs(1,2))) '\t' num2str(sum(td>thrs(1,3))) '\t' num2str(sum(td>thrs(1,4))) '\t' num2str(sum(td>thrs(1,5))) '\t'  num2str(sum(sts>thrs_s(1,1))) '\t'   num2str(sum(sts>thrs_s(1,2))) '\t' num2str(sum(sts>thrs_s(1,3))) '\t' num2str(sum(sts>thrs_s(1,4))) '\t' num2str(sum(sts>thrs_s(1,5))) '\n']);

  end;


% housekeeping
  disp(['... done with analyzing ' num2str(size(myparams,1)) ' sessions;']);
  disp(['... results were stored in ' mydata]);


% close output file, get back
  fclose(fid);
  cd(ori);
  return;


% ==========================================================================================================
%                                          Embedded function 1: subdir
% ==========================================================================================================

function varargout = subdir(varargin)
%SUBDIR Performs a recursive file search
%
% subdir
% subdir(name)
% files = subdir(...)
%
% This function performs a recursive file search.  The input and output
% format is identical to the dir function.
%
% Input variables:
%
%   name:   pathname or filename for search, can be absolute or relative
%           and wildcards (*) are allowed.  If ommitted, the files in the
%           current working directory and its child folders are returned    
%
% Output variables:
%
%   files:  m x 1 structure with the following fields:
%           name:   full filename
%           date:   modification date timestamp
%           bytes:  number of bytes allocated to the file
%           isdir:  1 if name is a directory; 0 if no
%
% Example:
%
%   >> a = subdir(fullfile(matlabroot, 'toolbox', 'matlab', '*.mat'))
%
%   a = 
%
%   67x1 struct array with fields:
%       name
%       date
%       bytes
%       isdir
%
%   >> a(2)
%
%   ans = 
%
%        name: '/Applications/MATLAB73/toolbox/matlab/audiovideo/chirp.mat'
%        date: '14-Mar-2004 07:31:48'
%       bytes: 25276
%       isdir: 0
%
% See also:
%
%   dir
% Copyright (c) 2009, Kelly Kearney
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
%    * Redistributions of source code must retain the above copyright 
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

%---------------------------
% Get folder and filter
%---------------------------

error(nargchk(0,1,nargin));
error(nargoutchk(0,1,nargout));

if nargin == 0
    folder = pwd;
    filter = '*';
else
    [folder, name, ext] = fileparts(varargin{1});
    if isempty(folder)
        folder = pwd;
    end
    if isempty(ext)
        if isdir(fullfile(folder, name))
            folder = fullfile(folder, name);
            filter = '*';
        end
    else
        filter = [name ext];
    end
end

%---------------------------
% Search all folders
%---------------------------

pathstr = genpath(folder);
seplocs = findstr(pathstr, pathsep);
loc1 = [1 seplocs(1:end-1)+1];
loc2 = seplocs(1:end)-1;
pathfolders = arrayfun(@(a,b) pathstr(a:b), loc1, loc2, 'UniformOutput', false);

Files = [];
for ifolder = 1:length(pathfolders)
    NewFiles = dir(fullfile(pathfolders{ifolder}, filter));
    if ~isempty(NewFiles)
        fullnames = cellfun(@(a) fullfile(pathfolders{ifolder}, a), {NewFiles.name}, 'UniformOutput', false); 
        [NewFiles.name] = deal(fullnames{:});
        Files = [Files; NewFiles];
    end
end

%---------------------------
% Output
%---------------------------
    
if nargout == 0
    if ~isempty(Files)
        fprintf('\n');
        fprintf('%s\n', Files.name);
        fprintf('\n');
    end
elseif nargout == 1
    varargout{1} = Files;
end
return;
