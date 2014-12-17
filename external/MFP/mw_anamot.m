function mw_anamot(mydata, silent);
%
% Little companion script that reads out result files from
% mw_mfp (so run this first) and summarizes motion in the form of
% total and scan-to-scan displacement at average cortical distance.
%
% To run, type 'mw_anamot' in the matlab window and select
% the directory that contains the mw_motion.mat files.
% Alternatively, pass a directory name and the optional "silent"
% parameter to suppress screen output.
%
% If you specify a directory, all sub-directories will be searched
% for mw_mfp result files and analyzed; results are written to a 
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
% as well as
%
%   Wilke M: Isolated assessment of translation or rotation
%   severely underestimates the effects of subject motion
%   in fMRI data;
%   PLoS ONE 2014, 9: e106498
% 
% freely available at http://www.dx.doi.org/10.1371/journal.pone.0106498
%
% It requires the following functions to run properly:
%
%   subdir (included, courtesy of the Mathworks File Exchange)
%   
% Idea and implementation by Marko Wilke, see file for version information.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details. You should have received a copy of the GNU General
% Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
%


% ==========================================================================================================
%                                          Check, set options
% ==========================================================================================================


% File version: 1.3, 2014-11-13


% get an nice working environment
  ori = pwd;


% process inputs
  if nargin == 0

	mydata = spm_select(1,'dir','Please select data directory to search',[],pwd);
  end;
  if nargin < 2 || isempty(silent),    silent = 0;   end;


% say hello
  if silent == 0,  clc;  disp(['... Welcome to ' mfilename]);  end;


% search this folder/file location
  if ~isdir(mydata),  [mydata, ~,~,~] = spm_fileparts(mydata);  if isempty(mydata),  mydata = pwd;  end;  end;
  cd(mydata);
  myparams = subdir(['mw_motion.mat']);
  if ~isempty(myparams)

	  if silent == 0,  disp(['... found ' num2str(size(myparams,1)) ' sessions, please wait...']);  end;
  else
	  if silent == 0,  error(['... Sorry, no motion parameters (mw_motion.mat) were found, aborting!']);
	  else,            disp(['... Sorry, no motion parameters (mw_motion.mat) were found, aborting!']);  return;
	  end;
  end;


% settings
  thrs = [1 2 3 4 5];    % assess these values [in mm] as indicators for total discplacement
  thrs_s = thrs./2;      % assess these values [in mm] as indicators for scan-to-scan discplacement
  dvs = 3;               % if automatic detection did not work, assume voxel size to be this value


% take a note
  outfile = [mydata filesep 'mw_anamot_results.txt'];
  if exist(outfile) ~= 2

	% create new file
	  if silent == 0,  disp(['... generating new mw_anamot results file, please wait...']);  end;
	  fid = fopen(outfile,'At+');
	  fprintf(fid, ['Source path' '\t' 'TD: mean' '\t' 'TD: STD' '\t' 'TD: median' '\t' 'TD: min' '\t' 'TD: max' '\t' ['TD:#>' num2str(thrs(1,1)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,2)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,3)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,4)) 'mm'] '\t' ['TD:#>' num2str(thrs(1,5)) 'mm'] '\t' 'STS: mean' '\t' 'STS: STD' '\t' 'STS: median' '\t' 'STS: min' '\t' 'STS: max' '\t' ['STS:#>' num2str(thrs_s(1,1)) 'mm'] '\t' ['STS:#>' num2str(thrs_s(1,2)) 'mm'] '\t' ['STS:#>' num2str(thrs_s(1,3)) 'mm'] '\t' ['STS:#>' num2str(thrs_s(1,4)) 'mm'] '\t' ['STS:#>' num2str(thrs_s(1,5)) 'mm'] '\t' 'Mean Voxel size [mm]' '\n']);

  else

	% if file exists already, no need to include headers again
	  if silent == 0,  disp(['... adding to existing  mw_anamot results file, please wait...']);  end;
	  fid = fopen(outfile,'At+');

  end;


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
		  vs = sqrt(sum(V.mat(1:3,1:3).^2));
		  vs = [sprintf('%0.1f', vs(1)) ' * ' sprintf('%0.1f', vs(2)) ' * ' sprintf('%0.1f', vs(3))];

	  catch,

		  vs = 'N/A';

	  end;


	% check if data is available
	  try,  sum(td);   catch,  td  = NaN;  end;
	  try,  sum(sts);  catch,  sts = NaN;  end;


	% simplified output to screen
	  if silent == 0

		  disp(['   ... session ' num2str(i) ' from ' p ' (with a voxel size of ' vs ' mm):']);

		  if ~isnan(td)

			disp(['       ... total displacement mean = ' sprintf('%0.2f', mean(td)) ' [' sprintf('%0.2f', std(td)) '], median = ' sprintf('%0.2f', median(td)) ' [' sprintf('%0.2f', min(td)) '-' sprintf('%0.2f', max(td)) ']']);

		  else 

			disp(['       ... total displacement data not available!']);

		  end;
		  if ~isnan(sts)

			disp(['       ... scan-to-scan displacement mean = ' sprintf('%0.2f', mean(sts)) ' [' sprintf('%0.2f', std(sts)) '], median = ' sprintf('%0.2f', median(sts)) ' [' sprintf('%0.2f', min(sts)) '-' sprintf('%0.2f', max(sts)) ']']);

		  else 

			disp(['       ... scan-to-scan displacement data not available!']);

		  end;

	  end;


	% more comprehensive output to file
	  p = strrep(p,'\','/');
	  fprintf(fid, [p '\t' sprintf('%0.2f', mean(td)) '\t' sprintf('%0.2f', std(td)) '\t' sprintf('%0.2f', median(td)) '\t' sprintf('%0.2f', min(td)) '\t' sprintf('%0.2f', max(td)) '\t' num2str(sum(td>thrs(1,1))) '\t' num2str(sum(td>thrs(1,2))) '\t' num2str(sum(td>thrs(1,3))) '\t' num2str(sum(td>thrs(1,4))) '\t' num2str(sum(td>thrs(1,5))) '\t' sprintf('%0.2f', mean(sts)) '\t' sprintf('%0.2f', std(sts)) '\t' sprintf('%0.2f', median(sts)) '\t' sprintf('%0.2f', min(sts)) '\t' sprintf('%0.2f', max(sts)) '\t'  num2str(sum(sts>thrs_s(1,1))) '\t'   num2str(sum(sts>thrs_s(1,2))) '\t' num2str(sum(sts>thrs_s(1,3))) '\t' num2str(sum(sts>thrs_s(1,4))) '\t' num2str(sum(sts>thrs_s(1,5)))  '\t'  vs '\n']);

  end;


% housekeeping
  if silent == 0

	  disp(['... done with analyzing ' num2str(size(myparams,1)) ' sessions;']);
	  disp(['... results were stored in ' mydata]);

  end;


% close output file, get back
  fclose(fid);
  cd(ori);
  return;


% ==========================================================================================================
%                                          Embedded function 1: subdir
% ==========================================================================================================

function varargout = subdir(varargin)
% SUBDIR Performs a recursive file search
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
