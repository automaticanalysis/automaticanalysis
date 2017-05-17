function cfg = tbx_cfg_mfp;

% SPM Configuration file for the mfp toolbox
%
% Based on sample scripts provided by several authors
% Included in the mfp-distribution in the hope that
% it might be helpful to users; 
%
% Implementation by Marko Wilke, see history file for version information.
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
%                                          Preludes: settings, inputs, etc.
% ==========================================================================================================


% File version: 1.2, 2014-11-13


% potentially add path
  if ~isdeployed,

	addpath(fullfile(spm('Dir'),'toolbox','mfp'));

  end


% ===== select parameters or directory =====
  rps            = cfg_files;
  rps.tag        = 'rps';
  rps.name       = 'Select realignment parameters';
  rps.help       = {'Select realignment parameters or directory to parse'};
  rps.filter     = 'any';
  rps.ufilter    = '^rp_.*txt$';
  rps.num        = [1 inf];


% ===== option 1: do_td =====
  td             = cfg_menu;
  td.name        = 'Calculate total displacement';
  td.tag         = 'td';
  td.help        = {'If set, combines translations and rotations into a single, comprehensive indicator of total displacement, using an indicator of average cortical distance (davg). This is shown as a figure and saved in each session''s directory. Defaults to yes, but can be set to no.',...
                 ' ',...
                 'NB: In contrast to earlier versions, a standard (fixed) davg as set in the options will be used by default; using the individually-determined value (legacy behavior :) is still available;'};
  td.labels      = {'No', 'Yes (standard davg)', 'Yes (individual davg)'};
  td.values      = {0 1 -1};
  td.val         = {1};


% ===== option 2: do_mm =====
  mm             = cfg_menu;
  mm.name        = 'Calculate motion mask';
  mm.tag         = 'mm';
  mm.help        = {'If set, generate a motion mask, illustrating the maximum amount of total displacement each voxel experienced as a function of rotation and translation during a given run; this is saved as an image volume (..._motmask.img) in each session''s directory. Defaults to yes, but can be set to no.'};
  mm.labels      = {'No', 'Yes'};
  mm.values      = {0 1};
  mm.val         = {1};


% ===== option 3: do_mfp =====
  mfp            = cfg_menu;
  mfp.name       = 'Calculate motion fingerprint';
  mfp.tag        = 'mfp';
  mfp.help       = {'If set, a motion fingerprint calculated, which should reflect the signal changes occurring in the brain as a function of motion in this individual subject. This generates a (temporary) sub-folder and saves the results in a text file (mw_mfp...txt) in each session''s directory. Defaults to yes, but can be set to no.'};
  mfp.labels     = {'No', 'Yes'};
  mfp.values     = {0 1};
  mfp.val        = {1};


% ===== option 4: keep =====
  keep           = cfg_menu;
  keep.name      = '.Include how many timecourses';
  keep.tag       = 'keep';
  keep.help      = {'If mfp is set, this option is used to define the number of mfp-timecourses to include in the resulting txt-file. While 9 mfp''s will always be determined, not necessarily all may be used later-on. Defaults to 3 but can be set to [1-9]'};
  keep.labels    = {'1', '2', '3', '4', '5', '6', '7', '8', '9'};
  keep.values    = {1 2 3 4 5 6 7 8 9};
  keep.val       = {3};


% ===== option 5: shifted =====
  shifted         = cfg_menu;
  shifted.name    = '.Include shifted timecourses';
  shifted.tag     = 'shifted';
  shifted.help    = {'If mfp is set, this option is used to define whether the shifted (-1 timepoint) timecourses are included in the resulting text file. This is also known as "Volterra-option 1". Defaults to yes but can be set to no'};
  shifted.labels = {'No', 'Yes'};
  shifted.values = {0 1};
  shifted.val    = {1};


% ===== option 6: squared =====
  squared         = cfg_menu;
  squared.name    = '.Include squared timecourses';
  squared.tag     = 'squared';
  squared.help    = {'If mfp is set, this option is used to define whether the squared timecourses are included in the resulting text file. This is also known as "Volterra-option 2". Defaults to no but can be set to yes'};
  squared.labels = {'No', 'Yes'};
  squared.values = {0 1};
  squared.val    = {0};


% ===== option 7: cleanup =====
  cleanup         = cfg_menu;
  cleanup.name    = 'Cleanup';
  cleanup.tag     = 'cleanup';
  cleanup.help    = {'This option determines the output behavior. It defaults to Silent Mode (remove directory, do not show figures or command line output during processing) but can be set to yes (remove directory, but show figure and command line output) and no (do not remove directory [potentially useful for debugging], show figure and command line output).'};
  cleanup.labels = {'No', 'Yes', 'Silent mode'};
  cleanup.values = {0 1 2};
  cleanup.val    = {2};


% ===== option 8: fout =====
  fout           = cfg_entry;
  fout.name      = 'Outfile';
  fout.tag       = 'fout';
  fout.help      = {'Name of additional output file or folder. If specified (provide an output directory path or a full filename), a graphics summary will additionally be saved there. Irrespective of this setting, a mw_motion.png in the same directory as the realignment parameters will always be saved.'};
  fout.strtype   = 's';
  fout.val       = {''};


% ===== option 9: anamot =====
  anamot         = cfg_menu;
  anamot.name    = 'Analyze motion';
  anamot.tag     = 'anamot';
  anamot.help    = {'Optional step: run mw_anamot on the result files to write out summary results for both total displacement and scan-to-scan displacement (useful for quality control purposes!). Results will be stored in the same directory in the form of a text-file (mw_anamot_results.txt). See mw_anamot for more options. Defaults to yes but can be set to no.'};
  anamot.labels  = {'No', 'Yes'};
  anamot.values  = {0 1};
  anamot.val     = {1};


% ==========================================================================================================
%                                          Provide information for toolbox integration
% ==========================================================================================================
  cfg          = cfg_exbranch;
  cfg.tag      = 'mfp_cfg';
  cfg.name     = 'Motion Fingerprint';
  cfg.val      = {rps td mm mfp keep shifted squared cleanup fout anamot};
  cfg.help     = {'Toolbox configuration utility of the motion fingerprint allgorithm.'};
  cfg.prog     = @mfp_run_tbx;
  cfg.modality = {'FMRI' 'PET' 'EEG'};
  cfg.vout     = @mfp_dep_tbx;


% ==========================================================================================================
%                                          Get going
% ==========================================================================================================
  function out = mfp_run_tbx(job)
  [td, sts, mfpfile] = mw_mfp(char(job.rps), job.td, job.mm, job.mfp, job.keep, job.shifted, job.squared, job.cleanup, char(job.fout));
  if job.anamot == 1

	mw_anamot(char(job.rps), 1);
  end;

  if ~isempty(mfpfile)

	out = struct;
	for i = 1:size(mfpfile,1)

		out.sess{i}.mfpfile = cellstr(deblank(mfpfile(i,:)));

	end;

  else

	out = [];

  end;
  return;

% ==========================================================================================================
%                                          Set up dependency
% ==========================================================================================================
  function dep = mfp_dep_tbx(job);
  for i = 1:size(job.rps,2)


	cdep(1)            = cfg_dep;
	cdep(1).sname      = ['Motion fingerprint output file (Session ' num2str(i) ')'];
	cdep(1).src_output = substruct('.','sess', '{}',{i}, '.','mfpfile');

	if i == 1,  dep = cdep;  else,  dep = [dep cdep];  end;

  end;
  return;
