function tsdiffana_tools = tbx_cfg_tsdiffana
% 'TSDiffAna' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2008-07-01 17:20:33.
% ---------------------------------------------------------------------
% imgs Image Series
% ---------------------------------------------------------------------
imgs         = cfg_files;
imgs.tag     = 'imgs';
imgs.name    = 'Image Series';
imgs.help    = {'Select a time series of images.'};
imgs.filter = 'image';
imgs.ufilter = '.*';
imgs.num     = [1 Inf];
% ---------------------------------------------------------------------
% imgs_unused Sessions/Subjects
% ---------------------------------------------------------------------
imgs_unused         = cfg_repeat;
imgs_unused.tag     = 'imgs_unused';
imgs_unused.name    = 'Sessions/Subjects';
imgs_unused.help    = {'Enter one or more time series. Each series will be processed independently.'};
imgs_unused.values  = {imgs };
imgs_unused.num     = [1 Inf];
% ---------------------------------------------------------------------
% vf Create Difference Images
% ---------------------------------------------------------------------
vf         = cfg_menu;
vf.tag     = 'vf';
vf.name    = 'Create Difference Images';
vf.help    = {'Select whether difference images should be written to disk.'};
vf.def = @(val)run_tsdiffana('defaults','vf',val{:});
vf.labels = {
             'No'
             'Yes'
             }';
vf.values = {
             false
             true
             }';
% ---------------------------------------------------------------------
% tsdiffana_timediff Analyse Time Series
% ---------------------------------------------------------------------
tsdiffana_timediff         = cfg_exbranch;
tsdiffana_timediff.tag     = 'tsdiffana_timediff';
tsdiffana_timediff.name    = 'Analyse Time Series';
tsdiffana_timediff.val     = {imgs_unused vf };
tsdiffana_timediff.help    = {'Run time series diagnostics.'};
tsdiffana_timediff.prog = @(job)run_tsdiffana('run','timediff',job);
tsdiffana_timediff.vout = @(job)run_tsdiffana('vout','timediff',job);
% ---------------------------------------------------------------------
% tdfn Timeseries Analysis Data Files
% ---------------------------------------------------------------------
tdfn         = cfg_files;
tdfn.tag     = 'tdfn';
tdfn.name    = 'Timeseries Analysis Data Files';
tdfn.help    = {'Select one or more files. If more than one file is selected, each report will be displayed on a separate page of the SPM graphics window.'};
tdfn.filter = 'mat';
tdfn.ufilter = '^timediff.mat$';
tdfn.num     = [1 Inf];
% ---------------------------------------------------------------------
% doprint Print to File
% ---------------------------------------------------------------------
doprint         = cfg_menu;
doprint.tag     = 'doprint';
doprint.name    = 'Print to File';
doprint.def = @(val)run_tsdiffana('defaults','doprint',val{:});
doprint.labels = {
                  'Yes'
                  'No'
                  }';
doprint.values = {
                  true
                  false
                  }';
% ---------------------------------------------------------------------
% tsdiffana_tsdiffplot Plot Analysis Results
% ---------------------------------------------------------------------
tsdiffana_tsdiffplot         = cfg_exbranch;
tsdiffana_tsdiffplot.tag     = 'tsdiffana_tsdiffplot';
tsdiffana_tsdiffplot.name    = 'Plot Analysis Results';
tsdiffana_tsdiffplot.val     = {tdfn doprint };
tsdiffana_tsdiffplot.help    = {'Display result from time series analyses.'};
tsdiffana_tsdiffplot.prog = @(job)run_tsdiffana('run','tsdiffplot',job);
% ---------------------------------------------------------------------
% tsdiffana_tools TSDiffAna
% ---------------------------------------------------------------------
tsdiffana_tools         = cfg_repeat;
tsdiffana_tools.tag     = 'tsdiffana_tools';
tsdiffana_tools.name    = 'TSDiffAna';
tsdiffana_tools.values  = {tsdiffana_timediff tsdiffana_tsdiffplot };
tsdiffana_tools.num     = [1 Inf];
tsdiffana_tools.forcestruct = true;
% ---------------------------------------------------------------------
% add path to this mfile
% ---------------------------------------------------------------------
addpath(fileparts(mfilename('fullpath')));
