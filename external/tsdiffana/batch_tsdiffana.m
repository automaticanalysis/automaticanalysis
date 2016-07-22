%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 301 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.name = 'Image Series';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.tag = 'imgs';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.filter = 'image';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.dir = '';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.num = [1 Inf];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.check = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.help = {'Select a time series of images.'};
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.def = [];
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.type = 'cfg_repeat';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.name = 'Sessions/Subjects';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.tag = 'imgs_unused';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1) = cfg_dep;
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tname = 'Values Item';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tgt_spec = {};
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).sname = 'Files: Image Series (cfg_files)';
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_output = substruct('()',{1});
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.num = [1 Inf];
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.forcestruct = false;
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.check = [];
matlabbatch{2}.menu_cfg{1}.menu_struct{1}.conf_repeat.help = {'Enter one or more time series. Each series will be processed independently.'};
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.type = 'cfg_menu';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.name = 'Create Difference Images';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.tag = 'vf';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.labels = {
                                                             'No'
                                                             'Yes'
                                                             }';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.values = {
                                                             false
                                                             true
                                                             }';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.check = [];
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.help = {'Select whether difference images should be written to disk.'};
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_menu.def = @(val)run_tsdiffana('defaults','vf',val{:});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'Analyse Time Series';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'tsdiffana_timediff';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'Repeat: Sessions/Subjects (cfg_repeat)';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).sname = 'Menu: Create Difference Images (cfg_menu)';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @(job)run_tsdiffana('run','timediff',job);
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = @(job)run_tsdiffana('vout','timediff',job);
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {'Run time series diagnostics.'};
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.type = 'cfg_files';
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.name = 'Timeseries Analysis Data Files';
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.tag = 'tdfn';
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.filter = 'mat';
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.ufilter = '^timediff.mat$';
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.dir = '';
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.num = [1 Inf];
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.check = [];
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.help = {'Select one or more files. If more than one file is selected, each report will be displayed on a separate page of the SPM graphics window.'};
matlabbatch{5}.menu_cfg{1}.menu_entry{1}.conf_files.def = [];
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.type = 'cfg_menu';
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.name = 'Print to File';
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.tag = 'doprint';
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.labels = {
                                                             'Yes'
                                                             'No'
                                                             }';
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.values = {
                                                             true
                                                             false
                                                             }';
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.check = [];
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.help = {};
matlabbatch{6}.menu_cfg{1}.menu_entry{1}.conf_menu.def = @(val)run_tsdiffana('defaults','doprint',val{:});
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'Plot Analysis Results';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'tsdiffana_tsdiffplot';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'Files: Timeseries Analysis Data Files (cfg_files)';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).sname = 'Menu: Print to File (cfg_menu)';
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @(job)run_tsdiffana('run','tsdiffplot',job);
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = [];
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{7}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {'Display result from time series analyses.'};
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.type = 'cfg_repeat';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.name = 'TSDiffAna';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.tag = 'tsdiffana_tools';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1) = cfg_dep;
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tname = 'Values Item';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tgt_spec = {};
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).sname = 'Exbranch: Analyse Time Series (cfg_exbranch)';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_output = substruct('()',{1});
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{2}(1) = cfg_dep;
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{2}(1).tname = 'Values Item';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{2}(1).tgt_spec = {};
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{2}(1).sname = 'Exbranch: Plot Analysis Results (cfg_exbranch)';
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{2}(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{2}(1).src_output = substruct('()',{1});
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.num = [1 Inf];
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.forcestruct = true;
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.check = [];
matlabbatch{8}.menu_cfg{1}.menu_struct{1}.conf_repeat.help = {};
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_fname = 'tbx_cfg_tsdiffana.m';
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_dir = {'/afs/fbi.ukl.uni-freiburg.de/home/vglauche/matlab/tsdiffana/'};
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_var(1) = cfg_dep;
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_var(1).tname = 'Root node of config';
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_var(1).tgt_spec = {};
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_var(1).sname = 'Repeat: TSDiffAna (cfg_repeat)';
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_var(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{9}.menu_cfg{1}.gencode_gen.gencode_var(1).src_output = substruct('()',{1});
