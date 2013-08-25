function long = cg_vbm8_long
% Configuration file for long correction between an image pair
%
% Christian Gaser
% $Id: cg_vbm8_longitudinal_multi.m 404 2011-04-11 10:03:40Z gaser $

mov = cfg_files;
mov.name = 'Longitudinal data for one subject';
mov.tag  = 'mov';
mov.filter = 'image';
mov.num  = [1 Inf];
mov.help   = {[...
'These are the data of the same subject.']};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {[...
'Images of the same subject.']};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj };
esubjs.num     = [1 Inf];
esubjs.help = {[...
'Specify data for each subject.']};

%------------------------------------------------------------------------

long = cfg_exbranch;
long.name = 'Process longitudinal data';
long.tag  = 'long';
long.val  = {esubjs};
long.prog = @cg_vbm8_longitudinal_multi_run;
long.vout = @vout_long;
long.help = {
'This option provides customized processing of longitudinal data.'};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_long(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Segmented longitudinal data (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
%------------------------------------------------------------------------
