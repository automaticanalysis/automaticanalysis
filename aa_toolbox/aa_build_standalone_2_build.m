function aa_build_standalone_2_build(config, outdir)
aa = aaClass('nopath','nogreet');
if isstruct(config) % aap (run-time compile)
    aap = config;
elseif ischar(config) % aap_parameters_defaults file (offline standalone)
    aap=xml_read(config,struct('ReadAttr',0));
end

% MATLAB Toolboxes: 
%   Statistics_Toolbox
%   Parallel_Computing_Toolbox
%   Image_Processing_Toolbox
tbx = {'-p',fullfile(matlabroot,'toolbox','stats'),...
    '-p',fullfile(matlabroot,'toolbox','images'),...
    '-p',fullfile(matlabroot,'toolbox','distcomp')};

% aa
aadir = cellstr(char(spm_select('FPList',aa.Path,'dir','.*'),spm_select('FPList',aa.Path,'.*')));
aadir = mcc_a(aadir);

% SPM
prepare_spm(aap);
spmtoolsdir = textscan(aap.directory_conventions.spmtoolsdir,'%s','delimiter',':'); 
spmtoolsdir = mcc_a(spmtoolsdir{1});

aas_makedir(aap,outdir);
mkdir(fullfile(outdir,[aa.Name strtok(aa.Version,'.')]));

mcc('-m', '-C', '-v',...
    '-o',aa.Name,...
    '-d',fullfile(outdir,[aa.Name strtok(aa.Version,'.')]),...
    '-N',tbx{:},...
    '-R','-singleCompThread',...
    aadir{:},...
    '-a',aap.directory_conventions.spmdir,...
    spmtoolsdir{:},...
    '-a',aap.directory_conventions.eeglabdir,...
    '-a',aap.directory_conventions.GIFTdir,...
    '-a',aap.directory_conventions.BrainWaveletdir,...
    '-a',aap.directory_conventions.ANTSdir,...
    '-a',aap.directory_conventions.DCMTKdir,...
    '-a',aap.directory_conventions.templatedir,...
    'aa_standalone.m');

end

function list_a = mcc_a(list)
list_a(2:2:numel(list)*2) = list(:);
list_a(1:2:numel(list_a)) = cellstr('-a');
end

function prepare_spm(aap)
% Prepare spm for deployment with aa
% Based on Guillaume Flandin's spm_make_standalone.m v6416

spmdir = aap.directory_conventions.spmdir;

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spmdir,'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% create code to insert toolbox config
%-Toolbox autodetection
%-Get the list of toolbox directories
tbxdir = fullfile(spmdir,'toolbox');
d  = dir(tbxdir); d = {d([d.isdir]).name};
dd = regexp(d,'^\.');
%(Beware, regexp returns an array if input cell array is of dim 0 or 1)
if ~iscell(dd), dd = {dd}; end
d  = {'' d{cellfun('isempty',dd)}};
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:length(d)
    d2 = fullfile(tbxdir,d{i});
    di = dir(d2); di = {di(~[di.isdir]).name};
    f2 = regexp(di,'.*_cfg_.*\.m$');
    if ~iscell(f2), f2 = {f2}; end
    fi = di(~cellfun('isempty',f2));
    if ~isempty(fi)
        ft = [ft(:); fi(:)];
    end
end
if ~isempty(ft)
    if isempty(ft)
        ftstr = '';
    else
        ft = cellfun(@(cft)strtok(cft,'.'),ft,'UniformOutput',false);
        ftstr  = sprintf('%s ', ft{:});
    end
    fprintf(fid,'values = {%s};\n', ftstr);
end
fclose(fid);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
%-Duplicate Contents.m in Contents.txt for use in spm('Ver')
%==========================================================================
sts = copyfile(fullfile(spm('Dir'),'Contents.m'),...
               fullfile(spm('Dir'),'Contents.txt'));
if ~sts, warning('Copy of Contents.m failed.'); end

end