function cg_vbm8_batch(namefile,writeonly,vbm8_defaults)
% wrapper for using batch mode (see cg_vbm8_batch.sh)
%
% namefile      - array of file names
% writeonly     - if "1" do not estimate segmentations
% vbm8_defaults - use this default file instead of cg_vbm8_defaults.m
%
%_______________________________________________________________________
% $Id: cg_vbm8_batch.m 404 2011-04-11 10:03:40Z gaser $

if nargin < 1
	fprintf('Syntax: cg_vbm8_batch(namefile)\n');
	return
end

if nargin < 2
	writeonly = 0;
end

spm_get_defaults;

if nargin < 3
    cg_vbm8_defaults;
else
    if isempty(vbm8_defaults)
        cg_vbm8_defaults;
    else
        fprintf('Use defaults in %s.\n',vbm8_defaults);
        [path, name] = fileparts(vbm8_defaults);
        oldpath = pwd;
        eval(['cd ' path])
        eval(name);
        eval(['cd ' oldpath])
    end
end
global defaults vbm8

% always deselect print option
vbm8.extopts.print = 0;

names = textread(namefile,'%s');
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end

if writeonly
	matlabbatch{1}.spm.tools.vbm8.write = vbm8;
else
	matlabbatch{1}.spm.tools.vbm8.estwrite = vbm8;
end

for i=1:n
	if writeonly
		matlabbatch{1}.spm.tools.vbm8.write.data{i} = names{i};
	else
		matlabbatch{1}.spm.tools.vbm8.estwrite.data{i} = names{i};
	end
end

try
  spm_jobman('initcfg');
  spm_jobman('run_nogui',matlabbatch);
end

spm_unlink(char(namefile))

exit
