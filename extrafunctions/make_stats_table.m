function make_stats_table(varargin)
%
% display and optionally save to jpeg an SPM thresholded stats table
%
% Usage:
%
%  make_stats_table(SPM, [fname, Ic, u, thresDesc, k, Im]);
%
% INPUT
%
%	SPM			- SPM struct containing thresholded results
%
% optional
%
%	fname		- jpeg name or path (omit or pass in [] to skip save)
%	Ic			- contrast # to examine
%	u			- threshold
%	thresDesc	- thresholding ('none' or 'FWE')
%	k			- cluster extent
%	Im			- masking option (currently ignored)
%
% if an optional parameter is missing, the following defaults are used:
%
%	no save (fname = [])
%	contrast #1 (Ic = 1)
%	significance = 0.001 (u = 0.001)
%	uncorrected (thresDesc = 'none')
%	voxel extent (k) = 0
%	no mask (Im = [])
%		
% sanity check

if (nargin < 1)
	fprintf('Usage: make_stats_table(SPM, [fname, Ic, u, thresDesc, k, Im])\n');
	return;
end

SPM = varargin{1};

% sanity check -- make sure the passed SPM has results to display 

if (~isfield(SPM,'xCon') || isempty(SPM.xCon))
	fprintf('%s: SPM struct does not contain contrasts. Aborting...\n', mfilename);
	return;
end

% defaults

fname = [];
SPM.Ic = 1;
SPM.u = 0.001;
SPM.thresDesc = 'none';
% 	SPM.u = 0.05;
% 	SPM.thresDesc = 'FWE';
SPM.k = 0;
SPM.Im = [];


if (nargin > 1)
	fname = varargin{2};
end

if (nargin > 2)
	SPM.Ic = varargin{3};
end

if (nargin > 3)
	SPM.u = varargin{4};
end

if (nargin > 4)
	SPM.thresDesc = varargin{5};
end

if (nargin > 5)
	SPM.k = varargin{6};
end

% Im is a ruse: currently we can only handle "no mask"

% 	if (nargin > 6)
%  		SPM.Im = varargin{7};
% 	end


% spm_list_display_noUI can only handle 'none' and 'FWE'

if (~strcmp(SPM.thresDesc,'none') && ~strcmp(SPM.thresDesc,'FWE'))
	fprintf('%s: Can only handle ''none'' and ''FWE'' thresholding. Aborting...\n', mfilename);
	return;
end

% extract info and make table

[~,xSPM] = spm_getSPM(SPM);
TabDat = spm_list('Table',xSPM);
hreg = spm_figure('CreateSatWin');
spm_list_display_noUI(TabDat,hreg);

% save?

if (~isempty(fname))

	set(hreg,'Renderer','opengl'); % I think this is the default

	% workaround for font rescaling weirdness

	set(findall(hreg,'Type','text'),'FontUnits','normalized');

	% tweak paper size to account for landscape layout

	set(hreg,'PaperUnits','inches','PaperPosition',[0 0 5 4]);

	% force jpg suffix otherwise there can be weirdness
	% also, bump up resolution to 200 bc numbers

	[p,n,~] = fileparts(fname);
	fname = fullfile(p,[n '.jpg']);

	print(hreg, '-djpeg', '-r200', fname);

	close(hreg);

end

end