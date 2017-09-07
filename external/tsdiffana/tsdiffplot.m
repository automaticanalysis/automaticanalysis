function h = tsdiffplot(tdfn,fg,flags, varargin)
% tsfiffplot - plots image difference etc info
% FORMAT tsdiffplot(tdfn,fg,flags, varargin)
% 
% tdfn       - time difference file name - mat file with diff parameters
% fg         - figure handle of figure to display in [spm graphics]
% flags      - zero or more of 
%              'r' - display realignment parameters
%              These are either passed as filename in next arg or
%              collected from the same directory as the above .mat file
%              or selected via the GUI if not present
  
if nargin < 1
  tdfn = spm_get(1, 'timediff.mat', 'Select time diff information');
end
if nargin < 2
  fg = [];
end
if isempty(fg)
  fg = spm_figure('GetWin', 'Graphics'); % use SPM figure window
  spm_figure('Clear');
end
if nargin < 3
  flags = '';
end
if isempty(flags)
  flags = ' ';
end

if any(flags == 'r')
  % plot realignment params
  if ~isempty(varargin)  % file name has been passed (maybe)
    fs(1).name = varargin{1};
  else
    % need to get realignment parameter file
    rwcard = 'realignment*.txt';
    [pn fn ext] = fileparts(tdfn);
    % check for realignment file in directory
    fs = dir(fullfile(pn, rwcard));
    if length(fs) > 1 || isempty(fs)
      % ask for realignment param file
      rfn = spm_get([0 1], rwcard, 'Realignment parameters');
      if ~isempty(rfn)
	fs(1).name = rfn;
      end
    end
  end
  if ~isempty(fs)
    % we do have movement parameters
    mparams = spm_load(fs(1).name);
    subpno = 5;
  else
    % we don't
    mparams = [];
    subpno = 4;
  end
else
  subpno = 4;
end

load(tdfn)
imgno = size(slicediff,1)+1;
zno =   size(slicediff,2);
mom = mean(globals);
sslicediff = slicediff/mom;

figure(fg);

h1 = axes('position', [.1 1-.7/subpno .8 .65*1/subpno]);
h2 = plot(2:imgno,td/mom);
axis([0.5 imgno+0.5 -Inf Inf]);
dxt = unique(round(get(h1, 'xtick')));
dxt = dxt(dxt > 1 & dxt <= imgno);
set(h1,'xtick',dxt);
h3 = xlabel('Difference between image and reference (1)');
h4 = ylabel('Scaled variance');
h  = [h1; h2; h3; h4];

h1 = axes('position', [.1 1-1.7/subpno .8 .65*1/subpno]);
h2 = plot(2:imgno,sslicediff, 'x');
axis([0.5 imgno+0.5 -Inf Inf]);
set(h1,'xtick',dxt);
h3 = xlabel('Difference between image and reference (1)');
h4 = ylabel('Slice by slice variance');
h  = [h; h1; h2; h3; h4];

h1 = axes('position', [.1 1-2.7/subpno .8 .65*1/subpno]);
h2 = plot(globals/mom);
axis([0.5 imgno+0.5 -Inf Inf]);
xt = unique(round(get(h1, 'xtick')));
xt = xt(xt > 0 & xt <= imgno);
set(h1,'xtick',xt);
h3 = xlabel('Image number');
h4 = ylabel('Scaled mean voxel intensity');
h  = [h; h1; h2; h3; h4];

h1 = axes('position', [.1 1-3.7/subpno .8 .65*1/subpno]);
mx = max(sslicediff);
mn = min(sslicediff);
avg = mean(sslicediff);
h2 = plot(avg, 'k');
axis([0 zno+1 -Inf Inf]);
hold on
h3 = plot(mn, 'b');
h4 = plot(mx, 'r');
hold off
h5 = xlabel('Slice number');
h6 = ylabel('Slice variance');
h7 = legend('Mean','Min','Max',0);
h  = [h; h1; h2; h3; h4; h5; h6; h7];

% realignment params
if any(flags == 'r')
  h1 = axes('position', [.1 1-4.7/subpno .8 .65*1/subpno]);
  h2 = plot(mparams(:,1:3));
  h3 = legend('x translation','y translation','z translation',0);
  h4 = xlabel('image');
  h5 = ylabel('translations in mm');
  h  = [h; h1; h2; h3; h4; h5];
end

% and label with first image at bottom
cp = get(gca,'Position');
h1 =  axes('Position', [0 0 1 1], 'Visible', 'off');
img1  = deblank(imgs(1,:));
h2 = text(0.5,cp(2)/2.5,{'First image:',img1},'HorizontalAlignment','center');
h  = [h; h1; h2];

return
