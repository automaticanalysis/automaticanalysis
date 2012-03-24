function tsdiffplot(tdfn,fg,flags, varargin)
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
  if length(varargin) > 0  % file name has been passed (maybe)
    fs(1).name = varargin{1};
  else
    % need to get realignment parameter file
    rwcard = 'realignment*.txt';
    [pn fn ext] = fileparts(tdfn);
    % check for realignment file in directory
    fs = dir(fullfile(pn, rwcard));
    if length(fs) > 1 | isempty(fs)
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

subplot(subpno,1,1);
plot(td/mom)
axis([0 imgno -Inf Inf]);
xlabel('Difference image number');
ylabel('Scaled variance');

subplot(subpno,1,2);
plot(sslicediff, 'x');
axis([0 imgno -Inf Inf]);
xlabel('Difference image number');
ylabel('Slice by slice variance');

subplot(subpno,1,3);
plot(globals/mom)
axis([0 imgno+1 -Inf Inf]);
xlabel('Image number');
ylabel('Scaled mean voxel intensity');

subplot(subpno,1,4);
mx = max(sslicediff);
mn = min(sslicediff);
avg = mean(sslicediff);
plot(avg, 'k');
axis([0 zno+1 -Inf Inf]);
hold on
plot(mn, 'b');
plot(mx, 'r');
hold off
xlabel('Slice number');
ylabel('Max/mean/min slice variance');

% realignment params
if any(flags == 'r')
  subplot(subpno,1,5);
  plot(mparams(:,1:3))
  legend('x translation','y translation','z translation',0);
  xlabel('image')
  ylabel('translations in mm')
end

% and label with first image at bottom
cp = get(gca,'Position');
wfa =  axes('Position', [0 0 1 1], 'Visible', 'off');
img1 = deblank(imgs(1,:));
text(0.5,cp(2)/2.5,{'First image:',img1},'HorizontalAlignment','center');

return
