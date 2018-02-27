function hs = tsdiffplot(tdfn,fgs,flags, varargin)
% tsfiffplot - plots image difference etc info
% FORMAT tsdiffplot(tdfn,fgs,flags, varargin)
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
    fgs = [];
end
if isempty(fgs)
    fgs(1) = spm_figure('Create', 'Graphics1'); spm_figure('Clear',fgs(1),'Graphics1');
    fgs(2) = spm_figure('Create', 'Graphics2'); spm_figure('Clear',fgs(2),'Graphics2');
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
imgno = numel(qa.global.mean);
zno =   size(qa.slice.mean,2);
mom = mean(qa.global.mean);
sslicediff = qa.slice.diff/mom;
slicemean_norm = qa.slice.mean - repmat(mean(qa.slice.mean,1),size(qa.slice.mean,1),1);

datatoplot = {...
    {@plot 2:imgno qa.global.diff/mom '-' 'Volume' 'Scaled variance'} {@imagesc 1:imgno 1:zno slicemean_norm' 'Volume' 'Scaled mean slice intensity'};...
    {@plot 2:imgno sslicediff '*' 'Volume' 'Slice by slice variance'} {@imagesc 2:imgno-1 1:zno log(qa.slice.fft)' 'Number of cycles in timecourse' 'FFT of slice intensity [log]'};...
    {@plot 1:imgno qa.global.mean/mom '-' 'Volume' 'Scaled mean voxel intensity'} {@plot 2:imgno-1 log(qa.global.fft) '-' 'Number of cycles in timecourse' 'FFT of mean intensity [log]'};...
    };

hs = [];
tickstep = round(imgno/100)*10;
dxt = tickstep:tickstep:imgno;

for m = 1:size(datatoplot,2)
    figure(fgs(m));
    
    for p = 1:size(datatoplot,1)
        h1 = axes('position', [.1 1-(.7+p-1)/subpno .6958 .65*1/subpno]);
        h2 = datatoplot{p,m}{1}(datatoplot{p,m}{2:4});
        axis([0.5 imgno+0.5 -Inf Inf]);
        set(h1,'xtick',dxt);
        xlabel(datatoplot{p,m}{5});
        ylabel(datatoplot{p,m}{6});
        if strcmp(func2str(datatoplot{p,m}{1}),'imagesc')
            colormap('jet')
            pos = get(h1,'position');
            colorbar;
            set(h1,'position',pos);
        end
        hs  = [hs; h2];
    end
    
    if m == 1
        axes('position', [.1 1-3.7/subpno .6958 .65*1/subpno]);
        mx = max(sslicediff);
        mn = min(sslicediff);
        avg = mean(sslicediff);
        h2 = plot(avg, 'k');
        axis([0 zno+1 -Inf Inf]);
        hold on
        h3 = plot(mn, 'b');
        h4 = plot(mx, 'r');
        hold off
        xlabel('Slice');
        ylabel('Slice variance');
        legend('Mean','Min','Max','Location','Best');
        hs  = [hs; h2; h3; h4];
    end
    
    % realignment params
    if any(flags == 'r')
        axes('position', [.1 1-4.7/subpno .6958 .65*1/subpno]);
        h2 = plot(mparams(:,1:3));
        legend('x translation','y translation','z translation',0);
        xlabel('image');
        ylabel('translations [mm]');
        hs  = [hs; h2];
    end
    
    % and label with first image at bottom
    cp = get(gca,'Position');
    axes('Position', [0 0 1 1], 'Visible', 'off');
    img1  = deblank(qa.imgs(1,:));
    text(0.5,cp(2)/2.5,{'First image:',img1},'HorizontalAlignment','center');
end
