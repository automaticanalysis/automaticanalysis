function cg_progress(action,varargin)
% Display progress and remaining time
% FORMAT cg_progress('Init',n_iterations,process_name,iteration_name)
% Initialises progress tool
%
% FORMAT cg_progress('Set',value)
% Sets iteration.
%
% FORMAT cg_progress('Clear')
% Clears the progress output and prints processing time.

global sum_time time_new time_old iter_old n_iterations arg2 arg3

if ~nargin, action = 'Init'; end

switch lower(action)
    % Initialise
    %-------------------------------------------------------------------
    case 'init'
        error(nargchk(2,5,nargin));
        n_iterations = varargin{1};
        if nargin > 2, arg2 = varargin{2}; else arg2 = 'Computing';  end
        if nargin > 3, arg3 = varargin{3}; else arg3 = 'Iterations'; end
        
        fprintf('%s\n',arg2);
        sum_time = 0;
        time_old = clock;
        iter_old = 0;

    % Set
    %-------------------------------------------------------------------
    case 'set'
        error(nargchk(2,3,nargin));
        iter = varargin{1};
        
        % estimate time for remaining iterations
        time_diff = etime(clock, time_old);
        sum_time = sum_time + time_diff;
        avg_time = sum_time/iter;
        str = sprintf('%.f%% (%s remaining)',iter/n_iterations*100,...
            time2str(avg_time*(n_iterations-iter)));
        fprintf('%-35s%-35s',repmat(sprintf('\b'),1,35),str);
        
        try

            h = axes('position',[0.5 0 0.1 0.05],'Units','normalized','Parent',...
                varargin{2},'Visible','off');

            text(0.5,0.5,sprintf('%-5s%-15s%-5s',repmat(sprintf(' '),5),str,repmat(sprintf(' '),5)),...
                'FontSize',spm('FontSize',8),...
                'HorizontalAlignment','Center',...
                'VerticalAlignment','middle',...
                'EraseMode','Background',...
                'BackgroundColor','white');

        end
  
        % save old values
        time_old = clock;
        iter_old = iter;
    
    % Clear
    %-------------------------------------------------------------------
    case 'clear'
        error(nargchk(1,2,nargin));
        fprintf('%-35s',repmat(sprintf('\b'),1,35));
        fprintf('Processing time for %d %s: %s\n',n_iterations,arg3,time2str(sum_time));

        try

            h = axes('position',[0.5 0 0.1 0.05],'Units','normalized','Parent',...
                varargin{1},'Visible','off');

            text(0.5,0.5,sprintf('Processing time for %d %s: %s\n',n_iterations,arg3,time2str(sum_time)),...
                'FontSize',spm('FontSize',8),...
                'HorizontalAlignment','Center',...
                'VerticalAlignment','middle',...
                'EraseMode','Background',...
                'BackgroundColor','white');

        end

    % Error
    %-------------------------------------------------------------------
    otherwise
        error('Unknown action string');
end

return

function str = time2str(t)

minutes = t/60;
hours   = t/3600;
days    = hours/24;

if days > 2
  str = sprintf('%d days %02.1f h', floor(days),24*(days-floor(days)));
elseif days > 1
  str = sprintf('%d day %02.1f h', floor(days),24*(days-floor(days)));
elseif hours > 1
  str = sprintf('%d:%02.0f h', floor(hours),60*(hours-floor(hours)));
elseif minutes > 1
  str = sprintf('%d:%02.0f min',floor(minutes),60*(minutes-floor(minutes)));
else
  str = sprintf('%d s',round(t));
end
