function varargout = run_tsdiffana(cmd, varargin)
% Wrapper function for tsdiffana routines

switch cmd,
    case 'run'
        subfun = varargin{1};
        job    = varargin{2};
        switch subfun
            case 'timediff'
                for k = 1:numel(job.imgs)
                    [p f e] = spm_fileparts(job.imgs{k}{1});
                    if job.vf
                        flags = 'mv';
                    else
                        flags = 'm';
                    end
                    imgs = char(job.imgs{k});
                    [td globals slicediff] = timediff(imgs,flags);
                    out.tdfn{k} = fullfile(p,'timediff.mat');
                    save(out.tdfn{k}, 'imgs', 'td', 'globals', 'slicediff');
                end
                varargout{1} = out;
            case 'tsdiffplot'
                fg = spm_figure('GetWin', 'Graphics');
                spm_figure('Clear');
                for k = 1:numel(job.tdfn)
                    h = tsdiffplot(job.tdfn{k}, fg);
                    spm_figure('NewPage', h);
                end
                if job.doprint
                    spm_figure('Print');
                end
        end
    case 'vout'
        subfun = varargin{1};
        job    = varargin{2};
        switch subfun
            case 'timediff'
                for k = 1:numel(job.imgs)
                    dep(k) = cfg_dep;
                    dep(k).sname      = sprintf('Timeseries Analysis Data File (%d)', k);
                    dep(k).src_output = substruct('.','tdfn','()',{k});
                    dep(k).tgt_spec   = cfg_findspec({{'filter','mat', ...
                        'strtype','e'}});
                end
                varargout{1}   = dep;
            case 'tsdiffplot'
        end
    case 'defaults'
        if nargin == 2
            varargout{1} = local_defs(varargin{1});
        else
            local_defs(varargin{1:2});
        end
end

function varargout = local_defs(defstr, defval)
persistent defs;
if isempty(defs)
    defs.vf = false;
    defs.doprint = true;
end
if nargin == 1
    varargout{1} = defs.(defstr);
else
    defs.(defstr) = defval;
end