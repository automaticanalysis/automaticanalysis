function data = ft_detectpeak(cfg,data)

% FT_DETECTPEAK calculates the amplitude and the latency of timelocked data
%
% You can specify the following configuration options:
%   cfg.parameter      = string indicating which parameter to analyse. default is set to 'avg', if it is present in the data
%   cfg.inflection     = string indicating whether we are looking for positive ('positive'), negative ('negative'), or any ('absolute') inflection (default = 'positive')
%   cfg.latency        = [begin end] in seconds, or 'all', 'prestim', 'poststim' (default = 'all')
%   cfg.neighbourwidth = width of neighbouring samples to compare the peak to, in seconds (default = 0.02)

% Copyright (C) 2020, Tibor Auer
%
% $Id$

% set the defaults
cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
cfg.inflection = ft_getopt(cfg, 'inflection', 'positive');
cfg.latency = ft_getopt(cfg, 'latency', 'all');
cfg.neighbourwidth = ft_getopt(cfg, 'neighbourwidth', 0.02);

if ischar(cfg.latency)
    switch cfg.latency
        case 'all'
            cfg.latency = [data.time(1) data.time(end)];
        case 'prestim'
            cfg.latency = [data.time(1) data.time(find((data.time<0),1,'last'))];
        case 'poststim'
            cfg.latency = [data.time(find((data.time>0),1,'first')) data.time(end)];
    end
end

% prepare data
time = data.time;
srate = round(1/mean(diff(time)));
soi = round(((cfg.latency-time(1))*srate)+1);
nb = cfg.neighbourwidth*srate;
dims = strsplit(data.dimord,'_');
dimTime = find(strcmp(dims,'time'));
[dat,perm] = shiftdata(data.(cfg.parameter),dimTime);

switch cfg.inflection
    case 'positive'
    case 'negative'
        dat = -dat;
    case 'absolute'
        dat = abs(dat);
end

% prepare output
cfg.previous = data.cfg;
data = keepfields(data,{'label','elec'});
data.cfg = cfg;
dims(dimTime) = [];
data.dimord = strjoin(dims(perm(2:end)),'_');

% calculate
for x = 1:size(dat,2)
    for y = 1:size(dat,3)
        if any(isnan(dat(:,x,y))), continue; end
        ts = dat(soi(1):soi(2),x,y);
        [junk,ord] = sort(ts(nb+1:end-nb),'descend'); ord = ord+nb;
        for p = ord'
            if (ts(p) > mean(ts(p-nb:p-1))) && (ts(p) > mean(ts(p+nb:p+nb))), break; end
        end
        amp(x,y) = ts(p);
        lat(x,y) = time(p+soi(1)-1);
    end
end

% save
data.amp = amp;
data.lat = lat;