function outData = ft_average_bands(cfg,ipfdata,inData)
% FT_AVERAGE_BANDS calculates average of frequency related metrics according to band specification
%
% You can specify the following configuration options:
%   cfg.parameter   = string or cell array indicating which parameter(s) to analyse. default is set to 'powspctrm', if it is present in the data
%   cfg.bandspec    = structure with band specification. it has two fields
%                       band, 1xN cell array with names of the bands as a subset of {'delta' 'theta' 'alpha' 'beta' 'low gamma' 'high gamma'}
%                           (default = {'alpha' 'beta'})
%                       bandbound, a 1xN cell array with lower and upper frequency bounds (as 1x2 numeric array) corresponding to cfg.bandspec.band 
%                           (default = {[8 13] [14 32]})
%
% Optionally, you can specify a band data indicating individual peaks and bandwidths of frequency bands to adjust band specification. 
%
%   An example of a band data structure containing the individual band information for 61 channels and 6 bands is
%
%               label: {61x1 cell}      the channel labels
%              dimord: 'chan_band'      defines how the numeric data should be interpreted
%                band: {1x6 cell}       names of the bands as a subset of {'delta' 'theta' 'alpha' 'beta' 'low gamma' 'high gamma'}
%            peakfreq: [61x6 double]    the peak of the power spectum within bands
%       peakbandwidth: [61x6 double]    the bandwidth of the peaks of the power spectum within bands
%             peakpow: [61x6 double]    the power of the peaks of the power spectum within bands
%                elec: [1×1 struct]     structure with electrode positions corresponding to 'label'
%            bandspec: [1×1 struct]     structure with band specification corresponding to 'band'. see FT_AVERAGE_BANDS
%                 cfg: [1x1 struct]     the configuration used by the function that generated this data structure
%
%   Required fields:
%       - label, dimord, band, peakfreq, peakbandwidth
%
%   Optional fields:
%       - elec, bandspec, peakpow

% Copyright (C) 2022, Tibor Auer
%
% $Id$

% set defaults and sanity check
ipfdataBands = {'alpha' 'beta'};
defaultbandspec.band = {'alpha' 'beta'};
defaultbandspec.bandbound = {[8 13.5] [13.5 32]};
defaultbands = {'delta' 'theta' 'alpha' 'beta' 'low gamma' 'high gamma'};

cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
if ~iscell(cfg.parameter), cfg.parameter = {cfg.parameter}; end

cfg.bandspec = ft_getopt(cfg, 'bandspec', defaultbandspec);
if ~isempty(setdiff(cfg.bandspec.band,defaultbands))
    ft_error('bandspecification MUST be a subset of {''delta'' ''theta'' ''alpha'' ''beta'' ''low gamma'' ''high gamma''}');
end
if numel(intersect(cfg.bandspec.band,defaultbandspec.band)) ~= 2
    ft_error('bandspecification MUST include {''alpha'' ''beta''}');
end

% adjust bands
bandspec = cfg.bandspec;

if ~isempty(ipfdata) % adjust bands according to ipfdata
    % prepare: ipfdata -> freqdata
    ipfdata.dimord = 'chan_freq';
    if isfield(ipfdata,'bandspec')
        ipfdata.freq = cellfun(@mean, ipfdata.bandspec.bandbound);
        ipfdata = rmfield(ipfdata,{'band' 'bandspec'});
    else
    end
    
    tmpcfg = [];
    tmpcfg.avgoverchan = 'yes';
    tmpcfg.frequency = cellfun(@(b) mean(get_bandbound(cfg.bandspec,b)), ipfdataBands);
    tmpcfg.avgoverfreq = 'no';
    tmpcfg.nanmean = 'yes';
    ipfdata = ft_selectdata(tmpcfg,ipfdata);
    % cfg.bandspec for invalid bandspec in ipfdata
    for b = find(isnan(ipfdata.peakfreq))
        cfgbandspec = get_bandbound(cfg.bandspec,ipfdataBands{b});
        ipfdata.peakfreq(b) = mean(cfgbandspec);
        ipfdata.peakbandwidth(b) = diff(cfgbandspec);
    end
    alphabetaDiff = diff(ipfdata.peakfreq);
    
    if alphabetaDiff > sum(ipfdata.peakbandwidth)/2 % enough room between peaks
        alphaBound = [...
            inData.freq(find(inData.freq>ipfdata.peakfreq(1)-ipfdata.peakbandwidth(1)/2,1,'first')) ... % first above the lower bound (round up)
            inData.freq(find(inData.freq>ipfdata.peakfreq(1)+ipfdata.peakbandwidth(1)/2,1,'first')) ... % first above the upper bound (round up)
            ];
        betaBounds = [...
            inData.freq(find(inData.freq>ipfdata.peakfreq(2)-ipfdata.peakbandwidth(2)/2,1,'first')) ... % first above the lower bound (round up)
            inData.freq(find(inData.freq>ipfdata.peakfreq(2)+ipfdata.peakbandwidth(2)/2,1,'first')) ... % first above the upper bound (round up)
            ];
    else % squashed bands
        halfBW = ipfdata.peakbandwidth*diff(ipfdata.peakfreq)/sum(ipfdata.peakbandwidth);
        
        alphaBound = [...
            inData.freq(find(inData.freq<ipfdata.peakfreq(1)-ipfdata.peakbandwidth(1)/2,1,'last')) ... % first below the lower bound (compensate for squash)
            inData.freq(find(inData.freq>ipfdata.peakfreq(1)+halfBW(1),1,'first')) ... % first above the (squashed) upper bound (round up)
            ];
        
        betaBounds = [...
            inData.freq(find(inData.freq>alphaBound(2),1,'first')) ... % first above the alpha upper bound
            inData.freq(find(inData.freq>ipfdata.peakfreq(2)+ipfdata.peakbandwidth(2)/2,1,'first')) ... % first above the upper bound
            ];
    end
    bandspec = set_bandbound(bandspec,'alpha',alphaBound);
    bandspec = set_bandbound(bandspec,'beta',betaBounds);
    
    if any(strcmp(bandspec.band,'theta'))
        thetaBound(2) = inData.freq(find(inData.freq<alphaBound(1),1,'last')); % first below the alpha lower bound
        thetaBound(1) = inData.freq(find(inData.freq<thetaBound(2)-diff(get_bandbound(cfg.bandspec,'theta')),1,'last')); % first below the [upper bound minus theta range]
        bandspec = set_bandbound(bandspec,'theta',thetaBound);
    end
    
    if any(strcmp(bandspec.band,'delta'))
        bandspec = set_bandbound(bandspec,'delta', [...
            inData.freq(1) ...
            inData.freq(find(inData.freq<thetaBound(1),1,'last')) ... % first below the theta lower bound
            ]);
    end
else % adjust bands only according to the freq data
    for band = {'delta' 'theta' 'alpha' 'beta'}
        if any(strcmp(bandspec.band,band{1}))
            bandBound = get_bandbound(cfg.bandspec,band{1});
            bandBound = [...
                inData.freq(find(inData.freq>bandBound(1),1,'first')) ... % first above the lower bound
                inData.freq(find(inData.freq<bandBound(2),1,'last')) ... % first below the upper bound
                ];
            bandspec = set_bandbound(bandspec,band{1},bandBound);
        end
    end
end

for band = {'low gamma' 'high gamma'} % adjust gamma bands only according to the freq data
    if any(strcmp(bandspec.band,band{1}))
        bandBound = get_bandbound(cfg.bandspec,band{1});
        bandBound = [...
            inData.freq(find(inData.freq>bandBound(1),1,'first')) ... % first above the lower bound
            inData.freq(find(inData.freq<bandBound(2),1,'last')) ... % first below the upper bound
            ];
        bandspec = set_bandbound(bandspec,band{1},bandBound);
    end
end

% average data within bands
outData = inData;
if isfield(inData,'avg'), dims = strsplit(inData.avg.dimord,'_');
else, dims = strsplit(inData.dimord,'_'); end
dimFreq = find(strcmp(dims,'freq'));

for p = cfg.parameter
    if isfield(inData,'avg')
        outData.avg.(p{1}) = [];
        if isfield(inData.avg,'dimord'), outData.avg.dimord = strrep(inData.avg.dimord,'freq','band'); end
    else
        outData.(p{1}) = []; 
        if isfield(inData,'dimord'), outData.dimord = strrep(inData.dimord,'freq','band'); end
    end    
end
outData.band = cfg.bandspec.band;
outData.bandspec = bandspec;

for b = 1:numel(outData.band)
    tmpcfg = [];
    tmpcfg.frequency = get_bandbound(outData.bandspec,outData.band{b});
    tmpcfg.avgoverfreq = 'yes';
    tmpcfg.nanmean = 'yes';
    dat = ft_selectdata(tmpcfg,inData);

    for p = cfg.parameter
        if isfield(inData,'avg'), outData.avg.(p{1}) = cat(dimFreq,outData.avg.(p{1}),dat.(p{1}));
        else, outData.(p{1}) = cat(dimFreq,outData.(p{1}),dat.(p{1})); end
    end
end

outData = rmfield(outData,'freq');
cfg.previous = inData.cfg;
outData.cfg = cfg;

end

function bound = get_bandbound(bandspec,band)
if ~any(strcmp(bandspec.band,band)), ft_error('band %s is not part of the band specification',band); end
bound = bandspec.bandbound{strcmp(bandspec.band,band)};
end

function bandspec = set_bandbound(bandspec,band,bound)
if ~any(strcmp(bandspec.band,band)), ft_error('band %s is not part of the band specification',band); end
bandspec.bandbound{strcmp(bandspec.band,band)} = bound;
end