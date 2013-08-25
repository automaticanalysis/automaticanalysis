% Example user script for JP_SPM8_SURFACERENDER. Set the img variable
% to the image (e.g., spmT*) you want to render, with options set in 
% the CFG structure. For a full list of options, see the help for
% JP_SPM8_SURFACERENDER.

clear all
close all

cfg = [];
img = '/path/to/tstat.img';

jp_spm8_surfacerender(img, 'hot', cfg);


%% Try again but change some default values

cfg.colorscale = [3 6]; % min and max of color scale
cfg.inflate = 13;       % amount of inflation
cfg.plots = [1 2];      % which views to show ([1 2] is left and right)

jp_spm8_surfacerender(img, 'hot', cfg);

