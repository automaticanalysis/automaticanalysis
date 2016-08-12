% Modify the values in this script for your sequence
hdr.RepetitionTime = 2; % TR, units are s
hdr.EchoTime = 0.027;         % TE, units are s  
hdr.EchoSpacing = 0;       % needed if doing field map undistortion, units are s   
hdr.SliceTiming = 1:32;    % Time of each of the slices, in ms. Not needed if slice time correction is not included, or if slice order is specified manually
hdr.PhaseEncodingDirection = '';  % not essential if fieldmap undistortion is 
fname='demodicomhdr.json';
savejson('hdr',hdr,fname);
