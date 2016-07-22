% Modify the values in this script for your sequence
hdr.RepetitionTime = 2000; % TR, units are ms
hdr.EchoTime = 27;         % TE, units are ms  
hdr.EchoSpacing = 0;       % needed if doing field map undistortion, units are ms   
hdr.SliceTiming = 1:32;    % slice order. 1:numslices for ascending sequential; numslices:-1:1 for descending sequential. For interleaved, be careful - need [1:2:numslices 2:2:numslices] when odd number of slices otherwise [2:2:numslices 1:2:numslices] when even numslices
hdr.PhaseEncodingDirection = '';  % not essential
fname='demodicomhdr.json';
savejson('hdr',hdr,fname);
