function makedicomhdr(hdr,fname)
% Specify acquisition timings according to your sequence !!!
% fname='manualdicomhdr.json';
% The followings are essential
% hdr.RepetitionTime =            2;      % TR, units are s
% hdr.EchoTime =                  0.027;  % TE, units are s  
%
% SliceTiming is needed if doing slice time correction (slicetiming)
% It can be specified as the slice acquisiton times relative to the TR, or
% as slice order assuming equidistant acquisition:
%   1:NumberOfSlices for ascending sequential
%   NumberOfSlices:-1:1 for descending sequential
%   For interleaved, be careful
%       [1:2:NumberOfSlices 2:2:NumberOfSlices] when odd slices first
%       [2:2:NumberOfSlices 1:2:NumberOfSlices] when even slices first
% hdr.SliceTiming =               1:32;
%
% The followings are needed if doing fieldmap-based undistortion (realignunwarp)
% hdr.EchoSpacing =               0.00051;  % , units are s   
% hdr.PhaseEncodingDirection =    'y-';       % not essential

if max(hdr.SliceTiming) > hdr.RepetitionTime % sliceorder is specified
    hdr.SliceTiming = (hdr.SliceTiming-1)*hdr.RepetitionTime/numel(hdr.SliceTiming);
end
savejson([],hdr,fname);