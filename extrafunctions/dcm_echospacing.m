function ESP = dcm_echospacing(dcm)
if ischar(dcm)
    h = dicominfo(dcm);
else
    h = dcm;
end

% Bandwidth per pixel phase encode
hstr = dec2hex(h.Private_0019_1028);
pBWpe = '';
for i = size(hstr,1):-1:1
    pBWpe = horzcat(pBWpe,hstr(i,:));
end
pBWpe = hex2num(pBWpe);

% AcquisitionMatrixText - fisrt is PE samples
hstr = char(h.Private_0051_100b)';
AMT = str2double(hstr(1:find(hstr=='*')-1));

% Echo spacing in ms
ESP = 1/(pBWpe * AMT) * 1000;