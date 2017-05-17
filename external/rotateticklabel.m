function th=rotateticklabel(h,rot)
%ROTATETICKLABEL rotates tick labels
%   TH=ROTATETICKLABEL(H,ROT) is the calling form where H is a handle to
%   the axis that contains the XTickLabels that are to be rotated. ROT is
%   an optional parameter that specifies the angle of rotation. The default
%   angle is 90. TH is a handle to the text objects created. For long
%   strings such as those produced by datetick, you may have to adjust the
%   position of the axes so the labels don't get cut off.
%
%   Of course, GCA can be substituted for H if desired.
%
%   TH=ROTATETICKLABEL([],[],'demo') shows a demo figure.
%
%   Known deficiencies: if tick labels are raised to a power, the power
%   will be lost after rotation.
%
%   See also datetick.

%   Written Oct 14, 2005 by Andy Bliss
%   Modified Dec 6, 2011 by Alejandro Vicente Grabovetsky
%   Modified Nov 28, 2013 by Tibor Auer
%   Copyright 2005 by Andy Bliss

%set the default rotation if user doesn't specify
if nargin==1
    rot=90;
end
%make sure the rotation is in the range 0:360 (brute force method)
while rot>360
    rot=rot-360;
end
while rot<0
    rot=rot+360;
end
%get current tick labels
a=get(h,'XTickLabel');
%erase current tick labels from figure
set(h,'XTickLabel',[]);
%get tick label positions
b=get(h,'XTick');
c=get(h,'YTick');
%make new tick labels
switch get(gca, 'XAxisLocation')
    case 'top'
        th=text(b,repmat(c(1),length(b),1)*1.03 - 1/2,a,'HorizontalAlignment','left','rotation',rot);
    case 'bottom'
        th=text(b,repmat(c(end),length(b),1)*1.03 + 1/2,a,'HorizontalAlignment','right','rotation',rot);
end
