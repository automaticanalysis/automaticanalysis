function [h1 hh1 obj] = showSurf(fnL,fnR,th,lims,surf,cm,subj)

if nargin==1
    if ~isempty(strfind(fnL,'_rh'))
        fnR = fnL;
        fnL = [];
    end
    
    if ~isempty(strfind(fnL,'_lh'))
        fnR = [];
    end   
end

if nargin<6
    cm = 'jet';
end

if isobject(gcf)
    obj.figno = get(gcf,'Number');
else
    obj.figno = gcf;
end
obj.newfig = 1;

go1=0;
go2=0;
if ~isempty(fnL)
    obj.input_lh = fnL;
    go1=1;
    dispOpt = -1;
    h = MRIread(fnL);
else
    obj.input_lh = [];
end
if ~isempty(fnR)
    obj.input_rh = fnR;
    go2=1;
    dispOpt = 1;
    h = MRIread(fnR);
else
    obj.input_rh = [];
end

if go1 && go2
%     dispOpt=2;
    dispOpt=4;
end

if nargin>2 && ~isempty(th);
    obj.overlaythresh = th;
else
    obj.overlaythresh=[0 0];
end
    
if nargin>3 && ~isempty(lims);
    obj.colorlims = lims;
else
    obj.colorlims = [-inf inf];
end

if nargin>4 && ~isempty(surf);
    surf = surf;
else
    surf='pi';
end

switch surf
    case 'white'
        shading = 'curv';
        shadingrange = [-2 2];
    case 'pi'
        shading = 'mixed';
        %shading = 'logcurv';
        shadingrange = [-2 2];
    case 'inflated'
        shading = 'logcurv';
        shadingrange = [-.75 .75];
    case 'pial'
        shading = 'curv';
        shadingrange = [-2 3];
    otherwise
        shading = 'curv';
        shadingrange = [-2 3];
end


obj.colomap = cm;  % Choose your color map:  See colmap.m for options
obj.direction = '+';  % sets direction of one sided threshold + = greater than -= less than
obj.reverse = 0;  % Option to reverse the image; i.e. m*-1
obj.round = 0;  % if = 1, rounds all values on the surface to nearest whole number.  Useful for masks

if nargin<7
    switch numel(h.vol)
        case 40962
            obj.fsaverage = 'fsaverage6';  %% Set which fsaverage to map to e.g. fsaverage, fsaverage3, fsaverage6
        case 163842
            obj.fsaverage = 'fsaverage';
        otherwise
            error('Number of vertices do not match fsaverage or fsaverage5');
    end
else
    obj.fsaverage = subj;
end

obj.surface = surf;          %% Set the surface: inflated, pial, or white
obj.shading = shading;          %% Set the shading information for the surface: curv, sulc, or thk
obj.shadingrange = shadingrange;    %% Set the min anx max greyscale values for the surface underlay (range of 0 to 1)

obj.Nsurfs = dispOpt; % %% Choose which hemispheres and surfaces to show:  4=L/R med/lat;  2= L/R lat; 1.9=L med/lat; 2.1 = R med/lat; -1= L lat; 1=R lat;

assignin('base','so',obj);

[h1 hh1] = surfPlot2(obj);

% save tmp.mat obj
% set(hh1,'FaceAlpha',.75);