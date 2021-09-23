function varargout = BrainNet(varargin)
%BrainNet Viewer, a graph-based brain network mapping tool, by Mingrui Xia
%-----------------------------------------------------------
%	Copyright(c) 2019
%	Beijing Normal University
%	Written by Mingrui Xia
%	Mail to Author:  <a href="mingruixia@gmail.com">Mingrui Xia</a>
%   Version 1.63;
%   Create Date 20110906;
%   Last edited 20190329
%-----------------------------------------------------------
%
% BrainNet MATLAB code for BrainNet.fig
%      BrainNet, by itself, creates a new BrainNet or raises the existing
%      singleton*.
%
%      H = BrainNet returns the handle to a new BrainNet or the handle to
%      the existing singleton*.
%
%      BrainNet('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BrainNet.M with the given input arguments.
%
%      BrainNet('Property','Value',...) creates a new BrainNet or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrainNet_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrainNet_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrainNet

% Last Modified by GUIDE v2.5 21-Jun-2019 10:08:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BrainNet_OpeningFcn, ...
    'gui_OutputFcn',  @BrainNet_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin
    if ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
    if isstruct(varargin{1})
        h = figure;
        set(h,'Position',[100,100,1280,720]);
        set(h,'PaperPositionMode','auto');
        InitVariables;
        global EC       
        if nargin > 1 && isstruct(varargin{2})
            for f = fieldnames(varargin{2})'
                for ff = fieldnames(varargin{2}.(f{1}))'
                    EC.(f{1}).(ff{1}) = varargin{2}.(f{1}).(ff{1});
                end               
            end
        end
        BrainNet_LoadFiles(varargin{1});
        a = FileView;
        varargout = {h};
        return
    end    
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before BrainNet is made visible.
function BrainNet_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BrainNet (see VARARGIN)


warning off
fprintf('Please cite:\nXia M, Wang J, He Y (2013) BrainNet Viewer: A Network Visualization Tool for Human Brain Connectomics. PLoS ONE 8: e68910.\nAn example:\n''The brain networks were visualized with the BrainNet Viewer (http://www.nitrc.org/projects/bnv/) (Xia et al., 2013)''.\n');

% Choose default command line output for BrainNet
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% if exist('BrainNet_Icon.png','file')==2
%     newIcon = javax.swing.ImageIcon('BrainNet_Icon.png');
%     figFrame = get(handles.NV_fig,'JavaFrame');
%     figFrame.setFigureIcon(newIcon);
% end
% movegui(handles.NV_fig,'center');
movegui(hObject,'center');

%opengl software;  %% Add by Mingrui, 20140815 enable to fix upside-down text

InitVariables;

% if exist('BrainNet_Background.tif','file')==2
%     imshow(imread('BrainNet_Background.tif'));
% end

% Modified according to Chris Rorden's suggestion to avoid using addtional image toolbox.
if exist('BrainNet_Background.tif','file')==2 && exist('imshow','file')
    imshow(imread('BrainNet_Background.tif'));
end

global OutputText %%% Added by Mingrui Xia, 20120413, give out aal information for selected vertex.
if exist('BrainNet_AAL_Label.mat','file')==2
    OutputText.AAL = load('BrainNet_AAL_Label.mat');
else
    OutputText.AAL = [];
end
if exist('BrainNet_Brodmann_Label.mat','file')==2 %%% Added by Mingrui Xia, 20120417, give out brodmann information for selected vertex.
    OutputText.Brodmann = load('BrainNet_Brodmann_Label.mat');
else
    OutputText.Brodmann = [];
end
dcm_obj = datacursormode(hObject);
set(dcm_obj,'UpdateFcn',@BrainNet_Output_txt)
% set(hObject,'toolbar','figure')
% set(gcf,'Renderer','OpenGL');
%set(handles.NV_fig,'position',get(0,'screensize'));

% UIWAIT makes BrainNet wait for user response (see UIRESUME)
% uiwait(handles.NV_fig);


% --- Outputs from this function are returned to the command line.
function varargout = BrainNet_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%set(handles.NV_axes,'Position',get(handles.NV_fig,'Position'))

% --------------------------------------------------------------------
function NV_m_file_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NV_m_exit_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a
if ~isempty(a)
    delete(a);
    a=[];
end
close(findobj('Tag','NV_fig'));

% --------------------------------------------------------------------
function NV_m_nm_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG

if (isfield(FLAG,'IsCalledByREST') && FLAG.IsCalledByREST==1) || ...
        (isfield(FLAG,'IsCalledByCMD') && FLAG.IsCalledByCMD==1)
    % YAN Chao-Gan 111023. For calling from outside. Edited Mingrui Xia,
    % 20121031 add call from commandline
    handles.NV_axes=gca(hObject);
end
axes(handles.NV_axes);
cla
global a
if FLAG.Loadfile==0
    if exist('BrainNet_background.tif','file')==2
        imshow(imread('BrainNet_background.tif'));
    end
else
    if ~isempty(a) && mean(ishandle(a))==1
        delete(a);
    end
    a=FileView;
end
FLAG.IsCalledByREST = 0;

function NetPrepare
global surf
global EC
switch EC.edg.draw
    case 1
        [surf.ncyl,surf.cylinder]=Nettrans(surf.net);
        % cylinder(:,1:3): node index and value
        % cylinder(:,4): size
        % cylinder(:,5): color
        % cylinder(:,6): opacity
        
    case 2
        [surf.ncyl,surf.cylinder]=Nettrans(surf.net,EC.edg.draw_threshold);
end
surf.cylinder(:,4) = surf.cylinder(:,3);
surf.cylinder(:,5) = surf.cylinder(:,3);
surf.cylinder(:,6) = surf.cylinder(:,3);
if EC.edg.size_abs == 1
    surf.cylinder(:,4) = abs(surf.cylinder(:,4));
end
switch EC.edg.size
    case 1
        surf.cylinder(:,4) = EC.edg.size_size;
    case 2
        switch EC.edg.size_value
            case 1
                if max(surf.cylinder(:,4)) ~=min(surf.cylinder(:,4))
                    k = 5.6 /(max(surf.cylinder(:,4))-min(surf.cylinder(:,4)));
                    b = 6 - k*max(surf.cylinder(:,4));
                else
                    k = 0;
                    b = 2;
                end
                surf.cylinder(:,4) = surf.cylinder(:,4) * k + b;
            case 2
                tmp = surf.cylinder(:,4);
                tmp(tmp<0.2) = 0.2;
                tmp(tmp>6) = 6;
                surf.cylinder(:,4) = tmp;
        end
    case 3
        if EC.edg.size_value == 1
            indu = surf.cylinder(:,4) >= EC.edg.size_threshold;
            tmp = surf.cylinder(:,4);
            tmp2 = tmp(indu);
            if max(tmp2) ~= min(tmp2)
                k = 5.6 /(max(tmp2)-min(tmp2));
                b = 6 - k*max(tmp2);
            else
                k = 0;
                b = 2;
            end
            tmp2 = tmp2 * k + b;
            tmp(indu) = tmp2;
            tmp(surf.cylinder(:,4) < EC.edg.size_threshold) = 0.4;
            surf.cylinder(:,4) = tmp;
        elseif EC.edg.size_value == 2
            tmp = surf.cylinder(:,4);
            tmp(tmp<EC.edg.size_threshold) = 0.4;
            tmp(tmp>6) = 6;
            surf.cylinder(:,4) = tmp;
        end
end
if EC.edg.color_abs == 1
    surf.cylinder(:,5) = abs(surf.cylinder(:,5));
end
switch EC.edg.color
    case 1
        surf.cylinder(:,5) = 1;
    case 2
        tmp = surf.cylinder(:,5);
        if EC.edg.color_map_type == 1
            lowend = min(tmp);
            highend = max(tmp);
        else
            lowend = EC.edg.color_map_low;
            highend = EC.edg.color_map_high;
        end
        
        if highend~=lowend
            k = 63 / (highend - lowend);
            b = 64 - k*highend;
        else
            k = 0;
            b = 64;
        end
        
        tmp = round(tmp * k + b);
        tmp(tmp < 1) = 1;
        tmp(tmp > 64) = 64;
        surf.cylinder(:,5) = tmp;
    case 3
        tmp = surf.cylinder(:,5);
        tmp2 = tmp;
        tmp(tmp2 <= EC.edg.color_threshold) = 64;
        tmp(tmp2 > EC.edg.color_threshold) = 1;
        surf.cylinder(:,5) = tmp;
    case 4
        n1 = surf.sphere(surf.cylinder(:,1),1:3);
        n2 = surf.sphere(surf.cylinder(:,2),1:3);
        length_cyl = sqrt(sum((n1 - n2) .* (n1 - n2), 2));
        tmp = length_cyl;
        tmp(length_cyl >= EC.edg.color_distance) = 1;
        tmp(length_cyl < EC.edg.color_distance) = 64;
        surf.cylinder(:,5) = tmp;
    case 5
        edg_mod = surf.sphere(surf.cylinder(:,1),4);
        tmp = surf.sphere(surf.cylinder(:,1),4) - surf.sphere(surf.cylinder(:,2),4);
        edg_mod(tmp~=0) = 21;
        
        % Edited by Mingrui, 20150310, fix a bug when module index were not
        % noncontinuous.
        tmp = unique(surf.sphere(:,4));
        edg_mod2 = ones(size(edg_mod))*21;
        for i = 1:length(tmp)
            edg_mod2(edg_mod == tmp(i)) = i;
        end
        
        surf.cylinder(:,5) = edg_mod2;
    case 6
        index = sub2ind(size(surf.net),surf.cylinder(:,1),surf.cylinder(:,2));
        tmp = EC.edg.color_custom_matrix(index);
        tmp2 = zeros(size(tmp));
        for i = 1:length(EC.edg.color_custom_index)
            tmp2(tmp == EC.edg.color_custom_index(i)) = i;
        end
        surf.cylinder(:,5) = tmp2;
        %                 if max(tmp) ~=min(tmp)
        %                     k = 63 /(max(tmp)-min(tmp));
        %                     b = 64 - k*max(tmp);
        %                 else
        %                     k = 0;
        %                     b = 64;
        %                 end
        %                 tmp = round(tmp * k + b);
        %                 tmp(tmp < 1) = 1;
        %                 tmp(tmp > 64) = 64;
        %                 surf.cylinder(:,5) = tmp;
end

if EC.edg.opacity_abs == 1
    surf.cylinder(:,6) = abs(surf.cylinder(:,6));
end
switch EC.edg.opacity
    case 1
        surf.cylinder(:,6) = EC.edg.opacity_same;
    case 2
        tmp = surf.cylinder(:,6);
        k = (EC.edg.opacity_max - EC.edg.opacity_min)/(max(tmp)-min(tmp));
        b = EC.edg.opacity_max - k*max(tmp);
        tmp = tmp * k + b;
        tmp(tmp < EC.edg.opacity_min) = EC.edg.opacity_min;
        tmp(tmp > EC.edg.opacity_max) = EC.edg.opacity_max;
        surf.cylinder(:,6) = tmp;
end

function NodePrepare
global surf
global EC
if EC.nod.color==2
    surf.sphere(:,6)=surf.sphere(:,4);
    % Modified by Mingrui, 20170303, use only linear mapping, add fixed
    % colormap range
    
    %     if min(surf.sphere(:,6))<0
    %         surf.sphere(:,6)=surf.sphere(:,6)-min(surf.sphere(:,6));
    %     end
    %     if min(surf.sphere(:,6))<1
    %         surf.sphere(:,6)=surf.sphere(:,6)+1;
    %     end
    %     while max(surf.sphere(:,6))/min(surf.sphere(:,6))>10
    %         surf.sphere(:,6)=log(surf.sphere(:,6));
    %         if min(surf.sphere(:,6))<1
    %             surf.sphere(:,6)=surf.sphere(:,6)+1;
    %         end
    %     end
    if EC.nod.color_map_type == 1
        lowend = min(surf.sphere(:,6));
        highend = max(surf.sphere(:,6));
    else
        lowend = EC.nod.color_map_low;
        highend = EC.nod.color_map_high;
    end
    %Editied by Mingrui, 20181123, fix a bug when the length for inputed
    %colormap is not 64
    colormap_length = length(EC.nod.CM);
    if highend~=lowend
        EC.nod.k = (colormap_length-1)/(highend-lowend);
        EC.nod.b = colormap_length-EC.nod.k*highend;
    else
        EC.nod.k=0;
        EC.nod.b=1;
    end
    
    surf.sphere(:,6)=surf.sphere(:,6)*EC.nod.k+EC.nod.b;
else
    surf.sphere(:,6)=surf.sphere(:,4);
end
surf.sphere(:,7) = surf.sphere(:,5); % Edited by Mingrui Xia, 20120702, fix edge length when using fix node size.
switch EC.nod.size
    case 1
        surf.sphere(:,7) = EC.nod.size_size;
    case {2,3}
        if EC.nod.size_value == 1
            if min(surf.sphere(:,7))<0
                surf.sphere(:,7)=surf.sphere(:,7)-min(surf.sphere(:,7));
            end
            if min(surf.sphere(:,7))<1
                surf.sphere(:,7)=surf.sphere(:,7)+1;
            end
            while max(surf.sphere(:,7))/min(surf.sphere(:,7))>10
                surf.sphere(:,7)=log(surf.sphere(:,7));
                if min(surf.sphere(:,7))<1
                    surf.sphere(:,7)=surf.sphere(:,7)+1;
                end
            end
            if max(surf.sphere(:,7))~=min(surf.sphere(:,7))
                EC.nod.k=5/(max(surf.sphere(:,7))-min(surf.sphere(:,7)));
                EC.nod.b=7-EC.nod.k*max(surf.sphere(:,7));
            else
                EC.nod.k=0;
                EC.nod.b=4;
            end
            surf.sphere(:,7)=surf.sphere(:,7)*EC.nod.k+EC.nod.b;
        end
        if EC.nod.size == 3
            tmp = surf.sphere(:,7);
            tmp(tmp < EC.nod.size_threshold) = 1;
            surf.sphere(:,7) = tmp;
        end
end
surf.sphere(:,7) = surf.sphere(:,7) * EC.nod.size_ratio;

function fv = ROIPrepare
global EC
global surf

% Added by Mingrui 20140925, draw cluster in statistical files
if ~strcmp(surf.test,'No')
    vol_tmp = surf.vol;
    vol_tmp(vol_tmp < EC.vol.threshold & vol_tmp > -EC.vol.threshold) = 0;
    [L, num] = bwlabeln(vol_tmp,EC.vol.rmm);
    vol_tmp2 = zeros(size(vol_tmp));
    n = 0;
    for x = 1:num
        theCurrentCluster = L == x;
        if length(find(theCurrentCluster)) >= EC.vol.clustersize
            n = n + 1;
            vol_tmp2(logical(theCurrentCluster)) = n;
        end
    end
else
    vol_tmp2 = surf.vol;
end

fv = cell(length(EC.vol.roi.draw),1);
for i = 1:length(fv)
    vol = vol_tmp2;
    vol(vol ~= EC.vol.roi.draw(i)) = 0;
    vol(vol == EC.vol.roi.draw(i)) = 1;
    if EC.vol.roi.smooth == 1
        %         vol = smooth3(vol,'gaussian');
        if EC.vol.roi.smooth_kernal == 1
            vol = smooth3(vol);
        else
            vol = smooth3(vol,'gaussian');
        end
        %         vol = smooth3(vol);
    end
    fv{i,1} = isosurface(vol);
    while size(fv{i,1}.vertices,1) == 0 %% Added by Mingrui Xia, expand ROI while it is too small
        if EC.vol.roi.smooth_kernal == 1
            vol = smooth3(vol);
        else
            vol = smooth3(vol,'gaussian');
        end% Modified by Mingrui Xia 20160921, using Gaussian kernel instead of the 'box' kernel
        vol(vol>0) = 1;
        fv{i,1} = isosurface(vol);
    end
    coord = fv{i,1}.vertices(:,1);
    fv{i,1}.vertices(:,1) = fv{i,1}.vertices(:,2);
    fv{i,1}.vertices(:,2) = coord;
    fv{i,1}.vertices = fv{i,1}.vertices';
    fv{i,1}.vertices(4,:) = 1;
    fv{i,1}.vertices = surf.hdr.mat * fv{i,1}.vertices;
    fv{i,1}.vertices(4,:) = [];
    fv{i,1}.vertices = fv{i,1}.vertices';
    clear coord
    %     if length(find(fv{i,1}.vertices(:,1) < 0)) > length(fv{i,1}.vertices(:,1))/2
    % Edited by Mingrui Xia, 20160324, modify the methods for judging if
    % the ROI belongs to left or right side.
    if (min(fv{i,1}.vertices(:,1)) + max(fv{i,1}.vertices(:,1))) < 0
        fv{i,1}.side = 1;
    else
        fv{i,1}.side = 2;
    end
end

function [ncyl,cylinder]=Nettrans(net,t)
% Modified by Mingrui, 20170309, draw two edges for bidirected edges
global EC
global surf
switch EC.edg.draw_abs
    case 0
        temp=net;
    case 1
        temp=abs(net);
end
temp = temp - diag(diag(temp));
if nargin == 2
    temp(temp < t) = 0;
end
if EC.edg.directed == 0
    index = find(triu(temp) ~= 0); % Eddited by Mingrui Xia, 20120702, fix a drawall bug.
    
    %     cylinder(:,7) = 1;
else % Add by Mingrui Xia, 20120621, draw directed network.
    index = find(temp ~= 0);
    %
    %     temp(temp < t) = 0;
    %     net_up = triu(temp);
    %     net_low = tril(temp);
    %     net_conj = net_low';
    %     net_conj(net_up == 0) = 0;
    %     net_conj(net_conj ~= 0) = (net_up(net_conj ~= 0) + net_conj(net_conj ~= 0)) / 2;
    %     ind = find(net_conj ~= 0);
    %     cylinder1 = zeros(length(ind),4);
    %     [cylinder1(:,1),cylinder1(:,2)]=ind2sub(size(net),ind);
    %     cylinder1(:,3) = net(ind);
    %     cylinder1(:,7) = 2;
    %
    %     net_conj = net_conj + net_conj';
    %     temp(net_conj ~= 0) = 0;
    %     ind = find(temp ~= 0);
    %     cylinder2 = zeros(length(ind),4);
    %     [cylinder2(:,2),cylinder2(:,1)]=ind2sub(size(net),ind);
    %     cylinder2(:,3) = net(ind);
    %     cylinder2(:,7) = 1;
    %
    %     cylinder = [cylinder1',cylinder2']';
    %     ncyl = size(cylinder,1);
end
ncyl = length(index);
cylinder = zeros(ncyl,6);
[cylinder(:,1),cylinder(:,2)]=ind2sub(size(net),index);
cylinder(:,3) = net(index);


% Add by Mingrui Xia, 20120109, draw inter hemisphere edges.
if EC.edg.interhemiedges == 1
    index = surf.sphere(cylinder(:,1)).*surf.sphere(cylinder(:,2))>0;
    cylinder(index,:) = [];
    ncyl=size(cylinder,1);
end

function [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf)
t=size(surf.tri,1);
v=size(surf.coord,2);
tmax=max(surf.tri,[],2);
tmin=min(surf.tri,[],2);
if min(tmin(t/2+1:t))-max(tmax(1:ceil(t/2)))==1
    cut=t/2;
    cuv=v/2;
else
    for i=1:t-1
        tmax(i+1)=max(tmax(i+1),tmax(i));
        tmin(t-i)=min(tmin(t-i),tmin(t-i+1));
    end
    cut=min([find((tmin(2:t)-tmax(1:t-1))==1, 1 ) t]); %%% Edited by Mingrui, 20120320, fix a bug.
    cuv=tmax(cut);
end
tl=1:cut;
tr=(cut+1):t;
vl=1:cuv;
vr=(cuv+1):v;
h=0.39;
r=max(surf.coord,[],2)-min(surf.coord,[],2);
w1=h/r(2)*r(1)*3/4;
h1=h/r(2)*r(1);

function w1=NoMeshPage(surf)
h=0.39;
r=max(surf.sphere,[],1)-min(surf.sphere,[],1)*2;
w1=h/r(2)*r(1)*3/4;

function Viewer=GenerateView8(h1,w1)
Viewer(1,:)=[0.055,0.62,0.2925,0.4,-90,0];
Viewer(2,:)=[0.3,0.58,0.4,0.39,0,90];
Viewer(3,:)=[0.6525,0.62,0.2925,0.4,90,0];
Viewer(4,:)=[0.055,0.29,0.2925,0.4,90,0];
Viewer(5,:)=[0.3,0.18,0.4,0.39,0,-90];
Viewer(6,:)=[0.6525,0.29,0.2925,0.4,-90,0];
Viewer(7,:)=[0.055,0.02,w1,h1,180,0];
Viewer(8,:)=[0.945-w1,0.02,w1,h1,0,0];

function Viewer=GenerateView6(w1)
% Viewer(1,:)=[0.05,0.58,0.34,0.4,-90,0];
% Viewer(2,:)=[0.4,0.58,0.2,0.39,0,90];
% Viewer(3,:)=[0.61,0.58,0.34,0.4,180,0];
% Viewer(4,:)=[0.05,0.15,0.34,0.4,90,0];
% Viewer(5,:)=[0.4,0.15,0.2,0.39,0,-90];
% Viewer(6,:)=[0.61,0.15,0.34,0.4,0,0];
Viewer(1,:)=[0.055,0.58,0.2925,0.4,-90,0];
Viewer(2,:)=[0.3,0.58,0.4,0.39,0,90];
Viewer(3,:)=[0.6525,0.58,0.2925,0.4,90,0];
Viewer(4,:)=[0.055,0.15,0.2925,0.4,90,0];
Viewer(5,:)=[0.3,0.18,0.4,0.39,0,-90];
Viewer(6,:)=[0.6525,0.15,0.2925,0.4,-90,0];

function Viewer=GenerateView4
% YAN Chao-Gan 111023 Added. For Medium View (4 views)
Viewer(1,:)=[0.045,0.45,0.4387,0.6,-90,0];
Viewer(2,:)=[0.52,0.45,0.4387,0.6,90,0];
Viewer(3,:)=[0.045,-0.01,0.4387,0.6,90,0];
Viewer(4,:)=[0.52,-0.01,0.4387,0.6,-90,0];

function Viewer=GenerateView4v
% Lijie Huang 130221 Added. For Medium View (4 views) and 2 Ventral Views
Viewer(1,:)=[0.045,0.6,0.4387,0.36,-90,0];
Viewer(2,:)=[0.52,0.6,0.4387,0.36,90,0];
Viewer(3,:)=[0.045,0.23,0.4387,0.36,90,0];
Viewer(4,:)=[0.52,0.23,0.4387,0.36,-90,0];
Viewer(5,:)=[0.045,0.03,0.4387,0.18,-90,-90];
Viewer(6,:)=[0.52,0.03,0.4387,0.18,90,-90];

function Viewer = GenerateView5
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
Viewer(1,:) = [0.055,0.5,0.2925,0.4,-90,0];
Viewer(2,:) = [0.3,0.29,0.4,0.39,0,90];
Viewer(3,:) = [0.6525,0.5,0.2925,0.4,90,0];
Viewer(4,:) = [0.055,0.145,0.2925,0.4,90,0];
Viewer(5,:) = [0.6525,0.145,0.2925,0.4,-90,0];

function Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv)
global EC
Surfmatrix.tri{1}=surf.tri(tl,:);
Surfmatrix.tri{2}=surf.tri;
Surfmatrix.tri{3}=surf.tri(tr,:)-cuv;
Surfmatrix.tri{4}=surf.tri(tl,:);
Surfmatrix.tri{5}=surf.tri;
Surfmatrix.tri{6}=surf.tri(tr,:)-cuv;
Surfmatrix.tri{7}=surf.tri;
Surfmatrix.tri{8}=surf.tri;
Surfmatrix.coord{1}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{2}=[surf.coord(1,:)',surf.coord(2,:)',surf.coord(3,:)'];
Surfmatrix.coord{3}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
Surfmatrix.coord{4}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{5}=[surf.coord(1,:)',surf.coord(2,:)',surf.coord(3,:)'];
Surfmatrix.coord{6}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
Surfmatrix.coord{7}=[surf.coord(1,:)',surf.coord(2,:)',surf.coord(3,:)'];
Surfmatrix.coord{8}=[surf.coord(1,:)',surf.coord(2,:)',surf.coord(3,:)'];
if isfield(surf,'T')
    Surfmatrix.T{1}=surf.T(vl);
    Surfmatrix.T{2}=surf.T;
    Surfmatrix.T{3}=surf.T(vr);
    Surfmatrix.T{4}=surf.T(vl);
    Surfmatrix.T{5}=surf.T;
    Surfmatrix.T{6}=surf.T(vr);
    Surfmatrix.T{7}=surf.T;
    Surfmatrix.T{8}=surf.T;
end
if isfield(surf,'boundary_value')
    Surfmatrix.boundary_value{1}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{2}=surf.boundary_value;
    Surfmatrix.boundary_value{3}=surf.boundary_value(vr);
    Surfmatrix.boundary_value{4}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{5}=surf.boundary_value;
    Surfmatrix.boundary_value{6}=surf.boundary_value(vr);
    Surfmatrix.boundary_value{7}=surf.boundary_value;
    Surfmatrix.boundary_value{8}=surf.boundary_value;
end

if EC.msh.color_type == 2
    Surfmatrix.color_table{1} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{2} = EC.msh.color_table;
    Surfmatrix.color_table{3} = EC.msh.color_table(vr,:);
    Surfmatrix.color_table{4} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{5} = EC.msh.color_table;
    Surfmatrix.color_table{6} = EC.msh.color_table(vr,:);
    Surfmatrix.color_table{7} = EC.msh.color_table;
    Surfmatrix.color_table{8} = EC.msh.color_table;
end

function Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv)
% YAN Chao-Gan 111023 Added. For Medium View (4 views)
global EC
Surfmatrix.tri{1}=surf.tri(tl,:);
Surfmatrix.tri{2}=surf.tri(tr,:)-cuv;
Surfmatrix.tri{3}=surf.tri(tl,:);
Surfmatrix.tri{4}=surf.tri(tr,:)-cuv;
Surfmatrix.coord{1}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{2}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
Surfmatrix.coord{3}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{4}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
if isfield(surf,'T')
    Surfmatrix.T{1}=surf.T(vl);
    Surfmatrix.T{2}=surf.T(vr);
    Surfmatrix.T{3}=surf.T(vl);
    Surfmatrix.T{4}=surf.T(vr);
end
if isfield(surf,'boundary_value')
    Surfmatrix.boundary_value{1}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{2}=surf.boundary_value(vr);
    Surfmatrix.boundary_value{3}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{4}=surf.boundary_value(vr);
end

if EC.msh.color_type == 2
    Surfmatrix.color_table{1} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{2} = EC.msh.color_table(vr,:);
    Surfmatrix.color_table{3} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{4} = EC.msh.color_table(vr,:);
end


function Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv)
% Lijie Huang 130221 Added. For Medium View (4 views) and 2 ventral views
global EC
Surfmatrix.tri{1}=surf.tri(tl,:);
Surfmatrix.tri{2}=surf.tri(tr,:)-cuv;
Surfmatrix.tri{3}=surf.tri(tl,:);
Surfmatrix.tri{4}=surf.tri(tr,:)-cuv;
Surfmatrix.tri{5}=surf.tri(tl,:);
Surfmatrix.tri{6}=surf.tri(tr,:)-cuv;
Surfmatrix.coord{1}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{2}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
Surfmatrix.coord{3}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{4}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
Surfmatrix.coord{5}=[surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{6}=[surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
if isfield(surf,'T')
    Surfmatrix.T{1}=surf.T(vl);
    Surfmatrix.T{2}=surf.T(vr);
    Surfmatrix.T{3}=surf.T(vl);
    Surfmatrix.T{4}=surf.T(vr);
    Surfmatrix.T{5}=surf.T(vl);
    Surfmatrix.T{6}=surf.T(vr);
end
if isfield(surf,'boundary_value')
    Surfmatrix.boundary_value{1}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{2}=surf.boundary_value(vr);
    Surfmatrix.boundary_value{3}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{4}=surf.boundary_value(vr);
    Surfmatrix.boundary_value{5}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{6}=surf.boundary_value(vr);
end

if EC.msh.color_type == 2
    Surfmatrix.color_table{1} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{2} = EC.msh.color_table(vr,:);
    Surfmatrix.color_table{3} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{4} = EC.msh.color_table(vr,:);
    Surfmatrix.color_table{5} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{6} = EC.msh.color_table(vr,:);
end

function Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global EC
Surfmatrix.tri{1} = surf.tri(tl,:);
Surfmatrix.tri{2} = surf.tri;
Surfmatrix.tri{3} = surf.tri(tr,:)-cuv;
Surfmatrix.tri{4} = surf.tri(tl,:);
Surfmatrix.tri{5} = surf.tri(tr,:)-cuv;
Surfmatrix.coord{1} = [surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{2} = [surf.coord(1,:)',surf.coord(2,:)',surf.coord(3,:)'];
Surfmatrix.coord{3} = [surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
Surfmatrix.coord{4} = [surf.coord(1,vl)',surf.coord(2,vl)',surf.coord(3,vl)'];
Surfmatrix.coord{5} = [surf.coord(1,vr)',surf.coord(2,vr)',surf.coord(3,vr)'];
if isfield(surf,'T')
    Surfmatrix.T{1} = surf.T(vl);
    Surfmatrix.T{2} = surf.T;
    Surfmatrix.T{3} = surf.T(vr);
    Surfmatrix.T{4} = surf.T(vl);
    Surfmatrix.T{5} = surf.T(vr);
end
if isfield(surf,'boundary_value')
    Surfmatrix.boundary_value{1}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{2}=surf.boundary_value;
    Surfmatrix.boundary_value{3}=surf.boundary_value(vr);
    Surfmatrix.boundary_value{4}=surf.boundary_value(vl);
    Surfmatrix.boundary_value{5}=surf.boundary_value(vr);
end

if EC.msh.color_type == 2
    Surfmatrix.color_table{1} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{2} = EC.msh.color_table;
    Surfmatrix.color_table{3} = EC.msh.color_table(vr,:);
    Surfmatrix.color_table{4} = EC.msh.color_table(vl,:);
    Surfmatrix.color_table{5} = EC.msh.color_table(vr,:);
end

function a=PlotMesh4(Viewer,Surfmatrix,alpha)
% Mingrui Xia 111026 Added. For Medium View (4 views)
global EC
global cam
a=zeros(1,4);
for i=1:4
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),'EdgeColor','none');
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']);eval(['shading ',EC.glb.shading,';']);axis off
    set(Brain,'FaceColor',EC.msh.color);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    
    if EC.msh.color_type == 2
        Brain.FaceColor = 'flat';
        Brain.FaceVertexCData = Surfmatrix.color_table{i};
    end
    
    if EC.msh.boundary == 4
        Brain.FaceColor = 'flat';
        if EC.msh.color_type == 1
            Brain.FaceVertexCData = repmat(EC.msh.color,length(Surfmatrix.boundary_value{i}),1);
        else
            Brain.FaceVertexCData = Surfmatrix.color_table{i};
        end
        Brain.FaceVertexCData(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
    end
    
    if alpha~=1
        eval(['lighting ',EC.glb.lighting,';']);axis tight; axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 2
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end

function a=PlotMesh4v(Viewer,Surfmatrix,alpha)
% Mingrui Xia 111026 Added. For Medium View (4 views)
global EC
global cam
a=zeros(1,6);
for i=1:6
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),'EdgeColor','none');
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']);eval(['shading ',EC.glb.shading,';']);axis off
    set(Brain,'FaceColor',EC.msh.color);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    
    if EC.msh.color_type == 2
        Brain.FaceColor = 'flat';
        Brain.FaceVertexCData = Surfmatrix.color_table{i};
    end
    
    if EC.msh.boundary == 4
        Brain.FaceColor = 'flat';
        if EC.msh.color_type == 1
            Brain.FaceVertexCData = repmat(EC.msh.color,length(Surfmatrix.boundary_value{i}),1);
        else
            Brain.FaceVertexCData = Surfmatrix.color_table{i};
        end
        Brain.FaceVertexCData(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
    end
    
    if alpha~=1
        eval(['lighting ',EC.glb.lighting,';']);axis tight; axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 2
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end

function a=PlotMesh8(Viewer,Surfmatrix,alpha)
global EC
global cam
a=zeros(1,8);
for i=1:8
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),'EdgeColor','none');
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    set(Brain,'FaceColor',EC.msh.color);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    
    if EC.msh.color_type == 2
        Brain.FaceColor = 'flat';
        Brain.FaceVertexCData = Surfmatrix.color_table{i};
    end
    
    if EC.msh.boundary == 4
        Brain.FaceColor = 'flat';
        if EC.msh.color_type == 1
            Brain.FaceVertexCData = repmat(EC.msh.color,length(Surfmatrix.boundary_value{i}),1);
        else
            Brain.FaceVertexCData = Surfmatrix.color_table{i};
        end
        Brain.FaceVertexCData(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
    end
    
    if alpha~=1
        eval(['lighting ',EC.glb.lighting,';']);axis tight; axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 3
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end

function a=PlotMesh6(Viewer,surf,alpha)
global EC
global cam
a=zeros(1,6);
for i=1:6
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),'EdgeColor','none');
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    set(Brain,'FaceColor',EC.msh.color);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    
    if EC.msh.color_type == 2
        Brain.FaceColor = 'flat';
        Brain.FaceVertexCData = EC.msh.color_table;
    end
    
    if EC.msh.boundary == 4
        Brain.FaceColor = 'flat';
        if EC.msh.color_type == 1
            Brain.FaceVertexCData = repmat(EC.msh.color,surf.vertex_number,1);
        else
            Brain.FaceVertexCData = EC.msh.color_table;
        end
        Brain.FaceVertexCData(surf.boundary_value == 1,:) = repmat(EC.msh.boundary_color,sum(surf.boundary_value),1);
    end
    
    
    if alpha~=1
        axis tight; eval(['lighting ',EC.glb.lighting,';']);axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
end

function a=PlotMesh1(surf,alpha)
global EC
global cam
% a = axes('position',[0.03,0.03,0.94,0.94]);
a = axes;
Brain=trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),'EdgeColor','none');
if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
    hold on
    Brain2 = trisurf(surf.tri,surf.coord2(1,:),surf.coord2(2,:),surf.coord2(3,:),'EdgeColor','none');
    hold off
end
whitebg(gcf,EC.bak.color);
set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
eval(['material ',EC.glb.material,';']);
% material([0.1 0.9 0.2 23,0.1]);
eval(['shading ',EC.glb.shading,';']);axis off
set(Brain,'FaceColor',EC.msh.color);
set(Brain,'FaceAlpha',EC.msh.alpha);

if EC.msh.color_type == 2
    Brain.FaceColor = 'flat';
    Brain.FaceVertexCData = EC.msh.color_table;
end

if EC.msh.boundary == 4
    Brain.FaceColor = 'flat';
    if EC.msh.color_type == 1
        Brain.FaceVertexCData = repmat(EC.msh.color,surf.vertex_number,1);
    else
        Brain.FaceVertexCData = EC.msh.color_table;
    end
    Brain.FaceVertexCData(surf.boundary_value == 1,:) = repmat(EC.msh.boundary_color,sum(surf.boundary_value),1);
end

if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
    set(Brain2,'FaceColor',EC.msh.color);
    set(Brain2,'FaceAlpha',EC.msh.alpha);
    if EC.msh.boundary == 4
        Brain2.FaceColor = 'flat';
        Brain2.FaceVertexCData = repmat(EC.msh.color,surf.vertex_number,1);
        Brain2.FaceVertexCData(surf.boundary_value == 1,:) = repmat(EC.msh.boundary_color,sum(surf.boundary_value),1);
    end
end
daspect([1 1 1])
% Edited  by Mingrui Xia, 20120806, add custom view for single brain.
switch EC.lot.view_direction
    case 1
        view(-90,0);
    case 2
        view(0,90);
    case 3
        view(180,0);
    case 4
        view(EC.lot.view_az,EC.lot.view_el);
end

if alpha~=1
    axis tight;  axis vis3d off;eval(['lighting ',EC.glb.lighting,';']);
    cam = camlight(EC.glb.lightdirection);
end

function a = PlotMesh5(Viewer,Surfmatrix,alpha)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global EC
global cam
a = zeros(1,5);
for i = 1:5
    a(i) = axes('position',Viewer(i,1:4));
    Brain = trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),'EdgeColor','none');
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    set(Brain,'FaceColor',EC.msh.color);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    if EC.msh.boundary == 4
        Brain.FaceColor = 'flat';
        Brain.FaceVertexCData = repmat(EC.msh.color,length(Surfmatrix.boundary_value{i}),1);
        Brain.FaceVertexCData(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
    end
    if alpha~=1
        eval(['lighting ',EC.glb.lighting,';']);axis tight; axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 3
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end

function a = PlotROI1(fv,a,alpha)
% Added by Mingrui, 20120529, DrawROI
global EC
global cam
hold on
% axis manual
for i = 1:length(fv)
    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
    set(roi,'FaceColor',EC.vol.roi.color(i,:));
    %     set(roi,'FaceAlpha',20);
end
hold off
if alpha ~= 1
    axis tight;
    axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
    cam = camlight(EC.glb.lightdirection);
end

function a = PlotROI8(fv,a,alpha)
% Added by Mingrui, 20120529, DrawROI
global EC
global cam
for j = 1:8
    axes(a(j));
    switch j
        case {1,4}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 1
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
        case {2,5,7,8}
            hold on
            for i = 1:length(fv)
                roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                set(roi,'FaceColor',EC.vol.roi.color(i,:));
            end
            hold off
        case{3,6}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 2
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
    end
    if alpha ~= 1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(j) = camlight(EC.glb.lightdirection);
    end
end

function a = PlotROI6(fv,a,alpha)
% Added by Mingrui, 20120529, DrawROI
global EC
global cam
for j = 1:6
    axes(a(j));
    hold on
    for i = 1:length(fv)
        roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
        set(roi,'FaceColor',EC.vol.roi.color(i,:));
    end
    hold off
    if alpha ~= 1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(j) = camlight(EC.glb.lightdirection);
    end
end

function a = PlotROI4(fv,a,alpha)
global EC
global cam
for j = 1:4
    axes(a(j));
    switch j
        case {1,3}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 1
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
        case{2,4}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 2
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
    end
    if alpha ~= 1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(j) = camlight(EC.glb.lightdirection);
    end
end

function a = PlotROI4v(fv,a,alpha)
global EC
global cam
for j = 1:6
    axes(a(j));
    switch j
        case {1,3, 5}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 1
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
        case{2,4, 6}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 2
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
    end
    if alpha ~= 1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(j) = camlight(EC.glb.lightdirection);
    end
end

function a = PlotROI5(fv,a,alpha)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global EC
global cam
for j = 1:5
    axes(a(j));
    switch j
        case {1,4}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 1
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
        case {2}
            hold on
            for i = 1:length(fv)
                roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                set(roi,'FaceColor',EC.vol.roi.color(i,:));
            end
            hold off
        case{3,5}
            hold on
            for i = 1:length(fv)
                if fv{i,1}.side == 2
                    roi =  trisurf(fv{i,1}.faces,fv{i,1}.vertices(:,1),fv{i,1}.vertices(:,2),fv{i,1}.vertices(:,3),'EdgeColor','none');
                    set(roi,'FaceColor',EC.vol.roi.color(i,:));
                end
            end
            hold off
    end
    if alpha ~= 1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(j) = camlight(EC.glb.lightdirection);
    end
end

function a=MapMesh1(surf,low,high,alpha)
global EC
global cam
a=axes('position',[0.07,0.1,0.8,0.8]);
Brain=trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),surf.T,'EdgeColor','none');
% Edited by Mingrui Xia, 20111116, add translucency show
set(Brain,'FaceAlpha',EC.msh.alpha);
colormap(EC.vol.CM);
caxis([low,high]);
if EC.msh.boundary ~= 1
    index = fix((Brain.FaceVertexCData-low)/(high-low)*length(EC.vol.CM))+1;
    RGB = squeeze(ind2rgb(index,EC.vol.CM));
    RGB(surf.boundary_value == 1,:) = repmat(EC.msh.boundary_color,sum(surf.boundary_value),1);
    Brain.FaceVertexCData = RGB;
end
switch EC.lot.view_direction
    case 1
        view(-90,0);
    case 2
        view(0,90);
    case 3
        view(180,0);
    case 4
        view(EC.lot.view_az,EC.lot.view_el);
end
daspect([1 1 1]);
whitebg(gcf,EC.bak.color);
set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
eval(['material ',EC.glb.material,';']);
eval(['shading ',EC.glb.shading,';']);axis off
if alpha~=1
    axis tight; eval(['lighting ',EC.glb.lighting,';']); axis vis3d off;
    cam= camlight(EC.glb.lightdirection);
end
cb=colorbar('location','East');
set(cb,'Position',[0.9 0.1 0.03 0.3]);

% Modified by Mingrui, 20150115, support matlab 2014b
tmp = version;
ind = find(tmp == '.');
if str2double(tmp(1:ind(2)-1))<8.4
    set(cb,'YAxisLocation','right');
    set(cb,'YTick',get(cb,'YLim'));
else
    cb.AxisLocation = 'out';
    cb.Ticks = cb.Limits;
end

function a=MapMesh8(Viewer,Surfmatrix,low,high,alpha)
%%% Edited by Mingrui Xia, 20111116, add translucency show
global EC
global cam
a=zeros(1,8);
for i=1:8
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),Surfmatrix.T{i},'EdgeColor','none');
    colormap(EC.vol.CM);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    caxis([low,high])
    if EC.msh.boundary ~= 1
        index = fix((Brain.FaceVertexCData-low)/(high-low)*length(EC.vol.CM))+1;
        RGB = squeeze(ind2rgb(index,EC.vol.CM));
        RGB(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
        Brain.FaceVertexCData = RGB;
    end
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    if alpha~=1
        eval(['lighting ',EC.glb.lighting,';']);axis tight; axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 3
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end
cb=colorbar('location','South');
set(cb,'Position',[0.35 0.085 0.3 0.03]);

% Modified by Mingrui, 20150115, support matlab 2014b
tmp = version;
ind = find(tmp == '.');
if str2double(tmp(1:ind(2)-1))<8.4
    set(cb,'XAxisLocation','bottom');
    set(cb,'XTick',get(cb,'XLim'));
else
    cb.AxisLocation = 'out';
    cb.Ticks = cb.Limits;
end

function a=MapMesh4(Viewer,Surfmatrix,low,high,alpha)
% YAN Chao-Gan 111023 Added. For Medium View (4 views)
%%% Edited by Mingrui Xia, 20111116, add translucency show
global EC
global cam
a=zeros(1,4);
for i=1:4
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),Surfmatrix.T{i},'EdgeColor','none');
    set(Brain,'FaceAlpha',EC.msh.alpha);
    colormap(EC.vol.CM);
    caxis([low,high])
    if EC.msh.boundary ~= 1
        index = fix((Brain.FaceVertexCData-low)/(high-low)*length(EC.vol.CM))+1;
        RGB = squeeze(ind2rgb(index,EC.vol.CM));
        RGB(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
        Brain.FaceVertexCData = RGB;
    end
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    if alpha~=1
        axis tight; eval(['lighting ',EC.glb.lighting,';']); axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 2
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end
cb=colorbar('location','South');
set(cb,'Position',[0.4 0.055 0.2 0.03]);
% Modified by Mingrui, 20150115, support matlab 2014b
tmp = version;
ind = find(tmp == '.');
if str2double(tmp(1:ind(2)-1))<8.4
    set(cb,'XAxisLocation','bottom');
    set(cb,'XTick',get(cb,'XLim'));
else
    cb.AxisLocation = 'out';
    cb.Ticks = cb.Limits;
end

function a=MapMesh4v(Viewer,Surfmatrix,low,high,alpha)
% Lijie Huang 130221 Added. For Medium View (4 views) and 2 ventral views
%%% Edited by Mingrui Xia, 20111116, add translucency show
global EC
global cam
a=zeros(1,6);
for i=1:6
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),Surfmatrix.T{i},'EdgeColor','none');
    set(Brain,'FaceAlpha',EC.msh.alpha);
    colormap(EC.vol.CM);
    caxis([low,high])
    if EC.msh.boundary ~= 1
        index = fix((Brain.FaceVertexCData-low)/(high-low)*length(EC.vol.CM))+1;
        RGB = squeeze(ind2rgb(index,EC.vol.CM));
        RGB(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
        Brain.FaceVertexCData = RGB;
    end
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    if alpha~=1
        axis tight; eval(['lighting ',EC.glb.lighting,';']); axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 2
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end
cb=colorbar('location','South');
set(cb,'Position',[0.4 0.6 0.2 0.03]);
% Modified by Mingrui, 20150115, support matlab 2014b
tmp = version;
ind = find(tmp == '.');
if str2double(tmp(1:ind(2)-1))<8.4
    set(cb,'XAxisLocation','bottom');
    set(cb,'XTick',get(cb,'XLim'));
else
    cb.AxisLocation = 'out';
    cb.Ticks = cb.Limits;
end

function a=MapMesh6(Viewer,surf,low,high,alpha)
%%% Edited by Mingrui Xia, 20111116, add translucency show
global EC
global cam
a=zeros(1,6);
for i=1:6
    a(i)=axes('position',Viewer(i,1:4));
    Brain=trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),surf.T,'EdgeColor','none');
    set(Brain,'FaceAlpha',EC.msh.alpha);
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    colormap(EC.vol.CM);
    caxis([low,high]);
    if EC.msh.boundary ~= 1
        index = fix((Brain.FaceVertexCData-low)/(high-low)*length(EC.vol.CM))+1;
        RGB = squeeze(ind2rgb(index,EC.vol.CM));
        RGB(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
        Brain.FaceVertexCData = RGB;
    end
    if alpha~=1
        axis tight; eval(['lighting ',EC.glb.lighting,';']); axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
end
cb=colorbar('location','South');
set(cb,'Position',[0.35 0.085 0.3 0.03]);
% Modified by Mingrui, 20150115, support matlab 2014b
tmp = version;
ind = find(tmp == '.');
if str2double(tmp(1:ind(2)-1))<8.4
    set(cb,'XAxisLocation','bottom');
    set(cb,'XTick',get(cb,'XLim'));
else
    cb.AxisLocation = 'out';
    cb.Ticks = cb.Limits;
end

function a = MapMesh5(Viewer,Surfmatrix,low,high,alpha)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global EC
global cam
a = zeros(1,5);
for i = 1:5
    a(i) = axes('position',Viewer(i,1:4));
    Brain = trisurf(Surfmatrix.tri{i},Surfmatrix.coord{i}(:,1),Surfmatrix.coord{i}(:,2),Surfmatrix.coord{i}(:,3),Surfmatrix.T{i},'EdgeColor','none');
    colormap(EC.vol.CM);
    set(Brain,'FaceAlpha',EC.msh.alpha);
    caxis([low,high])
    if EC.msh.boundary ~= 1
        index = fix((Brain.FaceVertexCData-low)/(high-low)*length(EC.vol.CM))+1;
        RGB = squeeze(ind2rgb(index,EC.vol.CM));
        RGB(Surfmatrix.boundary_value{i} == 1,:) = repmat(EC.msh.boundary_color,sum(Surfmatrix.boundary_value{i}),1);
        Brain.FaceVertexCData = RGB;
    end
    view(Viewer(i,5:6));
    daspect([1 1 1]);
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']); eval(['material ',EC.glb.material,';']); eval(['shading ',EC.glb.shading,';']);axis off
    if alpha~=1
        eval(['lighting ',EC.glb.lighting,';']);axis tight; axis vis3d off;
        cam(i) = camlight(EC.glb.lightdirection);
    end
    if EC.glb.lr == 1
        switch i
            case 1
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'L','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
            case 3
                text(0,max(Surfmatrix.coord{i}(:,2)),max(Surfmatrix.coord{i}(:,3)),'R','FontSize',20,'FontWeight','Bold','HorizontalAlignment','center');
        end
    end
end
cb = colorbar('location','South');
set(cb,'Position',[0.35 0.085 0.3 0.03]);

% Modified by Mingrui, 20150115, support matlab 2014b
tmp = version;
ind = find(tmp == '.');
if str2double(tmp(1:ind(2)-1))<8.4
    set(cb,'XAxisLocation','bottom');
    set(cb,'XTick',get(cb,'XLim'));
else
    cb.AxisLocation = 'out';
    cb.Ticks = cb.Limits;
end

function a=PlotNode6(Viewer,surf,a,n,EC)
global cam
global FLAG
if length(a)>1
    for i=1:6
        axes(a(i));
        hold on
        for j=1:surf.nsph
            switch EC.nod.draw
                case 1
                    DrawSphere(surf,j,i,6,EC);
                case 2
                    switch EC.nod.draw_threshold_type
                        case 1
                            if surf.sphere(j,5)>EC.nod.draw_threshold
                                DrawSphere(surf,j,i,6,EC);
                            end
                        case 2
                            if surf.sphere(j,4)>EC.nod.draw_threshold
                                DrawSphere(surf,j,i,6,EC);
                            end
                    end
                case 3
                    if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                        DrawSphere(surf,j,i,6,EC);
                    end
            end
        end
        
        if n>1
            axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
            cam(i) = camlight(EC.glb.lightdirection);
        end
        hold off
    end
    % Modified by Mingrui 20150603, Add colorbar for nodes
    % Modified by Mingrui 20170303, Add fixed color map
    % Modified by Mingrui 20170807, Add colorbar for nodes in modular option
    % Modified by Mingrui 20181115, Fix a bug in display node color with volume
    % mapping
    if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
        
        switch(EC.nod.color)
            case 2
                tmp = EC.nod.CM;
            case 3
                tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
        end
        colormap(tmp);
        cb=colorbar('location','South');
        if EC.nod.color_map_type == 1
            caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
        else
            caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
        end
        set(cb,'Position',[0.35 0.085 0.3 0.03]);
        tmp = version;
        ind = find(tmp == '.');
        if str2double(tmp(1:ind(2)-1))<8.4
            set(cb,'XAxisLocation','bottom');
            set(cb,'XTick',get(cb,'XLim'));
        else
            cb.AxisLocation = 'out';
            cb.Ticks = cb.Limits;
        end
    end
else
    a=zeros(1,6);
    for i=1:6
        a(i)=axes('position',Viewer(i,1:4));
        hold on
        view(Viewer(i,5:6));
        daspect([1 1 1]);
        whitebg(gcf,EC.bak.color);
        set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
        eval(['lighting ',EC.glb.lighting,';']);
        for j=1:surf.nsph
            switch EC.nod.draw
                case 1
                    DrawSphere(surf,j,i,6,EC);
                case 2
                    switch EC.nod.draw_threshold_type
                        case 1
                            if surf.sphere(j,5)>EC.nod.draw_threshold
                                DrawSphere(surf,j,i,6,EC);
                            end
                        case 2
                            if surf.sphere(j,4)>EC.nod.draw_threshold
                                DrawSphere(surf,j,i,6,EC);
                            end
                    end
                case 3
                    if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                        DrawSphere(surf,j,i,6,EC);
                    end
            end
        end
        if n>1
            axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
            cam(i) = camlight(EC.glb.lightdirection);
        end
        hold off
    end
    % Modified by Mingrui 20150603, Add colorbar for nodes
    % Modified by Mingrui 20170807, Add colorbar for nodes in modular option
    % Modified by Mingrui 20181115, Fix a bug in display node color with volume
    % mapping
    if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
        
        switch(EC.nod.color)
            case 2
                tmp = EC.nod.CM;
            case 3
                tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
        end
        colormap(tmp);
        cb=colorbar('location','South');
        if EC.nod.color_map_type == 1
            caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
        else
            caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
        end
        set(cb,'Position',[0.35 0.085 0.3 0.03]);
        tmp = version;
        ind = find(tmp == '.');
        if str2double(tmp(1:ind(2)-1))<8.4
            set(cb,'XAxisLocation','bottom');
            set(cb,'XTick',get(cb,'XLim'));
        else
            cb.AxisLocation = 'out';
            cb.Ticks = cb.Limits;
        end
    end
end

function a=PlotNode4(surf, a, flag, centerX, n, EC)
% Mingrui Xia 111026 Added. For Medium View (4 views)
global cam
global FLAG
for i=1:4
    axes(a(i));
    hold on
    switch i
        case {1,3}
            if flag==1
                for j=1:surf.nsph
                    if surf.sphere(j,1)<centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,4,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,4,EC);
                                end
                        end
                    end
                end
            else
                for j=1:surf.nsph
                    if surf.sphere(j,1)>centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,4,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,4,EC);
                                end
                        end
                    end
                end
            end
        case{2,4}
            if flag==2
                for j=1:surf.nsph
                    if surf.sphere(j,1)<centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,4,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                    case 3
                                        if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                end
                        end
                    end
                end
            else
                for j=1:surf.nsph
                    if surf.sphere(j,1)>centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,4,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,4,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,4,EC);
                                end
                        end
                    end
                end
            end
    end
    if n>1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(i) = camlight(EC.glb.lightdirection);
    end
    
    
    hold off
end
% Modified by Mingrui 20150603, Add colorbar for nodes
% Modified by Mingrui 20170807, Add colorbar for nodes in modular option
% Modified by Mingrui 20181115, Fix a bug in display node color with volume
% mapping
if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
    switch(EC.nod.color)
        case 2
            tmp = EC.nod.CM;
        case 3
            tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
    end
    colormap(tmp);
    cb=colorbar('location','South');
    if EC.nod.color_map_type == 1
        caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
    else
        caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
    end
    set(cb,'Position',[0.4 0.055 0.2 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function a=PlotNode4v(surf, a, flag, centerX, n, EC)
% Lijie Huang 130221 Added. For Medium View (4 views) and 2 ventral views
global cam
global FLAG
for i=1:6
    axes(a(i));
    hold on
    switch i
        case {1,3,5}
            if flag==1
                for j=1:surf.nsph
                    if surf.sphere(j,1)<centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,5,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,5,EC);
                                end
                        end
                    end
                end
            else
                for j=1:surf.nsph
                    if surf.sphere(j,1)>centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,5,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,5,EC);
                                end
                        end
                    end
                end
            end
        case{2,4,6}
            if flag==2
                for j=1:surf.nsph
                    if surf.sphere(j,1)<centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,5,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,5,EC);
                                end
                        end
                    end
                end
            else
                for j=1:surf.nsph
                    if surf.sphere(j,1)>centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,5,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,5,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,5,EC);
                                end
                        end
                    end
                end
            end
    end
    if n>1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(i) = camlight(EC.glb.lightdirection);
    end
    
    hold off
end
% Modified by Mingrui 20150603, Add colorbar for nodes
% Modified by Mingrui 20170807, Add colorbar for nodes in modular option
% Modified by Mingrui 20181115, Fix a bug in display node color with volume
% mapping
if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
    
    switch(EC.nod.color)
        case 2
            tmp = EC.nod.CM;
        case 3
            tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
    end
    colormap(tmp);
    cb=colorbar('location','South');
    if EC.nod.color_map_type == 1
        caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
    else
        caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
    end
    set(cb,'Position',[0.4 0.6 0.2 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function a=PlotNode1(surf,a,n,EC)
global cam
global FLAG
if ~isempty(a)
    axes(a);
    hold on
    for j=1:surf.nsph
        switch EC.nod.draw
            case 1
                DrawSphere(surf,j,2,1,EC);
            case 2
                switch EC.nod.draw_threshold_type
                    case 1
                        if surf.sphere(j,5)>EC.nod.draw_threshold
                            DrawSphere(surf,j,2,1,EC);
                        end
                    case 2
                        if surf.sphere(j,4)>EC.nod.draw_threshold
                            DrawSphere(surf,j,2,1,EC);
                        end
                end
            case 3
                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                    DrawSphere(surf,j,2,1,EC);
                end
        end
    end
    
    % Modified by Mingrui 20150603, Add colorbar for nodes
    % Modified by Mingrui 20170807, Add colorbar for nodes in modular option
    % Modified by Mingrui 20181115, Fix a bug in display node color with volume
    % mapping
    if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
        
        switch(EC.nod.color)
            case 2
                tmp = EC.nod.CM;
            case 3
                tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
        end
        colormap(tmp);
        cb=colorbar('location','eastoutside');
        if EC.nod.color_map_type == 1
            caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
        else
            caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
        end
        set(cb,'Position',[0.9 0.1 0.03 0.3]);
        %         set(gca, 'FontSize', 30);
        tmp = version;
        ind = find(tmp == '.');
        if str2double(tmp(1:ind(2)-1))<8.4
            set(cb,'YAxisLocation','right');
            set(cb,'YTick',get(cb,'YLim'));
        else
            cb.AxisLocation = 'out';
            cb.Ticks = cb.Limits;
        end
    end
    
    
    hold off
    if n>1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam = camlight(EC.glb.lightdirection);
    end
else
    a=axes;
    hold on
    daspect([1 1 1]);
    % Edited  by Mingrui Xia, 20120806, add custom view for single brain.
    switch EC.lot.view_direction
        case 1
            view(-90,0);
        case 2
            view(0,90);
        case 3
            view(180,0);
        case 4
            view(EC.lot.view_az,EC.lot.view_el);
    end
    
    whitebg(gcf,EC.bak.color);
    set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
    eval(['lighting ',EC.glb.lighting,';']);
    for j=1:surf.nsph
        switch EC.nod.draw
            case 1
                DrawSphere(surf,j,2,1,EC);
            case 2
                switch EC.nod.draw_threshold_type
                    case 1
                        if surf.sphere(j,5)>EC.nod.draw_threshold
                            DrawSphere(surf,j,2,1,EC);
                        end
                    case 2
                        if surf.sphere(j,4)>EC.nod.draw_threshold
                            DrawSphere(surf,j,2,1,EC);
                        end
                end
            case 3
                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                    DrawSphere(surf,j,2,1,EC);
                end
        end
    end
    if n>1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam = camlight(EC.glb.lightdirection);
    end
    
    % Modified by Mingrui 20150603, Add colorbar for nodes
    % Modified by Mingrui 20170807, Add colorbar for nodes in modular option
    % Modified by Mingrui 20181115, Fix a bug in display node color with volume
    % mapping
    if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
        
        switch(EC.nod.color)
            case 2
                tmp = EC.nod.CM;
            case 3
                tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
        end
        colormap(tmp);
        cb=colorbar('location','East');
        if EC.nod.color_map_type == 1
            caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
        else
            caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
        end
        set(cb,'Position',[0.9 0.1 0.03 0.3]);
        tmp = version;
        ind = find(tmp == '.');
        if str2double(tmp(1:ind(2)-1))<8.4
            set(cb,'YAxisLocation','right');
            set(cb,'YTick',get(cb,'YLim'));
        else
            cb.AxisLocation = 'out';
            cb.Ticks = cb.Limits;
        end
    end
    
    hold off
end

function PlotNode8(surf,a,flag,centerX,n,EC)
global cam
global FLAG
for i=1:8
    axes(a(i));
    hold on
    switch i
        case {1,4}
            if flag==1
                for j=1:surf.nsph
                    if surf.sphere(j,1)<centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,8,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,8,EC);
                                end
                        end
                    end
                end
            else
                for j=1:surf.nsph
                    if surf.sphere(j,1)>centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,8,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,8,EC);
                                end
                        end
                    end
                end
            end
        case{3,6}
            if flag==2
                for j=1:surf.nsph
                    if surf.sphere(j,1)<centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,8,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,8,EC);
                                end
                        end
                    end
                end
            else
                for j=1:surf.nsph
                    if surf.sphere(j,1)>centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,8,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4)>EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,8,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                                    DrawSphere(surf,j,i,8,EC);
                                end
                        end
                    end
                end
            end
        case {2,5,7,8}
            for j=1:surf.nsph
                switch EC.nod.draw
                    case 1
                        DrawSphere(surf,j,i,8,EC);
                    case 2
                        switch EC.nod.draw_threshold_type
                            case 1
                                if surf.sphere(j,5)>EC.nod.draw_threshold
                                    DrawSphere(surf,j,i,8,EC);
                                end
                            case 2
                                if surf.sphere(j,4)>EC.nod.draw_threshold
                                    DrawSphere(surf,j,i,8,EC);
                                end
                        end
                    case 3
                        if ~isempty(find(surf.cylinder(:,1:2)==j, 1))
                            DrawSphere(surf,j,i,8,EC);
                        end
                end
            end
    end
    if n>1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(i) = camlight(EC.glb.lightdirection);
    end
    hold off
end
% Modified by Mingrui 20150603, Add colorbar for nodes
% Modified by Mingrui 20170807, Add colorbar for nodes in modular option
% Modified by Mingrui 20181115, Fix a bug in display node color with volume
% mapping
if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
    
    switch(EC.nod.color)
        case 2
            tmp = EC.nod.CM;
        case 3
            tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
    end
    colormap(tmp);
    cb=colorbar('location','South');
    if EC.nod.color_map_type == 1
        caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
    else
        caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
    end
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function a = PlotNode5(surf,a,flag,centerX,n,EC)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global cam
global FLAG
for i = 1:5
    axes(a(i));
    hold on
    switch i
        case {1,4}
            if flag == 1
                for j = 1:surf.nsph
                    if surf.sphere(j,1) < centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,52,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2) == j, 1))
                                    DrawSphere(surf,j,i,52,EC);
                                end
                        end
                    end
                end
            else
                for j = 1:surf.nsph
                    if surf.sphere(j,1) > centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,52,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2) == j, 1))
                                    DrawSphere(surf,j,i,52,EC);
                                end
                        end
                    end
                end
            end
        case{3,5}
            if flag == 2
                for j = 1:surf.nsph
                    if surf.sphere(j,1) < centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,52,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2) == j, 1))
                                    DrawSphere(surf,j,i,52,EC);
                                end
                        end
                    end
                end
            else
                for j = 1:surf.nsph
                    if surf.sphere(j,1) > centerX
                        switch EC.nod.draw
                            case 1
                                DrawSphere(surf,j,i,52,EC);
                            case 2
                                switch EC.nod.draw_threshold_type
                                    case 1
                                        if surf.sphere(j,5) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                    case 2
                                        if surf.sphere(j,4) > EC.nod.draw_threshold
                                            DrawSphere(surf,j,i,52,EC);
                                        end
                                end
                            case 3
                                if ~isempty(find(surf.cylinder(:,1:2) == j, 1))
                                    DrawSphere(surf,j,i,52,EC);
                                end
                        end
                    end
                end
            end
        case {2}
            for j = 1:surf.nsph
                switch EC.nod.draw
                    case 1
                        DrawSphere(surf,j,i,52,EC);
                    case 2
                        switch EC.nod.draw_threshold_type
                            case 1
                                if surf.sphere(j,5) > EC.nod.draw_threshold
                                    DrawSphere(surf,j,i,52,EC);
                                end
                            case 2
                                if surf.sphere(j,4) > EC.nod.draw_threshold
                                    DrawSphere(surf,j,i,52,EC);
                                end
                        end
                    case 3
                        if ~isempty(find(surf.cylinder(:,1:2) == j, 1))
                            DrawSphere(surf,j,i,52,EC);
                        end
                end
            end
    end
    if n > 1
        axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
        cam(i) = camlight(EC.glb.lightdirection);
    end
    hold off
end
% Modified by Mingrui 20150603, Add colorbar for nodes
% Modified by Mingrui 20170807, Add colorbar for nodes in modular option
% Modified by Mingrui 20181115, Fix a bug in display node color with volume
% mapping
if (EC.nod.color == 2||EC.nod.color == 3)&&(FLAG.Loadfile<4||EC.edg.color~=2)&&(FLAG.Loadfile<8||EC.vol.type~=1)
    
    switch(EC.nod.color)
        case 2
            tmp = EC.nod.CM;
        case 3
            tmp = EC.nod.CM(1:length(EC.nod.ModularNumber),:);
    end
    colormap(tmp);
    cb = colorbar('location','South');
    if EC.nod.color_map_type == 1
        caxis([min(surf.sphere(:,4)),max(surf.sphere(:,4))]);
    else
        caxis([EC.nod.color_map_low,EC.nod.color_map_high]);
    end
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

% Modified by Mingrui, 20150318, fix a compatity bug with matlab 2014b.
function PlotLabel6(surf,i,j)
global EC
if EC.lbl_front == 0
    switch i % Edited by Mingrui Xia, 111028. Adjust label postion. Edited by Mingrui Xia, 20111113, label position plus radius times ratio.
        case {1,4,3,6}
            text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case 2
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
            %     case 3
            %         if surf.sphere(j,1)<=0
            %             text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            %         else
            %             text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            %         end
        case 5
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
            %     case 6
            %         if surf.sphere(j,1)<=0
            %             text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            %         else
            %             text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            %         end
    end
else
    switch i
        case {1,6}
            pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
            text(pos-surf.sphere(j,7)-2,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case {3,4}
            pos = max([max(surf.coord(1,:)),max(surf.sphere(:,1))]);
            text(pos+surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            
        case 2
            pos = max([max(surf.coord(3,:)),max(surf.sphere(:,3))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
        case 5
            pos = min([min(surf.coord(3,:)),min(surf.sphere(:,3))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,pos-surf.sphere(j,7)-2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,pos-surf.sphere(j,7)-2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
    end
end

function PlotLabel8(surf,i,j)
global EC
if EC.lbl_front == 0
    switch i
        case {1,3,4,6} %%% Edited by Mingrui Xia, 111028. Combine several same situation, adjust label postion.Edited by Mingrui Xia, 20111113, label position plus radius times ratio.
            text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            %     case {4,6}
            %         text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)*2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            %     case 3
            %         text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)*2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case 2
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
        case 7
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            else
                text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            end
        case 5
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
        case 8
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
    end
else
    switch i
        case {1,6}
            pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
            text(pos-surf.sphere(j,7)-2,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            
        case {3,4}
            pos = max([max(surf.coord(1,:)),max(surf.sphere(:,1))]);
            text(pos+surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            
        case 2
            pos = max([max(surf.coord(3,:)),max(surf.sphere(:,3))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
        case 7
            pos = max([max(surf.coord(2,:)),max(surf.sphere(:,2))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),pos+surf.sphere(j,7)+1,surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            else
                text(surf.sphere(j,1),pos+surf.sphere(j,7)+1,surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            end
        case 5
            pos = min([min(surf.coord(3,:)),min(surf.sphere(:,3))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,pos-surf.sphere(j,7)-2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)-surf.sphere(j,7)-2,pos-surf.sphere(j,7)-2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
        case 8
            pos = min([min(surf.coord(2,:)),min(surf.sphere(:,2))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),pos+surf.sphere(j,7)+1,surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),pos+surf.sphere(j,7)+1,surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
    end
end

function PlotLabel5(surf,i,j)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global EC
if EC.lbl_front == 0
    switch i
        case {1,3,4,5}
            text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3) + surf.sphere(j,7) * EC.nod.size_ratio + 2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            
        case 2
            if surf.sphere(j,1) <= 0
                text(surf.sphere(j,1),surf.sphere(j,2) + surf.sphere(j,7) * EC.nod.size_ratio + 2,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2) + surf.sphere(j,7) * EC.nod.size_ratio + 2,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
    end
else
    switch i
        case {1,5}
            pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
            text(pos-surf.sphere(j,7)-2,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case {3,4}
            pos = max([max(surf.coord(1,:)),max(surf.sphere(:,1))]);
            text(pos+surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case 2
            pos = max([max(surf.coord(3,:)),max(surf.sphere(:,3))]);
            if surf.sphere(j,1)<=0
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
            else
                text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
            end
    end
end

function PlotLabel4v(surf,i,j)
global EC
if EC.lbl_front == 0
    switch i
        case {1,2,3,4}
            text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case 5
            text(surf.sphere(j,1)-surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case 6
            text(surf.sphere(j,1)+surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
    end
else
    switch i
        case {1,4}
            pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
            text(pos-surf.sphere(j,7)-2,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case {2,3}
            pos = max([max(surf.coord(1,:)),max(surf.sphere(:,1))]);
            text(pos+surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            
        case 5
            pos = min([min(surf.coord(3,:)),min(surf.sphere(:,3))]);
            text(surf.sphere(j,1)-surf.sphere(j,7)-2,surf.sphere(j,2),pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case 6
            pos = min([min(surf.coord(3,:)),min(surf.sphere(:,3))]);
            text(surf.sphere(j,1)+surf.sphere(j,7)+1,surf.sphere(j,2),pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
    end
end

function PlotLabel1(surf,j)
%%% Added by Mingrui Xia, 111028. Plot label for single view.
%%% Edited by Mingrui Xia, 20111113, label position plus radius times ratio.
global EC
if EC.lbl_front == 0
    if EC.msh.doublebrain == 1 %% Edited by Mingrui Xia, 20130818, fix the bug that the label only appear in the first brain
        switch EC.lot.view_direction
            case {1,4}
                %         text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
                text(0,surf.sphere2(j,2),surf.sphere2(j,3)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            case 2
                if surf.sphere2(j,1)<=0
                    text(surf.sphere2(j,1),surf.sphere2(j,2)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.sphere2(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                else
                    text(surf.sphere2(j,1),surf.sphere2(j,2)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.sphere2(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                end
            case 3
                if surf.sphere2(j,1)<=0
                    text(surf.sphere2(j,1),surf.sphere2(j,2),surf.sphere2(j,3)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                else
                    text(surf.sphere2(j,1),surf.sphere2(j,2),surf.sphere2(j,3)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                end
                
        end
    else
        switch EC.lot.view_direction
            case {1,4}
                %         text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
                text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            case 2
                if surf.sphere(j,1)<=0
                    text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                else
                    text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,surf.sphere(j,3),surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                end
            case 3
                if surf.sphere(j,1)<=0
                    text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                else
                    text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                end
                
        end
    end
else
    if EC.msh.doublebrain == 1
        switch EC.lot.view_direction
            case 1
                pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
                text(pos-surf.sphere2(j,7)*EC.nod.size_ratio-2,surf.sphere2(j,2),surf.sphere2(j,3)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            case 2
                pos = max([max(surf.coord(3,:)),max(surf.sphere(:,3))]);
                if surf.sphere2(j,1)<=0
                    text(surf.sphere2(j,1),surf.sphere2(j,2)+surf.sphere2(j,7)*EC.nod.size_ratio+1,pos+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                else
                    text(surf.sphere2(j,1),surf.sphere2(j,2)+surf.sphere2(j,7)*EC.nod.size_ratio+1,pos+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                end
            case 3
                pos = min([min(surf.coord(3,:)),min(surf.sphere(:,3))]);
                if surf.sphere2(j,1)<=0
                    text(surf.sphere2(j,1),surf.sphere2(j,2)-surf.sphere2(j,7)*EC.nod.size_ratio-2,pos-surf.sphere2(j,7)*EC.nod.size_ratio-2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                else
                    text(surf.sphere2(j,1),surf.sphere2(j,2)-surf.sphere2(j,7)*EC.nod.size_ratio-2,pos-surf.sphere2(j,7)*EC.nod.size_ratio-2,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                end
            case 4
                text(surf.sphere2(j,1),surf.sphere2(j,2),surf.sphere2(j,3)+surf.sphere2(j,7)*EC.nod.size_ratio+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
                
        end
    else
        switch EC.lot.view_direction
            case 1
                pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
                text(pos-surf.sphere(j,7)-2,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
            case 2
                pos = max([max(surf.coord(3,:)),max(surf.sphere(:,3))]);
                if surf.sphere(j,1)<=0
                    text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                else
                    text(surf.sphere(j,1),surf.sphere(j,2)+surf.sphere(j,7)+1,pos+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                end
            case 3
                pos = max([max(surf.coord(2,:)),max(surf.sphere(:,2))]);
                if surf.sphere(j,1)<=0
                    text(surf.sphere(j,1),pos+surf.sphere(j,7)+1,surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','left')
                else
                    text(surf.sphere(j,1),pos+surf.sphere(j,7)+1,surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','right')
                end
            case 4
                text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        end
    end
end

function PlotLabel4(surf,i,j)
% Mingrui Xia 111026 Added. For Medium View (4 views)
%%% Edited by Mingrui Xia, 20111113, label position plus radius times ratio.
global EC
if EC.lbl_front == 0
    text(surf.sphere(j,1),surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+1,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
else
    switch i
        case {1,4}
            pos = min([min(surf.coord(1,:)),min(surf.sphere(:,1))]);
            text(pos-surf.sphere(j,7)-2,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+3,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
        case {2,3}
            pos = max([max(surf.coord(1,:)),max(surf.sphere(:,1))]);
            text(pos+surf.sphere(j,7)+1,surf.sphere(j,2),surf.sphere(j,3)+surf.sphere(j,7)+3,surf.label{j},'FontName',EC.lbl_font.FontName,'FontWeight',EC.lbl_font.FontWeight,'FontAngle',EC.lbl_font.FontAngle,'FontSize',EC.lbl_font.FontSize,'FontUnits',EC.lbl_font.FontUnits,'HorizontalAlignment','center');
    end
end

function DrawSphere(surf,j,i,n,EC)
% switch EC.nod.size
%     case 1
% %         t=EC.nod.size_size;
% t = surf.sphere(j,7);
%     case 2
%         switch EC.nod.size_value
%             case 1
%                 t=surf.sphere(j,7);
%             case 2
%                 if surf.sphere(j,5)<1
%                     t=1;
%                 elseif surf.sphere(j,5)>10
%                     t=10;
%                 else
%                     t=surf.sphere(j,5);
%                 end
%         end
%     case 3
%         if surf.sphere(j,5)>EC.nod.size_threshold
%             switch EC.nod.size_value
%                 case 1
%                     t=surf.sphere(j,7);
%                 case 2
%                     if surf.sphere(j,5)>10
%                         t=10;
%                     else
%                         t=surf.sphere(j,5);
%                     end
%             end
%         else
%             t=1;
%         end
% end
% t=t*EC.nod.size_ratio;
t = surf.sphere(j,7); % Edited by Mingrui Xia, 20120702, fix edge length when using fix node size.
switch EC.glb.detail  % Add by Mingrui Xia, 20120413, adjust graph detail
    case 1
        [x,y,z]=sphere(10);
    case 2
        [x,y,z]=sphere(50);
    case 3
        [x,y,z]=sphere(100);
end
if EC.msh.doublebrain == 1
    x=x.*t+surf.sphere2(j,1);
    y=y.*t+surf.sphere2(j,2);
    z=z.*t+surf.sphere2(j,3);
else
    x=x.*t+surf.sphere(j,1);
    y=y.*t+surf.sphere(j,2);
    z=z.*t+surf.sphere(j,3);
end
Node=mesh(x,y,z,'EdgeColor','none');
switch EC.nod.color
    case 1
        ci=1;
    case 2
        ci=int32(surf.sphere(j,6));
        
        %Editied by Mingrui, 20181123, fix a bug when the length for
        %inputed colormap is not 64
        colormap_length = length(EC.nod.CM);
        if ci<1
            ci=1;
        elseif ci>colormap_length
            ci = colormap_length;
        end
        
    case 3
        ci = find(EC.nod.ModularNumber == surf.sphere(j,4));
        %         ci=int32(surf.sphere(j,4));
        %         if ci<1
        %             ci=22;
        %         elseif ci>21
        %             ci=22;
        %         end
    case 4
        if surf.sphere(j,4)>EC.nod.color_threshold
            ci=1;
        else
            ci=64;
        end
end
set(Node,'FaceColor',EC.nod.CM(ci,:));
set(Node,'EdgeAlpha',0)
% if surf.label{j}~='-'
if ~isequal(surf.label{j},'-') % Editied by Mingrui Xia, 20131226, fix the bug that the label doesn't displayed when including '-'
    switch EC.lbl
        case 1
            if n==8
                PlotLabel8(surf,i,j);
            elseif n==6
                PlotLabel6(surf,i,j);
            elseif n==4 %%% Edited by Mingrui Xia, 111028. Add label plot for medium view and single view.
                PlotLabel4(surf,i,j);
            elseif n == 5
                PlotLabel4v(surf,i,j);
            elseif n == 52 % Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
                PlotLabel5(surf,i,j);
            else
                PlotLabel1(surf,j);
            end
        case 3
            switch EC.lbl_threshold_type
                case 1
                    if surf.sphere(j,5)>EC.lbl_threshold
                        if n==8
                            PlotLabel8(surf,i,j);
                        elseif n==6
                            PlotLabel6(surf,i,j);
                        elseif n==4 %%% Edited by Mingrui Xia, 111028. Add label plot for medium view and single view.
                            PlotLabel4(surf,i,j);
                        elseif n == 5
                            PlotLabel4v(surf,i,j);
                        elseif n == 52 % Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
                            PlotLabel5(surf,i,j);
                        else
                            PlotLabel1(surf,j);
                        end
                    end
                case 2
                    if surf.sphere(j,4)>EC.lbl_threshold
                        if n==8
                            PlotLabel8(surf,i,j);
                        elseif n==6
                            PlotLabel6(surf,i,j);
                        elseif n==4 %%% Edited by Mingrui Xia, 111028. Add label plot for medium view and single view.
                            PlotLabel4(surf,i,j);
                        elseif n == 5
                            PlotLabel4v(surf,i,j);
                        elseif n == 52 % Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
                            PlotLabel5(surf,i,j);
                        else
                            PlotLabel1(surf,j);
                        end
                    end
            end
    end
end

function [centerX,flag]=JudgeNode(surf,vl,vr)
if abs(min(surf.coord(1,vl))-max(surf.coord(1,vr)))<abs(max(surf.coord(1,vl))-min(surf.coord(1,vr)))
    centerX=(min(surf.coord(1,vl))+max(surf.coord(1,vr)))/2;
else
    centerX=(max(surf.coord(1,vl))+min(surf.coord(1,vr)))/2;
end
if  max(surf.coord(1,vl))<max(surf.coord(1,vr))
    flag=1;
else
    flag=2;
end

function PlotLine6(surf,a)
global cam
global EC
for i=1:6
    axes(a(i));
    hold on
    for j=1:surf.ncyl
        DrawLine(surf,j,6,mod(i,3));
    end
    axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
    cam(i) = camlight(EC.glb.lightdirection);
    hold off
end
% Add by Mingrui, 20150114, show colorbar when using colormap
if EC.edg.color == 2
    colormap(EC.edg.CM);
    cb=colorbar('location','South');
    if EC.edg.color_map_type == 1 % Modified by Mingrui, 20170303, support fixed color mapping
        if EC.edg.color_abs == 0
            caxis([min(surf.cylinder(:,3)),max(surf.cylinder(:,3))]);
        else
            caxis([min(abs(surf.cylinder(:,3))),max(abs(surf.cylinder(:,3)))]);
        end
    else
        caxis([EC.edg.color_map_low,EC.edg.color_map_high]);
    end
    
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    % Modified by Mingrui, 20150115, support matlab 2014b
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function PlotLine1(surf)
global cam
global EC
hold on
switch EC.lot.view_direction
    case {1,4}
        for j=1:surf.ncyl
            DrawLine(surf,j,1,1);
        end
    case 2
        for j=1:surf.ncyl
            DrawLine(surf,j,1,2);
        end
    case 3
        for j=1:surf.ncyl
            DrawLine(surf,j,1,3);
        end
end
axis tight; axis vis3d off;daspect([1 1 1]);
eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
cam = camlight(EC.glb.lightdirection);

% Add by Mingrui, 20150114, show colorbar when using colormap
if EC.edg.color == 2
    hold on
    colormap(EC.edg.CM);
    cb=colorbar('location','East');
    if EC.edg.color_map_type == 1 % Modified by Mingrui, 20170303, support fixed color mapping
        if EC.edg.color_abs == 0
            caxis([min(surf.cylinder(:,3)),max(surf.cylinder(:,3))]);
        else
            caxis([min(abs(surf.cylinder(:,3))),max(abs(surf.cylinder(:,3)))]);
        end
    else
        caxis([EC.edg.color_map_low,EC.edg.color_map_high]);
    end
    set(cb,'Position',[0.9 0.1 0.03 0.3]);
    % Modified by Mingrui, 20150115, support matlab 2014b
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'YAxisLocation','right');
        set(cb,'YTick',get(cb,'YLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end
hold off

function PlotLine8(surf,a,flag,centerX)
global cam
global EC
for i=1:8
    axes(a(i));
    hold on
    switch i
        case {1,4}
            if flag==1
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,5,1);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,5,1);
                    end
                end
            end
        case{3,6}
            if flag==2
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,5,1);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,5,1);
                    end
                end
            end
        case {2,5}
            for j=1:surf.ncyl
                DrawLine(surf,j,5,2);
            end
        case {7,8}
            for j=1:surf.ncyl
                DrawLine(surf,j,5,3);
            end
    end
    axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
    cam(i) = camlight(EC.glb.lightdirection);
    hold off
end

% Add by Mingrui, 20150114, show colorbar when using colormap
if EC.edg.color == 2
    colormap(EC.edg.CM);
    cb=colorbar('location','South');
    if EC.edg.color_map_type == 1 % Modified by Mingrui, 20170303, support fixed color mapping
        if EC.edg.color_abs == 0
            caxis([min(surf.cylinder(:,3)),max(surf.cylinder(:,3))]);
        else
            caxis([min(abs(surf.cylinder(:,3))),max(abs(surf.cylinder(:,3)))]);
        end
    else
        caxis([EC.edg.color_map_low,EC.edg.color_map_high]);
    end
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
    
end

function PlotLine5(surf,a,flag,centerX)
% Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
global cam
global EC
for i = 1:5
    axes(a(i));
    hold on
    switch i
        case {1,4}
            if flag == 1
                for j = 1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1) < centerX && surf.sphere(surf.cylinder(j,2),1) < centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            else
                for j = 1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1) > centerX && surf.sphere(surf.cylinder(j,2),1) > centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            end
        case{3,5}
            if flag == 2
                for j = 1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1) < centerX && surf.sphere(surf.cylinder(j,2),1) < centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            else
                for j = 1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1) > centerX && surf.sphere(surf.cylinder(j,2),1) > centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            end
        case {2}
            for j = 1:surf.ncyl
                DrawLine(surf,j,3,2);
            end
    end
    axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
    cam(i) = camlight(EC.glb.lightdirection);
    hold off
end

% Add by Mingrui, 20150114, show colorbar when using colormap
if EC.edg.color == 2
    colormap(EC.edg.CM);
    cb = colorbar('location','South');
    if EC.edg.color_map_type == 1 % Modified by Mingrui, 20170303, support fixed color mapping
        if EC.edg.color_abs == 0
            caxis([min(surf.cylinder(:,3)),max(surf.cylinder(:,3))]);
        else
            caxis([min(abs(surf.cylinder(:,3))),max(abs(surf.cylinder(:,3)))]);
        end
    else
        caxis([EC.edg.color_map_low,EC.edg.color_map_high]);
    end
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function PlotLine4(surf,a,flag,centerX)
global cam
global EC
for i=1:4
    axes(a(i));
    hold on
    switch i
        case {1,3}
            if flag==1
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,2,1);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,2,1);
                    end
                end
            end
        case{2,4}
            if flag==2
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,2,1);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,2,1);
                    end
                end
            end
    end
    axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
    cam(i) = camlight(EC.glb.lightdirection);
    
    
    hold off
end

% Add by Mingrui, 20150114, show colorbar when using colormap
if EC.edg.color == 2
    colormap(EC.edg.CM);
    cb=colorbar('location','South');
    if EC.edg.color_map_type == 1 % Modified by Mingrui, 20170303, support fixed color mapping
        if EC.edg.color_abs == 0
            caxis([min(surf.cylinder(:,3)),max(surf.cylinder(:,3))]);
        else
            caxis([min(abs(surf.cylinder(:,3))),max(abs(surf.cylinder(:,3)))]);
        end
    else
        caxis([EC.edg.color_map_low,EC.edg.color_map_high]);
    end
    set(cb,'Position',[0.4 0.055 0.2 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function PlotLine4v(surf,a,flag,centerX)
global cam
global EC
for i=1:6
    axes(a(i));
    hold on
    switch i
        case {1,3}
            if flag==1
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            end
        case{2,4}
            if flag==2
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,3,1);
                    end
                end
            end
        case 5
            if flag==1
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,3,2);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,3,2);
                    end
                end
            end
        case 6
            if flag==2
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)<centerX&&surf.sphere(surf.cylinder(j,2),1)<centerX
                        DrawLine(surf,j,3,2);
                    end
                end
            else
                for j=1:surf.ncyl
                    if surf.sphere(surf.cylinder(j,1),1)>centerX&&surf.sphere(surf.cylinder(j,2),1)>centerX
                        DrawLine(surf,j,3,2);
                    end
                end
            end
    end
    axis tight; axis vis3d off;eval(['material ',EC.glb.material,';']);eval(['lighting ',EC.glb.lighting,';']);
    cam(i) = camlight(EC.glb.lightdirection);
    hold off
end

% Add by Mingrui, 20150114, show colorbar when using colormap
if EC.edg.color == 2
    colormap(EC.edg.CM);
    cb=colorbar('location','South');
    if EC.edg.color_map_type == 1 % Modified by Mingrui, 20170303, support fixed color mapping
        if EC.edg.color_abs == 0
            caxis([min(surf.cylinder(:,3)),max(surf.cylinder(:,3))]);
        else
            caxis([min(abs(surf.cylinder(:,3))),max(abs(surf.cylinder(:,3)))]);
        end
    else
        caxis([EC.edg.color_map_low,EC.edg.color_map_high]);
    end
    set(cb,'Position',[0.4 0.6 0.2 0.03]);
    tmp = version;
    ind = find(tmp == '.');
    if str2double(tmp(1:ind(2)-1))<8.4
        set(cb,'XAxisLocation','bottom');
        set(cb,'XTick',get(cb,'XLim'));
    else
        cb.AxisLocation = 'out';
        cb.Ticks = cb.Limits;
    end
end

function DrawLine(surf,i,nv,iv)
global EC
%length_cyl=norm(surf.sphere(surf.cylinder(i,2),:)-surf.sphere(surf.cylinder(i,1),:)); % Fixed a bug by Mingrui Xia, 20120411, the length calculat was wrong.

if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
    sphere = surf.sphere2;
else
    sphere = surf.sphere;
end

length_cyl=norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))-sphere(surf.cylinder(i,2),7)-sphere(surf.cylinder(i,1),7);
switch EC.glb.detail % Add by Mingrui Xia, 20120413, adjust graph detail
    case 1
        det = 5;
    case 2
        det = 10;
    case 3
        det = 20;
        
end
% switch EC.edg.size
%     case 1
%         n=EC.edg.size_size;
%     case 2
%         switch EC.edg.size_value
%             case 1
%                 n=surf.cylinder(i,3);
%             case 2
%                 switch EC.edg.size_abs
%                     case 0
%                         n=surf.cylinder2(i,3);
%                     case 1
%                         n=abs(surf.cylinder2(i,3));
%                 end
%                 if n<0.1
%                     n=0.2;
%                 elseif n>3;
%                     n=3;
%                 end
%         end
%     case 3
%         switch EC.edg.size_abs
%             case 0
%                 if surf.cylinder2(i,3)>EC.edg.size_threshold
%                     switch EC.edg.size_value
%                         case 1
%                             n=surf.cylinder(i,3);
%                         case 2
%                             n=surf.cylinder2(i,3);
%                     end
%                     % n = 5;
%                 else
%                     n=1;
%                 end
%
%             case 1
%                 if surf.cylinder2(i,3)>EC.edg.size_threshold||surf.cylinder2(i,3)<-EC.edg.size_threshold
%                     switch EC.edg.size_value
%                         case 1
%                             n=surf.cylinder(i,3);
%                         case 2
%                             n=surf.cylinder2(i,3);
%                     end
%                 else
%                     n=1;
%                 end
%         end
% end
n = surf.cylinder(i,4) * EC.edg.size_ratio;
theta = (0:det) / det * 2 * pi;
sintheta = sin(theta);
sintheta(det + 1) = 0;

if EC.edg.directed == 1% Add by Mingrui Xia, 20120621, draw directed network.
    % Modified by Mingrui, 20170309, draw two edges for bidirected network
    %     if surf.cylinder(i,7) == 1
    %         n = [linspace(0,1.5,9),ones(1,round(length_cyl-9)) * 1]' * n;
    %     else
    %         n = [linspace(0,2,9),ones(1,round(length_cyl-18)) * 1,linspace(2,0,9)]' * n;
    %     end
    ind = find(surf.cylinder(:,1)==surf.cylinder(i,2)&surf.cylinder(:,2)==surf.cylinder(i,1));
    if isempty(ind)
        n = [linspace(0,2,9),ones(1,round(length_cyl-9)) * 1]' * n;
        x = n * cos(theta);
        y = n * sintheta;
        w = length(n);
        z = (0:w-1)'/(w-1) * ones(1,det + 1);
        Line = mesh(x,y,z * length_cyl);
        unit_Vx=[0 0 1];
        angle_X1X2=acos( dot( unit_Vx,sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3) )/( norm(unit_Vx)*norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))) )*180/pi;
        axis_rot=cross([0 0 1],(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3)) );
        if angle_X1X2~=0 && angle_X1X2~=180 % Modified by Mingrui, 20160918, fix the bug when angle is 180 degree
            rotate(Line,axis_rot,angle_X1X2,[0 0 0])
        end
        set(Line,'XData',get(Line,'XData')+sphere(surf.cylinder(i,1),1) + (sphere(surf.cylinder(i,2),1) -sphere(surf.cylinder(i,1),1))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
        set(Line,'YData',get(Line,'YData')+sphere(surf.cylinder(i,1),2) + (sphere(surf.cylinder(i,2),2) -sphere(surf.cylinder(i,1),2))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
        if angle_X1X2~=180
            set(Line,'ZData',get(Line,'ZData')+sphere(surf.cylinder(i,1),3) + (sphere(surf.cylinder(i,2),3) -sphere(surf.cylinder(i,1),3))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
        else
            set(Line,'ZData',get(Line,'ZData')+sphere(surf.cylinder(i,2),3) + (sphere(surf.cylinder(i,1),3) -sphere(surf.cylinder(i,2),3))/norm(sphere(surf.cylinder(i,1),1:3)-sphere(surf.cylinder(i,2),1:3))*sphere(surf.cylinder(i,1),7));
        end
        
        set(Line,'FaceColor',EC.edg.CM(surf.cylinder(i,5),:));
        set(Line,'EdgeColor','none');
        set(Line,'FaceAlpha',surf.cylinder(i,6));
        set(Line,'EdgeAlpha',0);
    elseif ind > i;
        n1 = [linspace(0,2,9),ones(1,round(length_cyl-9)) * 1]' * n;
        n2 = n1;
        
        x1 = n1 * cos(theta);
        y1 = n1 * sintheta;
        w1 = length(n1);
        z1 = (0:w1 - 1)'/(w1 - 1) * ones(1,det + 1);
        Line1 = mesh(x1,y1,z1 * length_cyl);
        unit_Vx = [0 0 1];
        angle_X1X2 = acos(dot(unit_Vx,sphere(surf.cylinder(i,2),1:3) - sphere(surf.cylinder(i,1),1:3) )/( norm(unit_Vx)*norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))) )*180/pi;
        axis_rot = cross([0 0 1],(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3)) );
        if angle_X1X2~=0 && angle_X1X2~=180 % Modified by Mingrui, 20160918, fix the bug when angle is 180 degree
            rotate(Line1,axis_rot,angle_X1X2,[0 0 0])
        end
        set(Line1,'XData',get(Line1,'XData')+sphere(surf.cylinder(i,1),1) + (sphere(surf.cylinder(i,2),1) -sphere(surf.cylinder(i,1),1))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
        set(Line1,'YData',get(Line1,'YData')+sphere(surf.cylinder(i,1),2) + (sphere(surf.cylinder(i,2),2) -sphere(surf.cylinder(i,1),2))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
        if angle_X1X2~=180
            set(Line1,'ZData',get(Line1,'ZData')+sphere(surf.cylinder(i,1),3) + (sphere(surf.cylinder(i,2),3) -sphere(surf.cylinder(i,1),3))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
        else
            set(Line1,'ZData',get(Line1,'ZData')+sphere(surf.cylinder(i,2),3) + (sphere(surf.cylinder(i,1),3) -sphere(surf.cylinder(i,2),3))/norm(sphere(surf.cylinder(i,1),1:3)-sphere(surf.cylinder(i,2),1:3))*sphere(surf.cylinder(i,1),7));
        end
        set(Line1,'FaceColor',EC.edg.CM(surf.cylinder(i,5),:));
        set(Line1,'EdgeColor','none');
        set(Line1,'FaceAlpha',surf.cylinder(i,6));
        set(Line1,'EdgeAlpha',0);
        
        x2 = n2 * cos(theta);
        y2 = n2 * sintheta;
        w2 = length(n2);
        z2 = (0:w2 - 1)'/(w2 - 1) * ones(1,det + 1);
        Line2 = mesh(x2,y2,z2 * length_cyl);
        unit_Vx = [0 0 1];
        angle_X1X2 = acos(dot(unit_Vx,sphere(surf.cylinder(ind,2),1:3) - sphere(surf.cylinder(ind,1),1:3) )/( norm(unit_Vx)*norm(sphere(surf.cylinder(ind,2),1:3)-sphere(surf.cylinder(ind,1),1:3))) )*180/pi;
        axis_rot = cross([0 0 1],(sphere(surf.cylinder(ind,2),1:3)-sphere(surf.cylinder(ind,1),1:3)) );
        if angle_X1X2~=0 && angle_X1X2~=180 % Modified by Mingrui, 20160918, fix the bug when angle is 180 degree
            rotate(Line2,axis_rot,angle_X1X2,[0 0 0])
        end
        set(Line2,'XData',get(Line2,'XData')+sphere(surf.cylinder(ind,1),1) + (sphere(surf.cylinder(ind,2),1) -sphere(surf.cylinder(ind,1),1))/norm(sphere(surf.cylinder(ind,2),1:3)-sphere(surf.cylinder(ind,1),1:3))*sphere(surf.cylinder(ind,1),7));
        set(Line2,'YData',get(Line2,'YData')+sphere(surf.cylinder(ind,1),2) + (sphere(surf.cylinder(ind,2),2) -sphere(surf.cylinder(ind,1),2))/norm(sphere(surf.cylinder(ind,2),1:3)-sphere(surf.cylinder(ind,1),1:3))*sphere(surf.cylinder(ind,1),7));
        if angle_X1X2~=180
            set(Line2,'ZData',get(Line2,'ZData')+sphere(surf.cylinder(ind,1),3) + (sphere(surf.cylinder(ind,2),3) -sphere(surf.cylinder(ind,1),3))/norm(sphere(surf.cylinder(ind,2),1:3)-sphere(surf.cylinder(ind,1),1:3))*sphere(surf.cylinder(ind,1),7));
        else
            set(Line2,'ZData',get(Line2,'ZData')+sphere(surf.cylinder(ind,2),3) + (sphere(surf.cylinder(ind,1),3) -sphere(surf.cylinder(ind,2),3))/norm(sphere(surf.cylinder(ind,1),1:3)-sphere(surf.cylinder(ind,2),1:3))*sphere(surf.cylinder(ind,1),7));
        end
        set(Line2,'FaceColor',EC.edg.CM(surf.cylinder(ind,5),:));
        set(Line2,'EdgeColor','none');
        set(Line2,'FaceAlpha',surf.cylinder(ind,6));
        set(Line2,'EdgeAlpha',0);
        
        edgedirect = [sphere(surf.cylinder(ind,1),1) - sphere(surf.cylinder(ind,2),1), sphere(surf.cylinder(ind,1),2) - sphere(surf.cylinder(ind,2),2),sphere(surf.cylinder(ind,1),3) - sphere(surf.cylinder(ind,2),3)];
        switch nv
            case 1
                switch iv
                    case 1
                        normdirect = [edgedirect(2),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 2
                        normdirect = [edgedirect(2),edgedirect(1)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'XData',get(Line1,'XData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 3
                        normdirect = [edgedirect(1),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'XData',get(Line1,'XData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                end
            case 2
                normdirect = [edgedirect(2),edgedirect(3)];
                if sum(normdirect(:))~=0
                movelength = null(normdirect/norm(normdirect));
                set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                end
            case 3
                switch iv
                    case 1
                        normdirect = [edgedirect(2),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 2
                        normdirect = [edgedirect(2),edgedirect(1)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'XData',get(Line1,'XData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                end
            case 4
                switch iv
                    case 1
                        normdirect = [edgedirect(2),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 2
                        normdirect = [edgedirect(2),edgedirect(1)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'XData',get(Line1,'XData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                end
            case 5
                switch iv
                    case 1
                        normdirect = [edgedirect(2),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 2
                        normdirect = [edgedirect(2),edgedirect(1)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'XData',get(Line1,'XData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 3
                        normdirect = [edgedirect(1),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'XData',get(Line1,'XData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                       
                end
            case 6
                switch iv
                    case {1,0}
                        normdirect = [edgedirect(2),edgedirect(3)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'ZData',get(Line1,'ZData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'ZData',get(Line2,'ZData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                    case 2
                        normdirect = [edgedirect(2),edgedirect(1)];
                        if sum(normdirect(:))~=0
                        movelength = null(normdirect/norm(normdirect));
                        set(Line1,'YData',get(Line1,'YData')+ movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line1,'XData',get(Line1,'XData')+ movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'YData',get(Line2,'YData')- movelength(1)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        set(Line2,'XData',get(Line2,'XData')- movelength(2)*surf.cylinder(i,4) * EC.edg.size_ratio * 1.3);
                        end
                end
        end
        
    end
else
    n = ones(100,1) * 0.5 * n ;
    x = n * cos(theta);
    y = n * sintheta;
    w = length(n);
    z = (0:w-1)'/(w-1) * ones(1,det + 1);
    Line = mesh(x,y,z * length_cyl);
    unit_Vx=[0 0 1];
    angle_X1X2=acos( dot( unit_Vx,sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3) )/( norm(unit_Vx)*norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))) )*180/pi;
    axis_rot=cross([0 0 1],(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3)) );
    if angle_X1X2~=0 && angle_X1X2~=180 % Modified by Mingrui, 20160918, fix the bug when angle is 180 degree
        rotate(Line,axis_rot,angle_X1X2,[0 0 0])
    end
    set(Line,'XData',get(Line,'XData')+sphere(surf.cylinder(i,1),1) + (sphere(surf.cylinder(i,2),1) -sphere(surf.cylinder(i,1),1))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
    set(Line,'YData',get(Line,'YData')+sphere(surf.cylinder(i,1),2) + (sphere(surf.cylinder(i,2),2) -sphere(surf.cylinder(i,1),2))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
    if angle_X1X2~=180
        set(Line,'ZData',get(Line,'ZData')+sphere(surf.cylinder(i,1),3) + (sphere(surf.cylinder(i,2),3) -sphere(surf.cylinder(i,1),3))/norm(sphere(surf.cylinder(i,2),1:3)-sphere(surf.cylinder(i,1),1:3))*sphere(surf.cylinder(i,1),7));
    else
        set(Line,'ZData',get(Line,'ZData')+sphere(surf.cylinder(i,2),3) + (sphere(surf.cylinder(i,1),3) -sphere(surf.cylinder(i,2),3))/norm(sphere(surf.cylinder(i,1),1:3)-sphere(surf.cylinder(i,2),1:3))*sphere(surf.cylinder(i,1),7));
    end
    % switch EC.edg.color
    %     case 1
    %         ci=1;
    %     case 2
    %         switch EC.edg.color_abs %%% Edited by Mingrui Xia, 20120413, fix colormap
    %             %         for edges
    %             case 0
    %                 %                 ci=int32(45*surf.cylinder3(i,3)-3.5);
    %                 k = 63/(max(surf.cylinder2(:,3)) - min(surf.cylinder2(:,3)));
    %                 b = 64 - 63/(max(surf.cylinder2(:,3)) - min(surf.cylinder2(:,3))) *...
    %                     max(surf.cylinder2(:,3));
    %                 ci = round(k*surf.cylinder2(i,3) + b);
    %                 if ci<1
    %                     ci=1;
    %                 elseif ci>64
    %                     ci=64;
    %                 end
    %             case 1
    %                 k = 63/(max(abs(surf.cylinder2(:,3))) - min(abs(surf.cylinder2(:,3))));
    %                 b = 64 - 63/(max(abs(surf.cylinder2(:,3))) - min(abs(surf.cylinder2(:,3))))...
    %                     * max(abs(surf.cylinder2(:,3)));
    %                 ci = round(k*abs(surf.cylinder2(i,3)) + b);
    %                 %                 ci=int32(45*surf.cylinder3(i,3)-3.5);
    %                 if ci<1
    %                     ci=1;
    %                 elseif ci>64
    %                     ci=64;
    %                 end
    %         end
    %     case 3
    %         switch EC.edg.color_abs
    %             case 0
    %                 if surf.cylinder2(i,3)>EC.edg.color_threshold
    %                     ci=1;
    %                 else
    %                     ci=64;
    %                 end
    %             case 1
    %                 if surf.cylinder2(i,3)>EC.edg.color_threshold||surf.cylinder2(i,3)<-EC.edg.color_threshold
    %                     ci=1;
    %                 else
    %                     ci=64;
    %                 end
    %         end
    %     case 4
    %         if length_cyl>EC.edg.color_distance
    %             ci=1;
    %         else
    %             ci=64;
    %         end
    %     case 5 % Added by Mingrui Xia, 20120809 add edge color according to nodal module
    %         if sphere(surf.cylinder(i,1),4) == sphere(surf.cylinder(i,2),4)
    %             ci = find(EC.nod.ModularNumber == sphere(surf.cylinder(i,1),4));
    %         else
    %             ci = 21;
    %         end
    %
    % end
    set(Line,'FaceColor',EC.edg.CM(surf.cylinder(i,5),:));
    set(Line,'EdgeColor','none');
    set(Line,'FaceAlpha',surf.cylinder(i,6));
    set(Line,'EdgeAlpha',0);
end

function DoubleMeshPrepare
% Added by Mingrui Xia, 20120717, show two brains in one figure
global surf
global EC
surf.coord2 = surf.coord;
switch EC.lot.view_direction
    case 1
        width = 0.1 * (max(surf.coord(2,:)) - min(surf.coord(2,:)));
        %         surf.coord2(1,:) = -surf.coord2(1,:);
        %         surf.coord2(2,:) = -surf.coord(2,:) + 2 * min(surf.coord(2,:)) - width;
        surf.coord2(1,:) = -surf.coord2(1,:);
        surf.coord2(2,:) = -surf.coord(2,:) + 2 * min(surf.coord(2,:)) - width;
    case 2
        surf.coord2(1,:) = surf.coord(1,:) + 1.1 * (max(surf.coord(1,:)) - min(surf.coord(1,:)));
    case 3
        surf.coord2(1,:) = surf.coord(1,:) - 1.1 * (max(surf.coord(1,:)) - min(surf.coord(1,:)));
end

function DoubleNodePrepare
% Added by Mingrui Xia, 20120717, show two brains in one figure
global surf
global EC
surf.sphere2 = surf.sphere;
switch EC.lot.view_direction
    case 1
        width = 0.1 * (max(surf.coord(2,:)) - min(surf.coord(2,:)));
        %         surf.sphere2(floor(surf.nsph/2) + 1:end,1) = -surf.sphere2(floor(surf.nsph/2) + 1:end,1);
        %         surf.sphere2(floor(surf.nsph/2) + 1:end,2) = -surf.sphere2(floor(surf.nsph/2) + 1:end,2) + 2 * min(surf.coord(2,:)) - width;
        surf.sphere2(floor(surf.nsph/2) + 1:end,1) = -surf.sphere2(floor(surf.nsph/2) + 1:end,1);
        surf.sphere2(floor(surf.nsph/2) + 1:end,2) = -surf.sphere2(floor(surf.nsph/2) + 1:end,2) + 2 * min(surf.coord(2,:)) - width;
    case 2
        surf.sphere2(floor(surf.nsph/2) + 1:end,1) = surf.sphere2(floor(surf.nsph/2) + 1:end,1) + 1.1 * (max(surf.coord(1,:)) - min(surf.coord(1,:)));
    case 3
        surf.sphere2(floor(surf.nsph/2) + 1:end,1) = surf.sphere2(floor(surf.nsph/2) + 1:end,1) - 1.1 * (max(surf.coord(1,:)) - min(surf.coord(1,:)));
end

function BoundaryPrepare
%Added by Mingrui, 20190327, draw boundary
global surf
global EC
switch EC.msh.boundary
    case 2
        if isfield(surf,'T')
            mask = surf.T;
            if sum(surf.T-fix(surf.T))~=0
                mask(mask~=0) = 1;
            end
            facevalue = zeros(size(surf.tri));
            facevalue(:) = mask(surf.tri(:));
            facevalue_judge = 1 - all(diff(facevalue,[],2)==0,2);
            face_remove = surf.tri(facevalue_judge==1,:);
            surf.boundary_value = zeros(surf.vertex_number,1);
            surf.boundary_value(unique(face_remove)) = 1;
            surf.boundary_value(surf.T==0) = 0;
        end
    case 3
        if isfield(surf,'T')
            mask = zeros(surf.vertex_number,1);
            mask(surf.T>EC.msh.boundary_value) = 1;
            facevalue = zeros(size(surf.tri));
            facevalue(:) = mask(surf.tri(:));
            facevalue_judge = 1 - all(diff(facevalue,[],2)==0,2);
            face_remove = surf.tri(facevalue_judge==1,:);
            surf.boundary_value = zeros(surf.vertex_number,1);
            surf.boundary_value(unique(face_remove)) = 1;
            surf.boundary_value(surf.T==0) = 0;
        end
    case 4
        surf.boundary_value = EC.msh.boundary_value;
end

function a=FileView
global FLAG
global surf
global EC

set(gcf,'Renderer',EC.glb.render);
a=[];
if FLAG.Loadfile==0
    H=errordlg('Please read files first!');
    uiwait(H);
    if exist('BrainNet_background.tif','file')==2
        imshow(imread('BrainNet_background.tif'));
    end
else
    switch EC.lot.view
        case 1
            switch FLAG.Loadfile
                case 1
                    if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
                        DoubleMeshPrepare;
                    end
                    BoundaryPrepare;
                    a=PlotMesh1(surf,0);
                case 2
                    NodePrepare;
                    if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
                        DoubleNodePrepare;
                    end
                    a=PlotNode1(surf,a,2,EC);
                case 3
                    NodePrepare;
                    if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
                        DoubleMeshPrepare;
                        DoubleNodePrepare;
                    end
                    BoundaryPrepare;
                    a=PlotMesh1(surf,1);
                    a=PlotNode1(surf,a,2,EC);
                case 6
                    NodePrepare;
                    NetPrepare;
                    if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
                        DoubleNodePrepare;
                    end
                    a=PlotNode1(surf,a,1,EC);
                    PlotLine1(surf);
                case 7
                    NodePrepare;
                    NetPrepare;
                    if EC.msh.doublebrain == 1 % Added by Mingrui Xia, 20120717, show two brains in one figure
                        DoubleMeshPrepare;
                        DoubleNodePrepare;
                    end
                    BoundaryPrepare;
                    a=PlotMesh1(surf,1);
                    a=PlotNode1(surf,a,1,EC);
                    PlotLine1(surf);
                case 9
                    if EC.vol.type == 1
                        if FLAG.MAP==2
                            MapPrepare;
                        end %%% Edited by Mingrui Xia,111027, move FLAG.IsCalledByREST judgement into function MapCMPrepare.
                        %                     if isfield(FLAG,'IsCalledByREST') && FLAG.IsCalledByREST==1 % YAN Chao-Gan 111023. REST will define the colormap from outside
                        %                         low=EC.vol.nx;high=EC.vol.px;
                        %                     else
                        [low high]=MapCMPrepare;
                        %                     end
                        BoundaryPrepare;
                        a=MapMesh1(surf,low,high,0);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        a = PlotMesh1(surf,1);
                        a = PlotROI1(fv,a,0);
                    end
                case 11 %%% Added by Mingrui Xia, 20111116, add 11 for volume and node mode
                    NodePrepare;
                    if EC.vol.type == 1
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        [low high]=MapCMPrepare;
                        a=MapMesh1(surf,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        a = PlotMesh1(surf,1);
                        a = PlotROI1(fv,a,1);
                    end
                    a=PlotNode1(surf,a,2,EC);
                case 15 %%% Added by Mingrui Xia, 20120210, add 15 for volume, node and edge mode
                    NodePrepare;
                    NetPrepare;
                    if EC.vol.type == 1
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        [low high]=MapCMPrepare;
                        a=MapMesh1(surf,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        a = PlotMesh1(surf,1);
                        a = PlotROI1(fv,a,1);
                    end
                    a=PlotNode1(surf,a,1,EC);
                    PlotLine1(surf);
            end
        case 2
            switch FLAG.Loadfile
                case 1
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    if cut<t
                        Viewer=GenerateView8(h1,w1);
                        Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh8(Viewer,Surfmatrix,0);
                    else
                        Viewer=GenerateView6(w1);
                        a=PlotMesh6(Viewer,surf,0);
                    end
                case 2
                    NodePrepare;
                    w1=NoMeshPage(surf);
                    Viewer=GenerateView6(w1);
                    a=PlotNode6(Viewer,surf,0,2,EC);
                case 3
                    NodePrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    if cut<t
                        [centerX,flag]=JudgeNode(surf,vl,vr);
                        Viewer=GenerateView8(h1,w1);
                        Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh8(Viewer,Surfmatrix,1);
                        PlotNode8(surf,a,flag,centerX,2,EC);
                    else
                        Viewer=GenerateView6(w1);
                        a=PlotMesh6(Viewer,surf,1);
                        PlotNode6(Viewer,surf,a,2,EC);
                    end
                case 6
                    NodePrepare;
                    NetPrepare;
                    w1=NoMeshPage(surf);
                    Viewer=GenerateView6(w1);
                    a=PlotNode6(Viewer,surf,0,1,EC);
                    PlotLine6(surf,a);
                case 7
                    NodePrepare;
                    NetPrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    if cut<t
                        [centerX,flag]=JudgeNode(surf,vl,vr);
                        Viewer=GenerateView8(h1,w1);
                        Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh8(Viewer,Surfmatrix,1);
                        PlotNode8(surf,a,flag,centerX,1,EC);
                        PlotLine8(surf,a,flag,centerX);
                    else
                        Viewer=GenerateView6(w1);
                        a=PlotMesh6(Viewer,surf,1);
                        PlotNode6(Viewer,surf,a,1,EC);
                        PlotLine6(surf,a);
                    end
                case 9
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);  %%% Edited by Mingrui Xia,111027, move FLAG.IsCalledByREST judgement into function MapCMPrepare.
                    %                     if isfield(FLAG,'IsCalledByREST') && FLAG.IsCalledByREST==1 % YAN Chao-Gan 111023. REST will define the colormap from outside
                    %                         low=EC.vol.nx;high=EC.vol.px;
                    %                     else
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        %                     end
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        if cut<t
                            Viewer=GenerateView8(h1,w1);
                            Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                            a=MapMesh8(Viewer,Surfmatrix,low,high,0);
                        else
                            Viewer=GenerateView6(w1);
                            a=MapMesh6(Viewer,surf,low,high,0);
                        end
                    else
                        BoundaryPrepare;
                        fv = ROIPrepare;
                        if cut<t
                            Viewer=GenerateView8(h1,w1);
                            Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                            a=PlotMesh8(Viewer,Surfmatrix,1);
                            a = PlotROI8(fv,a,2);
                        else
                            Viewer=GenerateView6(w1);
                            a=PlotMesh6(Viewer,surf,1);
                            a = PlotROI6(fv,a,2);
                        end
                    end
                case 11 %%% Added by Mingrui Xia, 20111116, add 11 for volume and node mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    NodePrepare;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        if cut<t
                            [centerX,flag]=JudgeNode(surf,vl,vr);
                            Viewer=GenerateView8(h1,w1);
                            Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                            a=MapMesh8(Viewer,Surfmatrix,low,high,1);
                            PlotNode8(surf,a,flag,centerX,2,EC);
                        else
                            Viewer=GenerateView6(w1);
                            a=MapMesh6(Viewer,surf,low,high,1);
                            PlotNode6(Viewer,surf,a,2,EC);
                        end
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        if cut<t
                            [centerX,flag]=JudgeNode(surf,vl,vr);
                            Viewer=GenerateView8(h1,w1);
                            Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                            a=PlotMesh8(Viewer,Surfmatrix,1);
                            a = PlotROI8(fv,a,1);
                            PlotNode8(surf,a,flag,centerX,2,EC);
                        else
                            Viewer=GenerateView6(w1);
                            a=PlotMesh6(Viewer,surf,1);
                            a = PlotROI6(fv,a,1);
                            PlotNode6(Viewer,surf,a,2,EC);
                        end
                    end
                case 15 %%% Added by Mingrui Xia, 20120210, add 15 for volume, node and edge mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    NodePrepare;
                    NetPrepare;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        if cut<t
                            [centerX,flag]=JudgeNode(surf,vl,vr);
                            Viewer=GenerateView8(h1,w1);
                            Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                            a=MapMesh8(Viewer,Surfmatrix,low,high,1);
                            PlotNode8(surf,a,flag,centerX,1,EC);
                            PlotLine8(surf,a,flag,centerX);
                        else
                            Viewer=GenerateView6(w1);
                            a=MapMesh6(Viewer,surf,low,high,1);
                            PlotNode6(Viewer,surf,a,1,EC);
                            PlotLine6(surf,a);
                        end
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        if cut<t
                            [centerX,flag]=JudgeNode(surf,vl,vr);
                            Viewer=GenerateView8(h1,w1);
                            Surfmatrix=GenertateSurfM8(surf,tl,tr,vl,vr,cuv);
                            a=PlotMesh8(Viewer,Surfmatrix,1);
                            a = PlotROI8(fv,a,1);
                            PlotNode8(surf,a,flag,centerX,1,EC);
                            PlotLine8(surf,a,flag,centerX);
                        else
                            Viewer=GenerateView6(w1);
                            a=PlotMesh6(Viewer,surf,1);
                            a = PlotROI6(fv,a,1);
                            PlotNode6(Viewer,surf,a,1,EC);
                            PlotLine6(surf,a);
                        end
                    end
            end
        case 3  %%% YAN Chao-Gan 111023 Added BEGIN. For medium view (4 views). Only finished for volume view mode
            switch FLAG.Loadfile
                case 9
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf); %%% Edited by Mingrui Xia,111027, move FLAG.IsCalledByREST judgement into function MapCMPrepare.
                    %                     if isfield(FLAG,'IsCalledByREST') && FLAG.IsCalledByREST==1 % YAN Chao-Gan 111023. REST will define the colormap from outside
                    %                         low=EC.vol.nx;high=EC.vol.px;
                    %                     else
                    
                    Viewer=GenerateView4;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        %                     end
                        if FLAG.MAP == 2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                        a=MapMesh4(Viewer,Surfmatrix,low,high,0);%%% YAN Chao-Gan 111023 Added END. %%%
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh4(Viewer,Surfmatrix,1);
                        a = PlotROI4(fv,a,2);
                    end
                case 1 %%% Mingrui Xia 110226 Added below.
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    Viewer=GenerateView4;
                    BoundaryPrepare;
                    Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                    a=PlotMesh4(Viewer,Surfmatrix,0);
                case 3
                    NodePrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    Viewer=GenerateView4;
                    Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                    a=PlotMesh4(Viewer,Surfmatrix,1);
                    PlotNode4(surf, a, flag, centerX, 2, EC);
                case 7
                    NodePrepare;
                    NetPrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    Viewer=GenerateView4;
                    Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                    a=PlotMesh4(Viewer,Surfmatrix,1);
                    PlotNode4(surf,a,flag,centerX,1,EC);
                    PlotLine4(surf,a,flag,centerX);
                case 11 %%% Added by Mingrui Xia, 20111116, add 11 for volume and node mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    NodePrepare;
                    Viewer=GenerateView4;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                        a=MapMesh4(Viewer,Surfmatrix,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh4(Viewer,Surfmatrix,1);
                        a = PlotROI4(fv,a,1);
                    end
                    PlotNode4(surf, a, flag, centerX, 2, EC);
                    
                case 15 %%% Added by Mingrui Xia, 20120210, add 15 for volume, node and edge mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    NodePrepare;
                    NetPrepare;
                    Viewer=GenerateView4;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                        a=MapMesh4(Viewer,Surfmatrix,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh4(Viewer,Surfmatrix,1);
                        a = PlotROI4(fv,a,1);
                    end
                    PlotNode4(surf, a, flag, centerX, 1, EC);
                    PlotLine4(surf,a,flag,centerX);
            end
        case 4  %%% Lijie Huang 130221 Added. Modified based on YAN Chao-Gan's `case 3` version. For medium view (4 views) and 2 ventral views. Only finished for volume view mode
            switch FLAG.Loadfile
                case 9
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf); %%% Edited by Mingrui Xia,111027, move FLAG.IsCalledByREST judgement into function MapCMPrepare.
                    Viewer=GenerateView4v;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        %                     end
                        if FLAG.MAP == 2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                        a=MapMesh4v(Viewer,Surfmatrix,low,high,0);%%% YAN Chao-Gan 111023 Added END. %%%
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh4v(Viewer,Surfmatrix,1);
                        a = PlotROI4v(fv,a,2);
                    end
                    
                case 1 %%% Mingrui Xia 110226 Added below.
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    Viewer=GenerateView4v;
                    BoundaryPrepare;
                    Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                    a=PlotMesh4v(Viewer,Surfmatrix,0);
                case 3
                    NodePrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    Viewer=GenerateView4v;
                    Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                    a=PlotMesh4v(Viewer,Surfmatrix,1);
                    PlotNode4v(surf, a, flag, centerX, 2, EC);
                case 7
                    NodePrepare;
                    NetPrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    Viewer=GenerateView4v;
                    Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                    a=PlotMesh4v(Viewer,Surfmatrix,1);
                    PlotNode4v(surf,a,flag,centerX,1,EC);
                    PlotLine4v(surf,a,flag,centerX);
                case 11 %%% Added by Mingrui Xia, 20111116, add 11 for volume and node mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    NodePrepare;
                    Viewer=GenerateView4v;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                        a=MapMesh4v(Viewer,Surfmatrix,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh4v(Viewer,Surfmatrix,1);
                        a = PlotROI4v(fv,a,1);
                    end
                    PlotNode4v(surf, a, flag, centerX, 2, EC);
                    
                case 15 %%% Added by Mingrui Xia, 20120210, add 15 for volume, node and edge mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv]=CutMesh(surf);
                    [centerX,flag]=JudgeNode(surf,vl,vr);
                    NodePrepare;
                    NetPrepare;
                    Viewer=GenerateView4v;
                    if EC.vol.type == 1
                        [low high]=MapCMPrepare;
                        if FLAG.MAP==2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                        a=MapMesh4v(Viewer,Surfmatrix,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix=GenertateSurfM4v(surf,tl,tr,vl,vr,cuv);
                        a=PlotMesh4v(Viewer,Surfmatrix,1);
                        a = PlotROI4v(fv,a,1);
                    end
                    PlotNode4v(surf, a, flag, centerX, 1, EC);
                    PlotLine4v(surf,a,flag,centerX);
            end
            
            % Added by Mingrui, 20170309, add layout for lateral, medial and dorsal view
        case 5
            switch FLAG.Loadfile
                case 1
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv] = CutMesh(surf);
                    Viewer = GenerateView5;
                    BoundaryPrepare;
                    Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                    a = PlotMesh5(Viewer,Surfmatrix,0);
                case 3
                    NodePrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv] = CutMesh(surf);
                    [centerX,flag] = JudgeNode(surf,vl,vr);
                    Viewer = GenerateView5;
                    Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                    a = PlotMesh5(Viewer,Surfmatrix,1);
                    PlotNode5(surf, a, flag, centerX, 2, EC);
                case 7
                    NodePrepare;
                    NetPrepare;
                    BoundaryPrepare;
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv] = CutMesh(surf);
                    [centerX,flag] = JudgeNode(surf,vl,vr);
                    Viewer = GenerateView5;
                    Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                    a = PlotMesh5(Viewer,Surfmatrix,1);
                    PlotNode5(surf,a,flag,centerX,1,EC);
                    PlotLine5(surf,a,flag,centerX);
                case 9
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv] = CutMesh(surf); %%% Edited by Mingrui Xia,111027, move FLAG.IsCalledByREST judgement into function MapCMPrepare.
                    Viewer = GenerateView5;
                    if EC.vol.type == 1
                        [low, high] = MapCMPrepare;
                        %                     end
                        if FLAG.MAP == 2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                        a = MapMesh5(Viewer,Surfmatrix,low,high,0);%%% YAN Chao-Gan 111023 Added END. %%%
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                        a = PlotMesh5(Viewer,Surfmatrix,1);
                        a = PlotROI5(fv,a,2);
                    end
                case 11 %%% Added by Mingrui Xia, 20111116, add 11 for volume and node mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv] = CutMesh(surf);
                    [centerX,flag] = JudgeNode(surf,vl,vr);
                    NodePrepare;
                    Viewer = GenerateView5;
                    if EC.vol.type == 1
                        [low, high] = MapCMPrepare;
                        if FLAG.MAP == 2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                        a = MapMesh5(Viewer,Surfmatrix,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                        a = PlotMesh5(Viewer,Surfmatrix,1);
                        a = PlotROI5(fv,a,1);
                    end
                    PlotNode5(surf, a, flag, centerX, 2, EC);
                case 15 %%% Added by Mingrui Xia, 20120210, add 15 for volume, node and edge mode
                    [t,tl,tr,vl,vr,h1,w1,cut,cuv] = CutMesh(surf);
                    [centerX,flag] = JudgeNode(surf,vl,vr);
                    NodePrepare;
                    NetPrepare;
                    Viewer = GenerateView5;
                    if EC.vol.type == 1
                        [low, high] = MapCMPrepare;
                        if FLAG.MAP == 2
                            MapPrepare;
                        end
                        BoundaryPrepare;
                        Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                        a = MapMesh5(Viewer,Surfmatrix,low,high,1);
                    else
                        fv = ROIPrepare;
                        BoundaryPrepare;
                        Surfmatrix = GenertateSurfM5(surf,tl,tr,vl,vr,cuv);
                        a = PlotMesh5(Viewer,Surfmatrix,1);
                        a = PlotROI5(fv,a,1);
                    end
                    PlotNode5(surf, a, flag, centerX, 1, EC);
                    PlotLine5(surf,a,flag,centerX);
            end
    end
end

function MapPrepare
global surf
global EC

% Added by Mingrui, 20140925, statistic for SPM or REST nifti files
% Edited by Mingrui, 20150302, remove judge for statistical map
if ~strcmp(surf.test,'No')
    vol_tmp = surf.vol;
    vol_tmp(vol_tmp < EC.vol.threshold & vol_tmp > -EC.vol.threshold) = 0;
    [L, num] = bwlabeln(vol_tmp,EC.vol.rmm);
    for x = 1:num
        theCurrentCluster = L == x;
        if length(find(theCurrentCluster)) < EC.vol.clustersize
            vol_tmp(logical(theCurrentCluster)) = 0;
        end
    end
else
    vol_tmp = surf.vol;
end

% Edited by Mingrui Xia, 20120726, selection for different mapping algorithm.
% 1 for Nearest Voxel
% 2 for Average Vertex
% 3 for Average Voxel (3x3)
% 4 for Gaussian
% 5 for Interpolated (default)
% 6 for Maximum Voxel
% 7 for Minimum Voxel
% 8 for Extremum Voxel
switch EC.vol.mapalgorithm
    case 1
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<1|index(1,:)>surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<1|index(2,:)>surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<1|index(3,:)>surf.hdr.dim(3)) = 1;
        index = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        surf.T(index~=1) = vol_tmp(index(index~=1));
    case 2
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<1|index(1,:)>surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<1|index(2,:)>surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<1|index(3,:)>surf.hdr.dim(3)) = 1;
        index = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        surf.T(index~=1) = vol_tmp(index(index~=1));
        tmpT = surf.T;
        for i = 1:surf.vertex_number
            [m,n] = find(surf.tri == i);
            neibour = unique(surf.tri(m,:));
            surf.T(i) = mean(tmpT(neibour));
        end
    case 3
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<1|index(1,:)>surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<1|index(2,:)>surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<1|index(3,:)>surf.hdr.dim(3)) = 1;
        index = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        kernal = ones(3,3,3)/27;
        tmpT = convn(vol_tmp,kernal,'same');
        surf.T(index~=1) = tmpT(index(index~=1));
    case 4
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<1|index(1,:)>surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<1|index(2,:)>surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<1|index(3,:)>surf.hdr.dim(3)) = 1;
        index = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        tmpT = smooth3(vol_tmp,'gaussian');
        surf.T(index~=1) = tmpT(index(index~=1));
    case 5
        surf.coord(4,:)=1;
        position = surf.hdr.mat\surf.coord;
        position(4,:) = [];
        index=round(position);
        surf.coord(4,:) = [];
        
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<=1|index(1,:)>=surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<=1|index(2,:)>=surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<=1|index(3,:)>=surf.hdr.dim(3)) = 1;
        
        index = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        for i = 1:surf.vertex_number
            if index(i)~=1
                cube = [floor(position(:,i))';ceil(position(:,i))'];
                portion = position(:,i)' - cube(1,:);
                cube(2,portion == 0) = cube(2,portion == 0) + 1;
                tmpT = vol_tmp(cube(1,1):cube(2,1),cube(1,2):cube(2,2),cube(1,3):cube(2,3));
                tmpT = (tmpT(:,:,2) - tmpT(:,:,1)) .* portion(3) + tmpT(:,:,1);
                tmpT = (tmpT(:,2) - tmpT(:,1)) .* portion(2) + tmpT(:,1);
                tmpT = (tmpT(2) - tmpT(1)) .* portion(1) + tmpT(1);
                surf.T(i) = tmpT;
            end
        end
    case 6
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<=1|index(1,:)>=surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<=1|index(2,:)>=surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<=1|index(3,:)>=surf.hdr.dim(3)) = 1;
        ind = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        for i = 1:surf.vertex_number
            if ind(i) ~=1
                tmpT = vol_tmp(index(1,i)-1:index(1,i)+1,index(2,i)-1:index(2,i)+1,index(3,i)-1:index(3,i)+1);
                surf.T(i) = max(tmpT(:));
            end
        end
    case 7
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<=1|index(1,:)>=surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<=1|index(2,:)>=surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<=1|index(3,:)>=surf.hdr.dim(3)) = 1;
        ind = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        for i = 1:surf.vertex_number
            if ind(i) ~=1
                tmpT = vol_tmp(index(1,i)-1:index(1,i)+1,index(2,i)-1:index(2,i)+1,index(3,i)-1:index(3,i)+1);
                surf.T(i) = min(tmpT(:));
            end
        end
    case 8
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<=1|index(1,:)>=surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<=1|index(2,:)>=surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<=1|index(3,:)>=surf.hdr.dim(3)) = 1;
        ind = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        for i = 1:surf.vertex_number
            if ind(i) ~=1
                tmpT = vol_tmp(index(1,i)-1:index(1,i)+1,index(2,i)-1:index(2,i)+1,index(3,i)-1:index(3,i)+1);
                mx = max(tmpT(:));
                mn = min(tmpT(:));
                if mx>=abs(mn)
                    surf.T(i) = mx;
                else
                    surf.T(i) = mn;
                end
            end
        end
    case 9 % Added by Mingrui Xia, 20130104, Most neighbour voxel value.
        surf.coord(4,:)=1;
        index=round(surf.hdr.mat\surf.coord);
        surf.coord(4,:) = [];
        index(4,:) = [];
        surf.T=zeros(1,surf.vertex_number);
        index(:,index(1,:)<=1|index(1,:)>=surf.hdr.dim(1)) = 1;
        index(:,index(2,:)<=1|index(2,:)>=surf.hdr.dim(2)) = 1;
        index(:,index(3,:)<=1|index(3,:)>=surf.hdr.dim(3)) = 1;
        ind = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
        %          for i = 1:surf.vertex_number
        %             if ind(i) ~=1
        %                 tmpT = reshape(vol_tmp(index(1,i)-1:index(1,i)+1,index(2,i)-1:index(2,i)+1,index(3,i)-1:index(3,i)+1),1,[]);
        %                 X = unique(tmpT);
        %                 D = histc(tmpT,X);
        %                 Y = max(D);
        %                 surf.T(i) = X(find(D == Y,1));
        %             end
        %          end
        rad=2;
        %%%%% dilate
        for i = 1:surf.vertex_number
            if ind(i) ~=1
                tmpT = reshape(vol_tmp(index(1,i)-rad:index(1,i)+rad,index(2,i)-rad:index(2,i)+rad,index(3,i)-rad:index(3,i)+rad),1,[]);
                X = setdiff(unique(tmpT),0);
                if isempty(X),continue;end
                [Y,x] = max(histc(tmpT,X));
                surf.T(i)=X(x);
                %                 mn = histc(tmpT(:),unique(tmpT));
                %                 if length(find(tmpT==0)) && mn(1)>0.6*length(tmpT),surf.T(i)=0;end
            end
        end
        %         iteration = 20;
        %         for j = 1:iteration
        %         tmpT = surf.T;
        %         for i = 1:surf.vertex_number
        %             [m,n] = find(surf.tri == i);
        %             neibour = unique(surf.tri(m,:));
        %             surf.T(i) = median(tmpT(neibour));
        %         end
        %         end
end
% surf.T(isnan(surf.T)) = 0; % Added by Mingrui, 20140410, replace NaN to 0;
% for i=1:surf.vertex_number
%     surf.T(i)=vol_tmp(index(1,i),index(2,i),index(3,i));
% end

function [low high]=MapCMPrepare
global EC
global FLAG %%% Added by Mingrui Xia, 111027, judge if is called by REST.
if isfield(FLAG,'IsCalledByREST') && FLAG.IsCalledByREST==1
    %     switch EC.vol.display %% Edited by Mingrui Xia, 20120414, fix a bug of compatibility with REST.
    %         case 1
    low=EC.vol.nx;
    high=EC.vol.px;
    %         case 2
    %             low=EC.vol.pn;
    %             high=EC.vol.px;
    %         case 3
    %             low=EC.vol.nx;
    %             high=EC.vol.nn;
    %     end
else
    switch EC.vol.color_map
        %%% Edited by Mingrui Xia, 130624, arrange colormap order by name.
        case 13
            EC.vol.CMt=jet(1000);
        case 12
            EC.vol.CMt=hsv(1000);
        case 11
            EC.vol.CMt=hot(1000);
        case 7
            EC.vol.CMt=cool(1000);
        case 19
            EC.vol.CMt=spring(1000);
        case 20
            EC.vol.CMt=summer(1000);
        case 4
            EC.vol.CMt=autumn(1000);
        case 22
            EC.vol.CMt=winter(1000);
        case 10
            EC.vol.CMt=gray(1000);
        case 5
            EC.vol.CMt=bone(1000);
        case 8
            EC.vol.CMt=copper(1000);
        case 15
            EC.vol.CMt=pink(1000);
        case 14
            EC.vol.CMt=lines(1000);
        case 17
            EC.vol.CMt=rand(1000,3);
            col = randi(3,[1000,1]); %%% Added by Mingrui Xia, 20120425, make random colorbar bright.
            row = [1:1000]';
            ind = sub2ind([1000,3],row,col);
            EC.vol.CMt(ind) = 1;
        case 1 %%% Added by Mingrui Xia, 111027. seven new colorbar.
            EC.vol.CMt = [zeros(1,500),ones(1,500);linspace(1,0,500),linspace(0,1,500);ones(1,500),zeros(1,500)]';
        case 6
            EC.vol.CMt = colorcube(1000);
        case 9
            EC.vol.CMt = flag(1000);
        case 16
            EC.vol.CMt = prism(1000);
        case 21
            EC.vol.CMt = white(1000);
        case 3
            EC.vol.CMt = [ones(1,1000); linspace(0,1,1000); zeros(1,1000)]';
        case 2
            EC.vol.CMt = [zeros(1,1000); linspace(1,0,1000); ones(1,1000)]';
        case 23 % Added by Mingrui Xia, 20120113, Xjviewer negative colorbar.
            EC.vol.CMt = [zeros(1,1000); linspace(0,1,1000); ones(1,1000)]';
        case 18 % Added by Mingrui Xia, 20120806, spectral colorbar of SurfStat.
            base = [
                0.2000 0.2000 0.2000
                0.4667 0.0000 0.5333
                0.5333 0.0000 0.6000
                0.0000 0.0000 0.6667
                0.0000 0.0000 0.8667
                0.0000 0.4667 0.8667
                0.0000 0.6000 0.8667
                0.0000 0.6667 0.6667
                0.0000 0.6667 0.5333
                0.0000 0.6000 0.0000
                0.0000 0.7333 0.0000
                0.0000 0.8667 0.0000
                0.0000 1.0000 0.0000
                0.7333 1.0000 0.0000
                0.9333 0.9333 0.0000
                1.0000 0.8000 0.0000
                1.0000 0.6000 0.0000
                1.0000 0.0000 0.0000
                0.8667 0.0000 0.0000
                0.8000 0.0000 0.0000
                0.8000 0.8000 0.8000
                ];
            n = length(base);
            X0 = linspace (1, n, 1000);
            EC.vol.CMt = interp1(1:n,base,X0);
            %             EC.vol.CMt = spectral(1000);
            
            % Add by Mingrui Xia, add support for annot file
        case 25
            EC.vol.CMt = EC.vol.CM_annot;
        case 26
            EC.vol.CMt = [zeros(1,500),ones(1,500);linspace(0.8,68/255,500),linspace(68/255,1,500);ones(1,500),zeros(1,500)]';
    end
    if EC.vol.adjustCM == 1
        switch EC.vol.display
            case 1 %%% Edited by Mingrui Xia, 111025, change algorithm for color mapping.
                
                EC.vol.CM=AdjustColorMap(EC.vol.CMt,EC.vol.null,EC.vol.nx,EC.vol.nn,EC.vol.pn,EC.vol.px);
                low=EC.vol.nx;
                high=EC.vol.px;
            case 2
                low=EC.vol.pn;
                high=EC.vol.px;
                %                 ind = round(linspace(1, size(EC.vol.CMt, 1), 999));
                ind = floor(linspace(1, size(EC.vol.CMt, 1) + 0.9999, 999)); % Edited by Mingrui Xia, 20120724 fix colormap adjusting
                EC.vol.CM(2:1000,:)=EC.vol.CMt(ind,:);
                EC.vol.CM(1,:)=EC.vol.null;
            case 3
                low=EC.vol.nx;
                high=EC.vol.nn;
                %                 ind = round(linspace(1, size(EC.vol.CMt, 1), 999));
                ind = floor(linspace(1, size(EC.vol.CMt, 1) + 0.9999, 999));% Edited by Mingrui Xia, 20120724 fix colormap adjusting
                EC.vol.CM(1:999,:)=EC.vol.CMt(ind,:);
                EC.vol.CM(1000,:)=EC.vol.null;
        end
    else
        EC.vol.CM = EC.vol.CMt;
        switch EC.vol.display
            case 1
                low=EC.vol.nx;
                high=EC.vol.px;
            case 2
                low=EC.vol.pn;
                high=EC.vol.px;
            case 3
                low=EC.vol.nx;
                high=EC.vol.nn;
        end
    end
end

function NewColorMap=AdjustColorMap(OriginalColorMap,NullColor,NMax,NMin,PMin,PMax)
%%% Added by Mingrui Xia, 20111025, to adjust colormap.
% Adjust the colormap to leave blank to values under threshold, the orginal color map with be set into [NMax NMin] and [PMin PMax]. Written by YAN Chao-Gan, 111023
% Input: OriginalColorMap - the original color map
%        NullColor - The values between NMin and PMin will be set to this color (leave blank)
%        NMax, NMin, PMin, PMax - set the axis of colorbar (the orginal color map with be set into [NMax NMin] and [PMin PMax])
% Output: NewColorMap - the generated color map, a 1000 by 3 matrix.
TempColorMap = OriginalColorMap;
% OriginalColorMap = zeros(1000,3);
% ind = round(linspace(1,size(TempColorMap,1),1000));
ind = floor(linspace(1, size(TempColorMap, 1) + 0.9999, 1000)); % Edited by Mingrui Xia, 20121031 fix colormap adjusting
OriginalColorMap = TempColorMap(ind,:);
NewColorMap = repmat(NullColor,[1000 1]);
% ColorLen=size(OriginalColorMap,1);
% NegativeColorSegment = fix(1000*(NMin-NMax)/(PMax-NMax)/(ColorLen/2));
% for iColor=1:fix(ColorLen/2)
%     NewColorMap((iColor-1)*NegativeColorSegment+1:(iColor)*NegativeColorSegment,:) = repmat(OriginalColorMap(iColor,:),[NegativeColorSegment 1]);
% end
%
% PositiveColorSegment = fix(1000*(PMax-PMin)/(PMax-NMax)/(ColorLen/2));
% for iColor=ColorLen:-1:ceil(ColorLen/2+1)
%     NewColorMap(end-(ColorLen-iColor+1)*PositiveColorSegment+1:end-(ColorLen-iColor)*PositiveColorSegment,:) = repmat(OriginalColorMap(iColor,:),[PositiveColorSegment 1]);
% end
NegativeColorSegment = fix(1000*(NMin-NMax)/(PMax-NMax)); %%% Edited by Mingrui Xia 20111025, use linear sampling method.
if NegativeColorSegment == 0
    NegativeColorSegment = 1;
end
NegativeIndex = round(linspace(1,500,NegativeColorSegment));
NewColorMap(1:NegativeColorSegment,:) = OriginalColorMap(NegativeIndex,:);
PositiveColorSegment = fix(1000*(PMax-PMin)/(PMax-NMax));
if PositiveColorSegment == 0
    PositiveColorSegment = 1;
end
PositiveIndex = round(linspace(501,1000,PositiveColorSegment));
NewColorMap(end-PositiveColorSegment+1:end,:) = OriginalColorMap(PositiveIndex,:);

% --------------------------------------------------------------------
function NV_m_save_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EC
[filename,pathname]=uiputfile({'*.tif','TIFF 24-bit';'*.jpg','JPEG 24-bit';...
    '*.bmp','BMP 24-bit';'*.eps','EPS color';'*.png','PNG 24-bit';...
    '*.fig','Matlab Figure';'*.wrl','VRML'},'Save Image');
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    [pathstr, name, ext] = fileparts(fpath);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inch');
    set(gcf,'Paperposition',[1 1 EC.img.width/EC.img.dpi EC.img.height/EC.img.dpi]);
    switch ext
        case '.png'
            print(gcf,fpath,'-dpng',['-r',num2str(EC.img.dpi)])
        case '.tif'
            print(gcf,fpath,'-dtiff',['-r',num2str(EC.img.dpi)])
        case '.jpg'
            print(gcf,fpath,'-djpeg',['-r',num2str(EC.img.dpi)])
            %             saveas(gcf,fpath);
        case '.bmp'
            print(gcf,fpath,'-dbmp',['-r',num2str(EC.img.dpi)])
        case '.eps'
            print(gcf,fpath,'-depsc',['-r',num2str(EC.img.dpi)])
            
            % added by Mingrui Xia, save as matlab figure.
            % Modified by Mingrui, 20150112, using lowercase to be compatible with
            % higher version
        case '.fig'
            saveas(gcf,fpath,'fig');
            
            % added by Mingrui, 20170904, save as vrml format.
        case '.wrl'
            vrml(gcf,fpath);
    end
    msgbox('Image has saved!','Success','help');
end

% --------------------------------------------------------------------
function NV_m_es_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_es (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
global a
H=BrainNet_Option; %%% Edited by Mingrui Xia, 111027, changed option panel function name.
uiwait(H);
while FLAG.EC==2
    FLAG.EC=0;
    if ~isempty(a) && mean(ishandle(a))==1
        delete(a);
    end
    axes(handles.NV_axes);
    cla;
    a=FileView;
    uiwait(H);
end
if FLAG.EC==1
    FLAG.EC=0;
    if ~isempty(a) && mean(ishandle(a))==1
        delete(a);
    end
    axes(handles.NV_axes);
    cla;
    a=FileView;
end

% --------------------------------------------------------------------
function NV_m_LF_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_LF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
global a
H=BrainNet_LoadFiles;
uiwait(H);
if FLAG.LF==1
    %     if FLAG.Loadfile~=1
    H=BrainNet_Option;
    uiwait(H);
    while FLAG.EC==2
        FLAG.EC=3;
        if ~isempty(a) && mean(ishandle(a))==1
            delete(a);
        end
        axes(handles.NV_axes);
        cla;
        a=FileView;
        uiwait(H);
    end
    if FLAG.EC==1
        FLAG.EC=0;
        if ~isempty(a) && mean(ishandle(a))==1
            delete(a);
        end
        axes(handles.NV_axes);
        cla;
        a=FileView;
    end
    FLAG.LF=0;
end

% --------------------------------------------------------------------
function NV_m_Visualize_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_Visualize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NV_m_tools_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NV_m_others_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_others (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NV_m_cf_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_cf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global FLAG % Mingrui Xia 111029. For calling from outside.
% if isfield(FLAG,'IsCalledByREST') && FLAG.IsCalledByREST==1
%     handles.NV_axes=gca(hObject);
% end
axes(handles.NV_axes);
global a
if ~isempty(a) && mean(ishandle(a))==1
    delete(a);
    a=[];
end
cla;
if exist('BrainNet_background.tif','file')==2
    imshow(imread('BrainNet_background.tif'));
end

% --------------------------------------------------------------------
function NV_m_mm_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H=BrainNet_MergeMesh;
uiwait(H);

% --------------------------------------------------------------------
function NV_m_help_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NV_m_about_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'BrainNet Viewer 1.63 Released 20181219';'By Mingrui Xia';'mingruixia@gmail.com'},'About...','help');

% --------------------------------------------------------------------
% function NV_m_batch_Callback(hObject, eventdata, handles)
% % hObject    handle to NV_m_batch (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% global surf
% global FLAG
% global a
% global EC
% FLAG.Loadfile=3;
% surf=[];
% [filename,pathname]=uigetfile({'*.txt','Text Files (*.txt)'},'Select List File');
% fpath=fullfile(pathname,filename);
% fid=fopen(fpath);
% b=textscan(fid,'%s');
% fclose(fid);
% List=b{1,1};
%
% [filename2,pathname2]=uigetfile({'*.nv','NetViewer Files (*.nv)';'*.mesh','BrainVISA Mesh (*.mesh)';'*.pial','FreeSurfer Mesh (*.pial)';'*.*','All Files (*.*)'});
% fpath=fullfile(pathname2,filename2);
% fid=fopen(fpath);
% surf.vertex_number=fscanf(fid,'%f',1);
% surf.coord=fscanf(fid,'%f',[3,surf.vertex_number]);
% surf.ntri=fscanf(fid,'%f',1);
% surf.tri=fscanf(fid,'%d',[3,surf.ntri])';
% fclose(fid);
%
% [filename3,pathname3]=uiputfile({'*.TIF','TIFF 24-bit';'*.BMP','BMP 24-bit';'*.EPS','EPS color';'*.','JPEG 24-bit';'*.PNG','PNG 24-bit'},'Save Image');
% if isequal(filename,0)||isequal(pathname,0)
%     return;
% else
%     fpath=fullfile(pathname3,filename3);
%     [pathstr, name, ext] = fileparts(fpath);
% end
% H=BrainNet_Option;
% uiwait(H);
% for i=1:length(List)
%     load([pathname List{i} '.mat']);
%     surf.nsph=size(node,1);
%     surf.sphere=node;
%     fid=fopen([pathname List{i} '.txt']);
%     b=textscan(fid,'%s');
%     fclose(fid);
%     surf.label=b{1,1};
%     axes(handles.NV_axes);
%     cla;
%     if ~isempty(a)
%         delete(a);
%         a=[];
%     end
%     a=FileView;
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperUnits', 'inch');
%     set(gcf,'Paperposition',[1 1 EC.Img_sz(1)/EC.Img_sz(3) EC.Img_sz(2)/EC.Img_sz(3)]);
%     switch ext
%         case '.PNG'
%             print(gcf,[pathname3 name List{i} ext],'-dpng',['-r',num2str(EC.Img_sz(3))])
%         case '.TIF'
%             print(gcf,[pathname3 name List{i} ext],'-dtiff',['-r',num2str(EC.Img_sz(3))])
%         case '.JPG'
%             print(gcf,[pathname3 name List{i} ext],'-djpeg',['-r',num2str(EC.Img_sz(3))])
%         case '.BMP'
%             print(gcf,[pathname3 name List{i} ext],'-dbmp',['-r',num2str(EC.Img_sz(3))])
%         case '.EPS'
%             print(gcf,[pathname3 name List{i} ext],'-depsc',['-r',num2str(EC.Img_sz(3))])
%     end
% end
% msgbox('Image has saved!','Success','help');

% --------------------------------------------------------------------
function LoadFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to LoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
global a
global EC
H=BrainNet_LoadFiles;
uiwait(H);
if FLAG.LF==1
    EC.vol.display = [];
    EC.vol.pn = [];
    EC.vol.px = [];
    EC.vol.nn = [];
    EC.vol.nx = [];
    %     if FLAG.Loadfile~=1
    H=BrainNet_Option;
    uiwait(H);
    while FLAG.EC==2
        FLAG.EC=3;
        if ~isempty(a) && mean(ishandle(a))==1
            delete(a);
        end
        axes(handles.NV_axes);
        cla;
        a=FileView;
        uiwait(H);
    end
    if FLAG.EC==1
        FLAG.EC=0;
        if ~isempty(a) && mean(ishandle(a))==1
            delete(a);
        end
        axes(handles.NV_axes);
        cla;
        a=FileView;
    end
    FLAG.LF=0;
end

% --------------------------------------------------------------------
function SaveImage_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EC
[filename,pathname]=uiputfile({'*.tif','TIFF 24-bit';'*.jpg','JPEG 24-bit';...
    '*.bmp','BMP 24-bit';'*.eps','EPS color';...
    '*.png','PNG 24-bit';'*.fig','Matlab Figure';'*.wrl','VRML'},...
    'Save Image');
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    [pathstr, name, ext] = fileparts(fpath);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inch');
    set(gcf,'Paperposition',[1 1 EC.img.width/EC.img.dpi EC.img.height/EC.img.dpi]);
    switch ext
        case '.png'
            print(gcf,fpath,'-dpng',['-r',num2str(EC.img.dpi)])
        case '.tif'
            print(gcf,fpath,'-dtiff',['-r',num2str(EC.img.dpi)])
        case '.jpg'
            print(gcf,fpath,'-djpeg',['-r',num2str(EC.img.dpi)])
        case '.bmp'
            print(gcf,fpath,'-dbmp',['-r',num2str(EC.img.dpi)])
        case '.eps'
            print(gcf,fpath,'-depsc',['-r',num2str(EC.img.dpi)])
            
            % add by Mingrui Xia, save as matlab figure.
            % Modified by Mingrui, 20150112, using lowercase to be compatible with
            % higher version
        case '.fig'
            saveas(gcf,fpath,'fig');
            
            % added by Mingrui, 20170904, save as vrml format.
        case '.wrl'
            vrml(gcf,fpath);
    end
    msgbox('Image has saved!','Success','help');
end

% --------------------------------------------------------------------
function SagittalView_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SagittalView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam
global FLAG
global EC
if EC.lot.view == 1
    axis normal
    switch FLAG.sagittal
        case 0
            view(-90,0);
            FLAG.sagittal=1;
        case 1
            view(90,0);
            FLAG.sagittal=0;
    end
    if ~isempty(cam)
        for i = 1:length(cam) % Edited by Mingrui Xia, 20120809 adjust camlight in multi-surface view.
            camlight(cam(i));
        end
    end
    axis tight
    daspect([1 1 1])
end

% --------------------------------------------------------------------
function AxialView_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to AxialView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam
global FLAG
global EC
if EC.lot.view == 1
    axis normal
    switch FLAG.axial
        case 0
            view(0,90);
            FLAG.axial=1;
        case 1
            %             view(0,-90);
            view(-180,-90); % Edited by Mingrui Xia, 20181016, adjust ventral view.
            FLAG.axial=0;
    end
    if ~isempty(cam)
        for i = 1:length(cam) % Edited by Mingrui Xia, 20120809 adjust camlight in multi-surface view.
            camlight(cam(i));
        end
    end
    axis tight
    daspect([1 1 1])
end

% --------------------------------------------------------------------
function CoronalView_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to CoronalView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam
global FLAG
global EC
if EC.lot.view == 1
    axis normal
    switch FLAG.coronal
        case 0
            view(180,00);
            FLAG.coronal=1;
        case 1
            view(0,0);
            FLAG.coronal=0;
    end
    if ~isempty(cam)
        for i = 1:length(cam) % Edited by Mingrui Xia, 20120809 adjust camlight in multi-surface view.
            camlight(cam(i));
        end
    end
    axis tight
    daspect([1 1 1])
end

% --------------------------------------------------------------------
function Presentation_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Presentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam
global FLAG
FLAG.Rotate=1;%% Edited by Mingrui Xia 140523, fix the bug that demo and stop button are not matching
while FLAG.Rotate == 1
    camorbit(5,0,'camera');
    for i = 1:length(cam)
        camlight(cam(i));
    end
    drawnow;
end

% --------------------------------------------------------------------
function Stop_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
FLAG.Rotate=0;

% --------------------------------------------------------------------
function NV_m_movie_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam
[filename,pathname]=uiputfile({'*.avi','AVI movie'},'Save Movie');
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    vidObj = VideoWriter(fpath);
    open(vidObj);
    for num=1:360
        camorbit(1,0,'camera');
        camlight(cam);
        drawnow;
        print(gcf,[pathname,'temp.bmp'],'-dbmp');
        tempimg=imread([pathname,'temp.bmp']);
        currFrame=im2frame(tempimg);
        currFrame.cdata=imresize(currFrame.cdata,[534,735]);
        writeVideo(vidObj,currFrame);
    end
    close(vidObj);
    delete([pathname,'temp.bmp']);
    msgbox('Movie Saved!','Success','help');
end

% --------------------------------------------------------------------
function NV_m_so_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_so (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EC
[filename,pathname]=uiputfile({'*.mat','MAT-files (*.mat)'},'Save Option');
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    save(fpath,'EC');
    msgbox('Configure Saved!','Success','help');
end

% --------------------------------------------------------------------
function NV_m_option_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NV_m_lo_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.mat','MAT-files (*.mat)'},'Load Configuration');
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    
    % Edited by Mingrui Xia, 20140925, adapt New version to older configuration
    tmp = load(fpath);
    BrainNet_Option('CheckEC',tmp);
    msgbox('Option Loaded!','Success','help');
end

% --------------------------------------------------------------------
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printpreview;

% --------------------------------------------------------------------
function uitoggletool1_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam
for i = 1:length(cam)
    camlight(cam(i));
end

% --------------------------------------------------------------------
function NV_m_manual_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist('BrainNet_Manual.pdf','file')==2
    open('BrainNet_Manual.pdf');
else
    msgbox('Cannot find the manual file!','Error','error');
end

% --------------------------------------------------------------------
function NV_m_ColormapEditor_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_ColormapEditor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormapeditor;

% --------------------------------------------------------------------
function NV_m_ApplyColormap_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_ApplyColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EC
global a
EC.vol.CM = get(gca,'Colormap');
EC.vol.color_map = 22;
MapRange = get(a(end),'CLim');
for i = 1:length(a)
    axes(a(i));
    colormap(EC.vol.CM);
    caxis(MapRange);
end

% --------------------------------------------------------------------
function uitoggletool5_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function NV_fig_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to NV_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global EC
clear global a
clear global FLAG
clear global File
display('Thank you for using BrainNet Viewer!');

% --------------------------------------------------------------------
function NV_m_SaveColormap_Callback(hObject, eventdata, handles)
% hObject    handle to NV_m_SaveColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%% Added by Mingrui Xia, 20120710 save current colormap
Colormap = get(gca,'Colormap');
[filename,pathname]=uiputfile({'*.txt','Text Files'},'Save Colormap');
if isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    fid = fopen(fpath,'wt');
    fprintf(fid,'[');
    for i = 1:size(Colormap,1)-1
        fprintf(fid,'%f %f %f;',Colormap(i,:));
    end
    fprintf(fid,'%f %f %f]',Colormap(end,:));
    fclose(fid);
end

% --------------------------------------------------------------------
function MV_m_ViewMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to MV_m_ViewMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global surf
if ~isempty(surf.net)
    H = figure;
    imagesc(surf.net);
    daspect([1,1,1]);
    colorbar;
end

% --------------------------------------------------------------------
function NV_t_ViewMatrix_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to NV_t_ViewMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global surf
if ~isempty(surf.net)
    H = figure;
    imagesc(surf.net);
    daspect([1,1,1]);
    colorbar;
end


% --------------------------------------------------------------------
function NV_menu_saveboundary_Callback(hObject, eventdata, handles)
% hObject    handle to NV_menu_saveboundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global surf
global File
if isfield(surf,'boundary_value')
    [filename,pathname]=uiputfile({'*.txt','Text Files'},'Save Boundary');
    if isequal(filename,0)||isequal(pathname,0)
        return;
    else
        fpath=fullfile(pathname,filename);
        [~,meshname] = fileparts(File.MF);
        fid = fopen(fpath,'wt');
        fprintf(fid,'%s\n',['# for ',meshname]);
        for i = 1:length(surf.boundary_value)
            fprintf(fid,'%d\n',surf.boundary_value(i));
        end
        fclose(fid);
    end
else
end


% --------------------------------------------------------------------
function PlotEdit_toggle_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to PlotEdit_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function PlotEdit_toggle_OffCallback(hObject, eventdata, handles)
% hObject    handle to PlotEdit_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotedit('off');


% --------------------------------------------------------------------
function PlotEdit_toggle_OnCallback(hObject, eventdata, handles)
% hObject    handle to PlotEdit_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotedit('on');

% --------------------------------------------------------------------
function InitVariables
global FLAG
FLAG.Loadfile=0;
FLAG.EC=0;
FLAG.Rotate=0;
FLAG.EC_change=0;
FLAG.axial=0;
FLAG.sagittal=0;
FLAG.coronal=0;
FLAG.LF=0;
FLAG.MAP=0;
FLAG.IsCalledByCMD = 0;

global EC
EC=[];
EC.bak.color=[1,1,1];

EC.msh.color_type = 1; %Added by Mingrui, 20190521, load color for each vertex; 1 for same; 2 for defined

EC.msh.color=[0.95,0.95,0.95];
EC.msh.color_table = EC.msh.color;
EC.msh.color_table_tmp = EC.msh.color;

EC.msh.alpha=0.3;
EC.msh.doublebrain = 0; %%% Added by Mingrui Xia, 20120717, show two brains in one figure, 0 for one brain, 1 for two brains.

% Added by Mingrui, 20190327, parameters for boundary
EC.msh.boundary = 1;
% 1 for none, 2 for auto, 3 for threshold, and 4 for file
EC.msh.boundary_value = 0;
EC.msh.boundary_color = [0 0 0];

EC.nod.draw=1;
EC.nod.draw_threshold_type=1;
EC.nod.draw_threshold=0;
EC.nod.size=2;
EC.nod.size_size=3;
EC.nod.size_value=1;
EC.nod.size_threshold=0;
EC.nod.size_ratio=1;
EC.nod.color=1;
EC.nod.CM = [[37,89,152]/255;zeros(63,3)];

% Add by Mingrui, 20170317, custom colormap for node
EC.nod.CMt = EC.nod.CM;
EC.nod.cmstring = 'jet(64);';

EC.nod.color_map=1;

% Add by Mingrui, 20170303, fixed range of color mapping
EC.nod.color_map_type = 1; % 1 for auto, 2 for fix
EC.nod.color_map_low = 0;
EC.nod.color_map_high = 1;

EC.nod.color_threshold=0;
EC.nod.CMm=[0.956862745098039,1,0.545098039215686,0,0.247058823529412,0.913725490196078,1,0.803921568627451,0,0.129411764705882,0.611764705882353,1,1,0.298039215686275,0.0117647058823529,0.403921568627451,0.474509803921569,0.376470588235294,0.619607843137255,0,1;
    0.262745098039216,0.756862745098039,0.764705882352941,0.737254901960784,0.317647058823529,0.117647058823529,0.596078431372549,0.862745098039216,0.588235294117647,0.588235294117647,0.152941176470588,0.341176470588235,0.921568627450980,0.686274509803922,0.662745098039216,0.227450980392157,0.333333333333333,0.490196078431373,0.619607843137255,0,1;
    0.211764705882353,0.0274509803921569,0.290196078431373,0.831372549019608,0.709803921568628,0.388235294117647,0,0.223529411764706,0.533333333333333,0.952941176470588,0.690196078431373,0.133333333333333,0.231372549019608,0.313725490196078,0.956862745098039,0.717647058823529,0.282352941176471,0.545098039215686,0.619607843137255,0,1]';
EC.nod.ModularNumber = [];
% EC.nod.CMm = hsv(64);
% tmp = [];
% for i = 1:8
%     tmp = vertcat(tmp,EC.nod.CMm(i:8:64,:));
%     %     tmp = [tmp,EC.nod.CMm(i:8:256,:)'];
% end
% EC.nod.CMm = tmp;
% clear tmp

EC.lbl=2;
EC.lbl_threshold=0;
EC.lbl_threshold_type=1;

EC.edg.draw=1;
EC.edg.draw_threshold=0;
EC.edg.draw_abs=1;
EC.edg.size=2;
% EC.edg.size_size=0.3;
EC.edg.size_size = 1;
EC.edg.size_value=1;
EC.edg.size_threshold=0;
EC.edg.size_ratio=1;
EC.edg.size_abs=1;
EC.edg.color=1;
% 1 for same
% 2 for colormap
% 3 for threshold
% 4 for distance
% 5 for nodal module
% 6 for custom definition

% Add by Mingrui, 20170303, fixed range of color mapping
EC.edg.color_map_type = 1; % 1 for auto, 2 for fix
EC.edg.color_map_low = 0;
EC.edg.color_map_high = 1;

EC.edg.CM = [[1,0.9,0.73];ones(63,3)*0.38];

% Add by Mingrui 20170317, custom colormap for edge
EC.edg.CMt = EC.edg.CM;
EC.edg.cmstring = 'jet(64);';

EC.edg.color_map=1;
EC.edg.color_threshold=0;
EC.edg.color_distance=0;
EC.edg.color_abs=1;

% Add by Mingrui, 20150120, using custom matrix to define edge color.
EC.edg.CM_custom = [0.956862745098039,1,0.545098039215686,0,0.247058823529412,0.913725490196078,1,0.803921568627451,0,0.129411764705882,0.611764705882353,1,1,0.298039215686275,0.0117647058823529,0.403921568627451,0.474509803921569,0.376470588235294,0.619607843137255,0,1;
    0.262745098039216,0.756862745098039,0.764705882352941,0.737254901960784,0.317647058823529,0.117647058823529,0.596078431372549,0.862745098039216,0.588235294117647,0.588235294117647,0.152941176470588,0.341176470588235,0.921568627450980,0.686274509803922,0.662745098039216,0.227450980392157,0.333333333333333,0.490196078431373,0.619607843137255,0,1;
    0.211764705882353,0.0274509803921569,0.290196078431373,0.831372549019608,0.709803921568628,0.388235294117647,0,0.223529411764706,0.533333333333333,0.952941176470588,0.690196078431373,0.133333333333333,0.231372549019608,0.313725490196078,0.956862745098039,0.717647058823529,0.282352941176471,0.545098039215686,0.619607843137255,0,1]';
EC.edg.color_custom_matrix = [];
EC.edg.color_custom_index = [];



EC.edg.interhemiedges = 0; % Add by Mingrui Xia, 20120109, draw inter hemisphere edges.
EC.edg.directed = 0; % Add by Mingrui Xia, 20120621, draw directed network.

% Add by Mingrui, 20140925, dispaly edge by opacity
EC.edg.opacity = 1; % 1 for same opacity, 2 for value mapping
EC.edg.opacity_same = 0.7; % default 0.7 for same
EC.edg.opacity_max = 0.9;
EC.edg.opacity_min = 0.1;
EC.edg.opacity_abs = 1;

EC.img.width=2000;
EC.img.height=1500;
EC.img.dpi=300;

EC.lbl_font.FontName='Arial';
EC.lbl_font.FontWeight='bold';
EC.lbl_font.FontAngle='normal';
EC.lbl_font.FontSize=11;
EC.lbl_font.FontUnits='points';

% Added by Mingrui Xia, 20170919, add option to show label in front.
EC.lbl_front = 0;


EC.lot.view=1;
% 1 for single view
% 2 for full view
% 3 for lateral and medial view
% 4 for lateral, medial and ventral view
% 5 for lateral, medial and dorsal view

EC.lot.view_direction=1;
% Added by Mingrui Xia, 20120806, add custom view for single brain.
EC.lot.view_az = -90;
EC.lot.view_el = 0;


EC.vol.display = [];
EC.vol.pn = [];
EC.vol.px = [];
EC.vol.nn = [];
EC.vol.nx = [];
EC.vol.null=[0.95,0.95,0.95];
EC.vol.CM=ones(1000,3)*0.75;
EC.vol.color_map = 13;
EC.vol.cmstring = 'jet(1000);';
EC.vol.adjustCM = 1;

% Added by Mingrui Xia, 20160407, support annot format
EC.vol.CM_annot = rand(64,3);

% Added by Mingrui Xia, 20140916, add statistic for volume mapping
EC.vol.threshold = 0;
EC.vol.p = 0.05;
EC.vol.clustersize = 0;
EC.vol.rmm = 18;

% Added by Mingrui Xia, 20120726, selection for different mapping algorithm.
% 1 for Nearest Voxel
% 2 for Average Vertex
% 3 for Average Voxel
% 4 for Gaussian
% 5 for Interpolated (default)
% 6 for Maximum Voxel
% 7 for Minimum Voxel
% 8 for Extremum Voxel
% 9 for Most Neighbour Voxel, Added by Mingrui Xia, 20130104
EC.vol.mapalgorithm = 5;



EC.glb.material = 'dull';% Added by Mingrui Xia, 20120316, modify object material, shading, light.
EC.glb.material_ka = '0.5';
EC.glb.material_kd = '0.5';
EC.glb.material_ks = '0.5';
EC.glb.shading = 'interp';
EC.glb.lighting = 'phong';
EC.glb.lightdirection = 'right';
EC.glb.render = 'OpenGL'; % Added by Mingrui Xia, 20120413, selection for rendering methods
EC.glb.detail = 3; % Add by Mingrui Xia, 20120413, adjust graph detail

% Add by Mingrui, 20170309, display LR
EC.glb.lr = 1; % 1 for display, 0 for not display

% Added by Mingrui, 20120528, for ROI draw
EC.vol.type = 1; % 1 for volume to surface, 2 for ROI
EC.vol.roi.drawall = 1;
EC.vol.roi.draw = [];
EC.vol.roi.color = hsv(100);
EC.vol.roi.color = [EC.vol.roi.color(1:10:91,:)',EC.vol.roi.color(2:10:92,:)',EC.vol.roi.color(3:10:93,:)',EC.vol.roi.color(4:10:94,:)',EC.vol.roi.color(5:10:95,:)',EC.vol.roi.color(6:10:96,:)',EC.vol.roi.color(7:10:97,:)',EC.vol.roi.color(8:10:98,:)',EC.vol.roi.color(9:10:99,:)',EC.vol.roi.color(10:10:100,:)']';
EC.vol.roi.color = repmat(EC.vol.roi.color,11,1);
EC.vol.roi.colort = EC.vol.roi.color;
% EC.vol.roi.color = ones(2000,3)*0.7;
% EC.vol.roi.colort = EC.vol.roi.color;
EC.vol.roi.smooth = 1;
EC.vol.roi.drawcus = '';
EC.vol.roi.drawt = [];

% Added by Mingrui, 20170317, choice for smooth kernal
EC.vol.roi.smooth_kernal = 1;
% 1 for box
% 2 for Gaussian


global a
if ~isempty(a) && mean(ishandle(a))==1
    delete(a);
end
a=[];

global File
File.MF=[];
File.NI=[];
File.NT=[];
File.VF=[];