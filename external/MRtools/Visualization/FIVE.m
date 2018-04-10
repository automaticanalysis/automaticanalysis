function [Obj VOI peak] = FIVE(inputImage,linkHands)
%%% This is an image viewing and plotting utility.  Just launch it with
%%% "FIVE" at the command line.  There are a great many features in
%%% this viewer, and you can learn all about them at:
%%% http://nmr.mgh.harvard.edu/harvardagingbrain/People/AaronSchultz/Aarons_Scripts.html
%%%
%%% For launching FIVE you can do one of the following:
%%%     FIVE;  %% Launch FIVE and then load images
%%%     FIVE({'OverLay1.nii' 'Overlay2.img'});  %% This will launch FIVE with the specified overlays.
%%%     FIVE({{'Underlay.nii'}});  %% this will launch FIVE with the specified underlay image
%%%     FIVE({{'Underlay.nii'} {'OverLay1.nii' 'Overlay2.img'}});  %% This will launch FIVE with the specified underlay and load up the specified overlays.
%%%
%%% Ignore the linkHands input.  This is for specialized call back functions.
%%%
%%% Some FIVE features require additinal packages, for instance get
%%% peaks requires Donald McLaren's peak_nii package, and the VOI plotting
%%% options will require that the analysis be done with GLM_Flex.
%%%
%%% Written by Aaron P. Schultz - aschultz@martinos.org
%%%
%%% Copyright (C) 2011,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.
%%%
%%% system(['sh run_FIVE.sh /usr/pubsw/common/matlab/8.0/']);

% Update the colormaps, and colormapping system.

depth = 256;
[trash, cmaps] = colmap('jet',depth);
hand = [];
links = [];
origDat = [];
contrasts = [];
DM = [];
plotGo = 0;
Flist = [];
mniLims = [];
Des = [];
DataHeaders = [];
Outliers = [];

lastItem = [];
tabHand = [];
% global Obj;
Obj = [];

VOI = [];
peak = [];
modelType = [];

CachedClusterLoc = [];

ConExp = 0;
ConLayer = [];
ConHeader = [];

ssConExp = 0;
ssConLayer = [];
ssConHeader = [];
ssData = [];

clickFlag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the underlay
MH = spm_vol([fileparts(which('spm')) '/canonical/single_subj_T1.nii']); %% template underlay all image data is mapped to this orientation
if nargin > 0
    try
        if iscell(inputImage{1});
            [a b c] = fileparts(inputImage{1}{1});
            if strcmpi('.mgz',c)
                m = MRIread(inputImage{1}{1}); 
                m.descrip = []; 
                MRIwrite(m,[a b '.nii'],'float');
                inputImage{1}{1} = [a b '.nii'];
                %hh.fname = [m.fspec(1:end-3) 'nii'];
                %hh.dim = m.volsize;
                %hh.mat = m.vox2ras1;
                %hh.pinfo = [1;0;352];
                %hh.dt = [16 0];
                %hh.n = [1 1];
                %hh.descrip = ['MGZ Volume'];
                %hh.private = [];
            end
            
            h = spm_vol(inputImage{1}{1});
            [I mmm] = SliceAndDice3(h,MH,[],[],[3 NaN],[]);
            h.dim = size(I);
            h.mat = mmm;
            
        else
            [I h] = openIMG(which('defaultUnderlay.nii'));
        end
    catch
        [I h] = openIMG(which('defaultUnderlay.nii'));
    end
else
    [I h] = openIMG(which('defaultUnderlay.nii'));
end

Obj = initializeUnderlay(I,h);
movego = 0;
loc = [0 0 0];
mc = round([loc 1] * inv(Obj(1).h.mat)');
Obj(1).point =  round([loc 1] * inv(Obj(1).h.mat)');
Obj(1).lastpoint = [round(size(Obj(1).I)/2) 1]-1;
Obj(1).Range = [min(Obj(1).I(:)) max(Obj(1).I(:)) max(Obj(1).I(:))-min(Obj(1).I(:))];

Obj(1).FOV = sort([ [1 1 1 1]*Obj(1).h.mat'; [Obj(1).axLims 1]*Obj(1).h.mat']);
Obj(1).col = 1;

if any((size(Obj(1).I) - Obj(1).point(1:3))<0); 
    Obj(1).point = [round(size(Obj(1).I)/2) 1];
    tmp = Obj(1).point*Obj(1).h.mat';
    loc = tmp(1:3);
end
% if ~isempty(searchCellStr('aschultz',{UserTime})); keyboard; end

setupFrames(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the figure window
height = []; width = []; rat = [];
ax1 = []; ax2 = []; ax3 = []; ax4 = [];
con = []; menu = []; hcmenu = []; item = [];
setupFigure;
Obj(1).con = con;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the figure menus
paramenu1 = []; paramenu2 = []; paramenu3 = []; paramenu4 = []; S = [];
setupParamMenu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in a label atlas
RNH = spm_vol([which('aal_MNI_V4.img')]);
[RNI Rxyz] = spm_read_vols(RNH);
RNames = load('aal_MNI_V4_List.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw Underlay
ch = [];
drawFresh(ax1,1);
drawFresh(ax2,2);
drawFresh(ax3,3);
Obj(1).hand = hand;
for ii = 1:length(hand);
    set(hand{ii}, 'uicontextmenu',hcmenu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot CrossHairs
[ch(1,1), ch(1,2)] = crossHairs(ax1,[0 0]);
[ch(2,1), ch(2,2)] = crossHairs(ax2,[0 0]);
[ch(3,1), ch(3,2)] = crossHairs(ax3,[0 0]);
Obj(1).ch = ch;
set(ch, 'uicontextmenu',hcmenu);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Link CrossHair Positions
hListener1 = addlistener(ch(1,1),'XData','PostSet',@AutoUpdate);
hListener2 = addlistener(ch(2,1),'XData','PostSet',@AutoUpdate);
hListener3 = addlistener(ch(1,2),'YData','PostSet',@AutoUpdate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synchronize Views if requested
if nargin>1 && ~isempty(linkHands)
    links{end+1}=linkprop([ch(1,1),linkHands(1,1)],'XData');
    links{end+1}=linkprop([ch(1,2),linkHands(1,2)],'YData');
    links{end+1}=linkprop([ch(2,1),linkHands(2,1)],'XData');
    links{end+1}=linkprop([ch(2,2),linkHands(2,2)],'YData');
    links{end+1}=linkprop([ch(3,1),linkHands(3,1)],'XData');
    links{end+1}=linkprop([ch(3,2),linkHands(3,2)],'YData');
    set(gcf,'UserData',links);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open Prespecified Overlays
Count = 2;
if nargin > 0;    
%     if isa(inputImage,'uint8')
%         [stuff ras] = ReadFileFromBlob(bin);
%     else
        try
            if iscell(inputImage{1})
                if length(inputImage)>1
                    openOverlay(inputImage{2});
                end
            else
                openOverlay(inputImage);
            end
        catch
            openOverlay(inputImage);
        end
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yy = Obj(1).clim(1):(Obj(1).clim(2)-Obj(1).clim(1))/255:Obj(1).clim(2);
axes(ax4); cla
imagesc(yy,1,reshape(colmap(cmaps{Obj(1).col},256),256,1,3));
axis tight;
set(ax4,'YDir','Normal','YAxisLocation','right','YTick',  unique([1 get(ax4,'YTick')]));
set(ax4,'YTickLabel',round(min(yy):spm_range(yy)/(numel(get(ax4,'YTick'))-1):max(yy)));
set(ax4,'fontsize',6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Figure Visible
set(findobj(gcf,'Fontsize', 10),'fontsize',12); shg
set(findobj(gcf,'Fontsize', 12),'fontsize',13); shg
set(pane,'Visible','on');

%%
    function AutoUpdate(varargin)
        
        vn = get(con(21,1),'Value');
        
        y = get(ch(1,1),'XData');
        x = get(ch(2,1),'XData');
        z = get(ch(1,2),'YData');
        xyz = [x(1) y(1) z(1)];
        
        xyz((mniLims(1,:)-xyz)>0) = mniLims(1,(mniLims(1,:)-xyz)>0);
        xyz((mniLims(2,:)-xyz)<0) = mniLims(2,(mniLims(2,:)-xyz)<0);
        x = xyz(1); y = xyz(2); z = xyz(3);
        
        for ii = 1:length(Obj)
            Obj(ii).point =  ceil([x y z 1] * inv(Obj(ii).h.mat)');
            
            Obj(ii).point(Obj(ii).point(1:3)<1) = 1;
            ind = find((Obj(ii).axLims-Obj(ii).point(1:3))<0);
            Obj(ii).point(ind) = Obj(ii).axLims(ind);
            
            if ii == 1; 
                p = ceil([x y z 1] * inv(RNH.mat)');
               try
                   nm = RNames.ROI(RNI(p(1),p(2),p(3)));
                   set(paramenu1(2),'Label',nm.Nom_L); 
               catch
                   set(paramenu1(2),'Label','undefined');
               end
            end
        end
        
        setupFrames(1:length(Obj),0);
        updateGraphics([1 2 3],0);
        set(con(15,1),'String',num2str(round([x y z])));
        mc = round([x y z 1] * inv(Obj(get(con(21,1),'Value')).h.mat)');
        
        vn = get(con(21,1),'Value');
        
        try
            set(con(22,1),'String',num2str(Obj(get(con(21,1),'Value')).I(mc(1),mc(2),mc(3))));
        catch
            set(con(22,1),'String','NaN');
        end
        
        
        shg
    end

    function goTo(varargin)
        if varargin{1}==con(15,1);
            cl = str2num(get(con(15,1),'String'));
            set(ch(1,1), 'XData', [cl(2) cl(2)]);
            set(ch(1,2), 'YData', [cl(3) cl(3)]);
            set(ch(2,1), 'XData', [cl(1) cl(1)]);
            set(ch(2,2), 'YData', [cl(3) cl(3)]);
            set(ch(3,1), 'XData', [cl(1) cl(1)]);
            set(ch(3,2), 'YData', [cl(2) cl(2)]);
        end
    end

    function buttonUp(varargin)
        movego = 0;
    end

    function buttonDown(varargin)
        %get(gcf,'SelectionType')
        
        if strcmpi(get(gcf,'SelectionType'),'normal') || clickFlag==1;
            movego = 1;
            
            co1 = gco;
            co2 = gca;
            switch co2
                case ax1
                    cp = get(co2, 'CurrentPoint');
                    cp = round(cp(1,1:2));
                    y = cp(1);
                    z = cp(2);
                    set(ch(1,1), 'XData', [y y]);
                    set(ch(1,2), 'YData', [z z]);
                    set(ch(2,2), 'YData', [z z]);
                    set(ch(3,2), 'YData', [y y]);
                case ax2
                    cp = get(co2, 'CurrentPoint');
                    cp = round(cp(1,1:2));
                    x = cp(1);
                    z = cp(2);
                    
                    set(ch(2,1), 'XData', [x x]);
                    set(ch(2,2), 'YData', [z z]);
                    set(ch(1,2), 'YData', [z z]);
                    set(ch(3,1), 'XData', [x x]);
                case ax3
                    cp = get(co2, 'CurrentPoint');
                    cp = round(cp(1,1:2));
                    y = cp(1);
                    x = cp(2);
                    
                    set(ch(3,1), 'XData', [y y]);
                    set(ch(3,2), 'YData', [x x]);
                    set(ch(1,1), 'XData', [x x]);
                    set(ch(2,1), 'XData', [y y]);
                otherwise
            end
        end
        
        if strcmpi(get(gcf,'SelectionType'),'extend');
            plotVOI(lastItem);
        end
        
        if strcmpi(get(gcf,'SelectionType'),'alt');
            if ConExp ~=0 
                updateConnMap;
            elseif ssConExp ~= 0;
                ssUpdateConnMap;
            end
        end
            
    end

    function buttonMotion(varargin)
        if movego == 1;
            co1 = gco;
            
            co2 = gca;
            
            switch co2
                case ax1
                    
                    cp = get(co2, 'CurrentPoint');
                    cp = round(cp(1,1:2));
                    y = cp(1);
                    z = cp(2);
                    
                    if Obj(1).FOV(1,2)>y; y=Obj(1).FOV(1,2); end;
                    if Obj(1).FOV(2,2)<y; y=Obj(1).FOV(2,2); end;
                    if Obj(1).FOV(1,3)>z; z=Obj(1).FOV(1,3); end;
                    if Obj(1).FOV(2,3)<z; z=Obj(1).FOV(2,3); end;
                    
                    set(ch(1,1), 'XData', [y y]);
                    set(ch(1,2), 'YData', [z z]);
                    set(ch(2,2), 'YData', [z z]);
                    set(ch(3,2), 'YData', [y y]);

                case ax2
                    cp = get(co2, 'CurrentPoint');
                    cp = round(cp(1,1:2));
                    x = cp(1);
                    z = cp(2);
                    
                    if Obj(1).FOV(1,1)>x; x=Obj(1).FOV(1,1); end;
                    if Obj(1).FOV(2,1)<x; x=Obj(1).FOV(2,1); end;
                    if Obj(1).FOV(1,3)>z; z=Obj(1).FOV(1,3); end;
                    if Obj(1).FOV(2,3)<z; z=Obj(1).FOV(2,3); end;
                    
                    
                    set(ch(2,1), 'XData', [x x]);
                    set(ch(2,2), 'YData', [z z]);
                    set(ch(1,2), 'YData', [z z]);
                    set(ch(3,1), 'XData', [x x]);

                case ax3
                    cp = get(co2, 'CurrentPoint');
                    cp = round(cp(1,1:2));
                    y = cp(1);
                    x = cp(2);
                    
                    if Obj(1).FOV(1,1)>x; x=Obj(1).FOV(1,1); end;
                    if Obj(1).FOV(2,1)<x; x=Obj(1).FOV(2,1); end;
                    if Obj(1).FOV(1,2)>y; y=Obj(1).FOV(1,2); end;
                    if Obj(1).FOV(2,2)<y; y=Obj(1).FOV(2,2); end;
                    
                    set(ch(3,1), 'XData', [y y]);
                    set(ch(3,2), 'YData', [x x]);
                    set(ch(1,1), 'XData', [x x]);
                    set(ch(2,1), 'XData', [y y]);

                otherwise
            end
        end
    end

    function openOverlay(varargin)
        if nargin == 0
            nnn = spm_select(inf,'image');
        else
            if iscell(varargin{1}) 
                if ~isa(varargin{1}{1},'uint8')
                    nnn = char(varargin{1});
                else
                    nnn = {};
                    for ii = 1:numel(varargin{1})
                        [mm hh] = ReadFileFromBlob(varargin{1}{ii});
                        hh.n = [ii 1];
                        spm_write_vol(hh,mm);
                        nnn{ii,1} = [hh.fname ',' num2str(ii)];
                    end
                    nnn = char(nnn);
                end
            elseif ischar(varargin{1})
                nnn = varargin{1};
            elseif isa(varargin{1},'uint8')
                [mm hh] = ReadFileFromBlob(varargin{1});
                spm_write_vol(hh,mm);
                nnn = hh.fname;
            else
                nnn = spm_select(inf,'image');
            end
            
            for kk = 1:size(nnn,1)
%                 if isa(nnn,'uint8')
%                     fid = fopen([pwd '/tmp.nii']);
%                     fwrite(fid,nnn,'uint8');
%                     fclose(fid);
%                     n = [pwd '/tmp.nii'];
% %                 elseif isa(nnn{kk},'uint8')
% %                     fid = fopen(['/tmp' '/tmp.nii'],'w');
% %                     fwrite(fid,nnn{kk},'uint8');
% %                     fclose(fid);
% %                     n = ['/tmp' '/tmp.nii'];
%                 else
%                     n = strtrim(nnn(kk,:));
%                 end
                
                n = strtrim(nnn(kk,:));
                
                
                ind = find(n==filesep);
                if isempty(ind);
                    nn = n;
                    n2 = n;
                else
                    nn = n(ind(end)+1:end);
                    if numel(ind)>1
                        n2 = n(ind(end-1)+1:end);
                    else
                        n2 = n;
                    end
                end
                if mean(n==filesep)==0
                    n = [pwd filesep n];
                end
                
                hh = spm_vol(n);
                
                set(con(21,1),'TooltipString',n)
                
                a = world_bb(MH); b = world_bb(hh); tmp = a-b; tmp(1,:) = tmp(1,:)*-1;
                
                if hh.dt(1)>=16
                    intOrd = 3;
                else
                    intOrd = 0;
                end
                
                if all(tmp(:)>=0)
                    [m mmm] = SliceAndDice3(hh,      MH,    [],    hh,[intOrd NaN],[]);
                else                    
                    [m mmm] = SliceAndDice3(hh,      MH,    [],    Obj(1).h,[intOrd NaN],[]);
                end

                %[m mmm] = SliceAndDice3(hh,MH,[],Obj(1).h,[0 NaN],[]);
                
                hh.dim = size(m);
                hh.mat = mmm;
                
                Obj(Count).Name = n2;
                Obj(Count).FullPath = n;
                Obj(Count).DispName = n2;
                set(gcf,'Name', ['FIVE: ' Obj(Count).DispName]);
                Obj(Count).h = hh;
                Obj(Count).I = double(m);

                set(con(21,1),'String', [get(con(21,1),'String'); n2],'Value',Count);
                
                if ~isempty(searchCellStr('{T_', {hh.descrip}));
                    tmp = hh.descrip;
                    i1 = find(tmp=='[');
                    i2 = find(tmp==']');
                    
                    Obj(Count).DF = str2num(tmp(i1+1:i2-1));
                    th = spm_invTcdf(1-.001,Obj(Count).DF);
                    Obj(Count).Thresh = [th ceil(max(m(:))*1000)/1000];
                    Obj(Count).clim = [floor(min(m(:))*1000)/1000 ceil(max(m(:))*1000)/1000];
                    Obj(Count).PVal = .001;
                    Obj(Count).col = 2;
                    Obj(Count).Trans = 1;
                    
                    
                    set(con(4),'String',[num2str(th) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
                    set(con(5),'String',[num2str(floor(min(m(:))*1000)/1000) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
                    
                    set(con(8),'String',num2str(Obj(Count).DF));
                    set(con(9),'String','++0.001');
                elseif ~isempty(searchCellStr('{F_', {hh.descrip}))
                    tmp = hh.descrip;
                    i1 = find(tmp=='[');
                    i2 = find(tmp==']');
                    
                    t1 = regexp(tmp(i1(1)+1:i2(1)-1),',','split');
                    df(1) = str2num(t1{1});
                    df(2) = str2num(t1{2});
                    
                    Obj(Count).DF = df;
                    th = spm_invFcdf(1-.001,df(1), df(2));
                    Obj(Count).Thresh = [th ceil(max(m(:))*1000)/1000];
                    Obj(Count).clim = [floor(min(m(:))*1000)/1000 ceil(max(m(:))*1000)/1000];
                    Obj(Count).PVal = .001;
                    Obj(Count).col = 2;
                    Obj(Count).Trans = 1;
                    
                    set(con(4),'String',[num2str(th) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
                    set(con(5),'String',[num2str(floor(min(m(:))*1000)/1000) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
                    
                    set(con(8),'String',num2str(Obj(Count).DF));
                    set(con(9),'String','++0.001');
                else
                    set(con(4),'String',[num2str(floor(min(m(:))*1000)/1000) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
                    set(con(5),'String',[num2str(floor(min(m(:))*1000)/1000) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
                    set(con(8),'String','NaN');
                    set(con(9),'String','NaN');
                    
                    Obj(Count).DF = NaN;
                    Obj(Count).PVal = NaN;
                    Obj(Count).Thresh = [floor(min(m(:))*1000)/1000 ceil(max(m(:))*1000)/1000];
                    Obj(Count).clim = [floor(min(m(:))*1000)/1000 ceil(max(m(:))*1000)/1000];
                    Obj(Count).col = 2;
                    Obj(Count).Trans = 1;
                end
                
                set(con(23),'String','0');
                Obj(Count).ClusterThresh = 0;
                Obj(Count).Exclude = [];
                Obj(Count).MaskInd = [];
                Obj(Count).mask = ones(size(Obj(Count).I),'uint8');
               
                                
                y = get(ch(1,1),'XData');
                x = get(ch(2,1),'XData');
                z = get(ch(1,2),'YData');
                loc = [x(1) y(1) z(1)];
                
                Obj(Count).point = round([loc 1] * inv(Obj(Count).h.mat)');
                Obj(Count).lastpoint = (round([loc 1] * inv(Obj(Count).h.mat)'))-1;
                Obj(Count).axLims = size(Obj(Count).I);
                
                setupFrames(Count,0);
                Obj(Count).pos = axLim(Obj(Count).I,Obj(Count).h);
                                
                set(con(1,1),'Value',3);
                
                drawFresh(ax1,1,Count);
                drawFresh(ax2,2,Count);
                drawFresh(ax3,3,Count);
                
                uistack(ch(1,1),'top'); uistack(ch(1,2),'top');
                uistack(ch(2,1),'top'); uistack(ch(2,2),'top');
                uistack(ch(3,1),'top'); uistack(ch(3,2),'top');
                
                Obj(Count).Range = [min(Obj(Count).I(:)) max(Obj(Count).I(:)) max(Obj(Count).I(:))-min(Obj(Count).I(:))];
                
                for ii = 1:length(hand);
                    set(hand{ii}, 'uicontextmenu',hcmenu);
                end
                

                yy = Obj(Count).clim(1):(Obj(Count).clim(2)-Obj(Count).clim(1))/255:Obj(Count).clim(2);
                axes(ax4); cla

                try
                    imagesc(yy,1,reshape(colmap(cmaps{Obj(Count).col},256),256,1,3));
                    axis tight;
                    set(ax4,'YDir','Normal','YAxisLocation','right','YTick',  unique([1 get(ax4,'YTick')]));
                    set(ax4,'YTickLabel',round(min(yy):spm_range(yy)/(numel(get(ax4,'YTick'))-1):max(yy)));
                    set(ax4,'fontsize',6);
                end
                Count = Count+1;
            end
        end
    end

    function UpdateThreshold(varargin)
        % Get the Volume Number.
        vn = get(con(21,1),'Value');
        
        % Get the ThresholdParamters
        a = get(con(4),'String');
        c = get(con(8),'String');
%         d = get(con(9),'String');
        
        % Parse the Threshold Parameters
        if isempty(a);
            a = [-inf inf];
        else
            if searchCellStr(',',{a})
                ind = find(a==',');
                a = [str2num(strtrim(a(1:ind-1))) str2num(strtrim(a(ind+1:end)))];
            else
                a = str2num(a);
            end
        end
        
        
        if numel(a) == 2
            % Correct the Threshold Parameter End Points
            if  a(1)==-inf
                a(1) = floor(min(Obj(vn).I(:))*1000)/1000;
                set(con(4),'string',[num2str(a(1)) ', ' num2str(a(2))])
            end
            if a(2)==inf
                a(2) = ceil(max(Obj(vn).I(:))*1000)/1000;
                set(con(4),'string',[num2str(a(1)) ', ' num2str(a(2))])
            end
            
            if a(1)>0 && a(2)>0
                if numel(varargin)~=0 && varargin{1} ~= con(9,1)
                    if ~strcmpi(get(con(8,1),'String'),'NA')
                        df =  str2num(get(con(8,1),'String'));
                        if isnan(df)
                            p = NaN;
                        else
                            if numel(df) == 1;
                                if isinf(df)
                                    p = 1-spm_Ncdf(a(1),0,1);
                                else
                                    p = 1-spm_Tcdf(a(1),df);
                                end
                                set(con(9,1),'String',['++' num2str(p)]);
                            end
                            if numel(df) == 2;
                                p = 1-spm_Fcdf(a(1),df(1),df(2));
                                set(con(9,1),'String',['++' num2str(p)]);
                            end
                        end
                    end
                end
            end
            
            
            if a(1)<0 && a(2)<0
                if numel(varargin)~=0 && varargin{1} ~= con(9,1)
                    if ~strcmpi(get(con(8,1),'String'),'NA')
                        df =  str2num(get(con(8,1),'String'));
                        if isnan(df)
                            p = NaN;
                        else
                            if numel(df) == 1;
                                if isinf(df)
                                    p = 1-spm_Ncdf(abs(a(2)),0,1);
                                else
                                    p = 1-spm_Tcdf(abs(a(2)),df);
                                end
                                set(con(9,1),'String',['--' num2str(p)]);
                            end
                            if numel(df) == 2;
                                p = 1-spm_Fcdf(abs(a(2)),df(1),df(2));
                                set(con(9,1),'String',['--' num2str(p)]);
                            end
                        end
                    end
                end
            end
            
            Obj(vn).Thresh = a;
            
        elseif numel(a) == 4;
            if  a(1)==-inf
                a(1) = floor(min(Obj(vn).I(:))*1000)/1000;
                set(con(4),'String',[num2str(a(1)) ' ' num2str(a(2)) ', ' num2str(a(3)) ' ' num2str(a(4))]);
            end
            if a(4)==inf
                a(4) = ceil(max(Obj(vn).I(:))*1000)/1000;
                set(con(4),'String',[num2str(a(1)) ' ' num2str(a(2)) ', ' num2str(a(3)) ' ' num2str(a(4))]);
            end

            
            if numel(Obj(vn).Thresh)==numel(a);
                check = Obj(vn).Thresh-a;
                check = find(check~=0);
                check = setdiff(check,[1 4]);
                
                if numel(check)==1;
                    other = setdiff(2:3,check);
                    a(other) = -a(check);
                    set(con(4),'String',[num2str(a(1)) ' ' num2str(a(2)) ', ' num2str(a(3)) ' ' num2str(a(4))]);
                end

                if numel(varargin)>0 && varargin{1} ~= con(9,1)
                    if ~(strcmpi(get(con(8,1),'String'),'NA') || ~isempty(searchCellStr('nan', {lower(get(con(8,1),'String'))})))
                        df =  str2num(get(con(8,1),'String'));
                        if numel(df) == 1;
                            if isinf(df)
                                p = (1-spm_Ncdf(abs(a(3)),0,1))*2;
                            else
                                p = (1-spm_Tcdf(abs(a(3)),df))*2; 
                            end
                            set(con(9,1),'String',['+-' num2str(p)]);
                        end
                        if numel(df) == 2;
                            p = 1-spm_Fcdf(abs(a(3)),df(1),df(2));
                            set(con(9,1),'String',['+-' num2str(p)]);
                        end
                    end
                end
            else
                if abs(a(2))~=abs(a(3))
                    ind1 = find(abs(a(2:3))==min(abs(a(2:3))));
                    ind2 = find(abs(a(2:3))==max(abs(a(2:3))));
                    a(ind1+1) = -a(ind2+1);
                    set(con(4),'String',[num2str(a(1)) ' ' num2str(a(2)) ', ' num2str(a(3)) ' ' num2str(a(4))]);
                    tmp = get(con(8,1),'String');
                end
                

                
                if numel(varargin)>0 && varargin{1} ~= con(9,1)
                    if ~strcmpi(get(con(8,1),'String'),'nan')
                        df =  str2num(get(con(8,1),'String'));
                        if numel(df) == 1;
                            if isinf(df)
                                p = 1-spm_Ncdf(abs(a(2))/2,0,1);
                            else
                                p = 1-spm_Tcdf(abs(a(2))/2,df);
                            end
                            set(con(9,1),'String',['+-' num2str(p)]);
                        end
                        if numel(df) == 2;
                            p = 1-spm_Fcdf(abs(a(2)),df(1),df(2));
                            set(con(9,1),'String',['+-' num2str(p)]);
                        end
                    end
                end
            end
            Obj(vn).Thresh = a;
        end
                
        loc = str2num(get(con(15),'String'));
        Obj(vn).point =  round([loc 1] * inv(Obj(vn).h.mat)');
        
        
        if str2num(get(con(23),'String'))~=0
            ExtentThresh;
        end
        
        setupFrames(vn,1);
        updateGraphics({1:3 vn},1);
        %Obj(vn).pos = axLim(Obj(vn).I,Obj(vn).h);
    end

    function UpdatePVal(varargin)
        nv = get(con(21,1),'Value');
        
        df = get(con(8,1),'String');
        try
            df = str2num(df);
            if any(isnan(df));
                return;
            end
        catch
            set(con(8,1),'String','NA')
            return
        end
        
        p = get(con(9,1),'String');
        
        if mean(p(1:2)=='++')==1
            direc = 1;
            p = str2num(p(3:end));
        elseif mean(p(1:2)=='--')==1
            direc = -1;
            p = str2num(p(3:end));
        elseif mean(p(1:2)=='-+')==1 || mean(p(1:2)=='+-')==1
            direc = 0;
            p = str2num(p(3:end))/2;
        else
            p = str2num(p);
            direc = 1;
        end
        
        Obj(nv).PVal = p;
        Obj(nv).DF = df;
        
        if numel(df)==1
            if isinf(df)
                T = spm_invNcdf(1-abs(p),0,1);
            else
                T = spm_invTcdf(1-abs(p),df);
            end
        elseif numel(df)==2
            T = spm_invFcdf(1-abs(p),df(1),df(2));
        end
        
        curr = get(con(4),'String');
        cur = regexp(curr,',','split');
        
        if     direc ==  1
            a = ceil(T*1000)/1000;
            set(con(4),'String',[num2str(a) ', inf']);
            %set(con(4),'String',[num2str(a) ',' cur{2}]);
            
            UpdateThreshold;
        elseif direc == -1
            a = floor(-T*1000)/1000;
            set(con(4),'String',['-inf ,' num2str(a)]);
            %set(con(4),'String',[cur{1} ',' num2str(a)]);
            
            UpdateThreshold;
        elseif direc ==  0
            a = [];
            %keyboard;
            a(1) = floor(min(Obj(nv).I(:))*1000)/1000;
            a(4) = ceil(max(Obj(nv).I(:))*1000)/1000;
            a(2:3) = [floor(-T*1000)/1000 ceil(T*1000)/1000];
            set(con(4),'String',[num2str(a(1)) ' ' num2str(a(2)) ', ' num2str(a(3)) ' ' num2str(a(4))]);
            b = a([1 4]);
            ind1 = find(abs(b)==min(abs(b)));
            ind2 = find(abs(b)==max(abs(b)));
            b(ind2) = -b(ind1);

            UpdateThreshold;
        end
    end

    function changeColorMap(varargin)
        vn = get(con(21,1),'Value');
        map = get(con(1,1),'Value');
        if map == 1;
            set(con(1,1),'Value',Obj(vn).col+1);
            return;
        end
        Obj(vn).col = map-1;
        
        setupFrames(vn,1);
        updateGraphics({[1 2 3] vn},1);
        
        %%% Update colorbar.
        yy = Obj(vn).clim(1):(Obj(vn).clim(2)-Obj(vn).clim(1))/255:Obj(vn).clim(2);
        axes(ax4); cla
        imagesc(yy,1,reshape(colmap(cmaps{Obj(vn).col},256),256,1,3));
        axis tight;
        set(ax4,'YDir','Normal','YAxisLocation','right','YTick',  unique([1 get(ax4,'YTick')]));
        set(ax4,'YTickLabel',round(min(yy):spm_range(yy)/(numel(get(ax4,'YTick'))-1):max(yy)));      
        set(ax4,'fontsize',6);

    end

    function adjustTrans(varargin)
        vn = get(con(21,1),'Value');
        if numel(hand{1})==1
            return
        end
        const = get(varargin{1},'Value');
        const = const+.0001;
        if const>1; const = 1; end;
        
        Obj(vn).Trans = const;
        
        a = get(hand{1}(vn),'AlphaData');
        a(a~=0)=const;
        set(hand{1}(vn), 'AlphaData', a);
        
        a = get(hand{2}(vn),'AlphaData');
        a(a~=0)=const;
        set(hand{2}(vn), 'AlphaData', a);
        
        a = get(hand{3}(vn),'AlphaData');
        a(a~=0)=const;
        set(hand{3}(vn), 'AlphaData', a);
    end

    function resizeFig(varargin)
        a = get(gcf,'Position');
        
        hei = (height/width)*a(3);
        wid = (width/height)*a(4);
        
        if hei-a(4) > wid-a(3)
            a(3) = wid;
        elseif hei-a(4) < wid-a(3)
            a(4) = hei;
        end
        
        set(gcf,'Position',a);
    end

    function switchObj(varargin)
        if ~isempty(varargin)
            if iscell(varargin{1})
                nv = varargin{1}{1};    
            else
                nv = get(con(21,1),'Value');
            end
        else
            nv = get(con(21,1),'Value');
        end

        
        try
            set(con(21,1),'TooltipString',Obj(nv).FullPath)
            set(gcf,'Name', ['FIVE: ' Obj(nv).DispName]);
        catch
        end
        
        set(con(2,1),'Value', Obj(nv).Trans);
        
        set(con(23),'String', num2str(Obj(nv).ClusterThresh));
        
        set(con(5,1),'String',[num2str(Obj(nv).clim(1)) ', ' num2str(Obj(nv).clim(2))]);
        set(con(8,1), 'String',num2str(Obj(nv).DF));
        
        if numel(Obj(nv).Thresh) == 4;
            set(con(4,1),'String',[num2str(Obj(nv).Thresh(1)) ' ' num2str(Obj(nv).Thresh(2)) ', ' num2str(Obj(nv).Thresh(3)) ' ' num2str(Obj(nv).Thresh(4))]);
            set(con(9,1), 'String',['+-' num2str(Obj(nv).PVal*2)]);
        elseif sum(sign(Obj(nv).Thresh))>0
            set(con(4,1),'String',[num2str(Obj(nv).Thresh(1)) ', ' num2str(Obj(nv).Thresh(2))]);
            set(con(9,1), 'String',['++' num2str(Obj(nv).PVal)]);
        elseif sum(sign(Obj(nv).Thresh))<0
            set(con(4,1),'String',[num2str(Obj(nv).Thresh(1)) ', ' num2str(Obj(nv).Thresh(2))]);
            set(con(9,1), 'String',['--' num2str(Obj(nv).PVal)]);
        elseif sum(sign(Obj(nv).Thresh))==0
            set(con(4,1),'String',[num2str(Obj(nv).Thresh(1)) ', ' num2str(Obj(nv).Thresh(2))]);
            set(con(9,1), 'String',num2str(Obj(nv).PVal));
        end
        set(con(1,1),'Value', Obj(nv).col+1);
        
        yy = Obj(nv).clim(1):(Obj(nv).clim(2)-Obj(nv).clim(1))/255:Obj(nv).clim(2);
        axes(ax4); cla
        imagesc(yy,1,reshape(colmap(cmaps{Obj(nv).col},256),256,1,3));
        axis tight;
        set(ax4,'YDir','Normal','YAxisLocation','right','YTick',  unique([1 get(ax4,'YTick')]));
        set(ax4,'YTickLabel',round(min(yy):spm_range(yy)/(numel(get(ax4,'YTick'))-1):max(yy)));
        set(ax4,'fontsize',6);

        
        
        tmp1 = (Obj(nv).clim(2)-Obj(nv).Range(1))/Obj(nv).Range(3); if tmp1>1; tmp1=1; end; if tmp1<0; tmp1=0; end;
        set(con(27),'Value', tmp1);
        tmp1 = (Obj(nv).clim(1)-Obj(nv).Range(1))/Obj(nv).Range(3); if tmp1>1; tmp1=1; end; if tmp1<0; tmp1=0; end;
        set(con(28),'Value', tmp1);
    end

    function removeVolume(varargin)
        nv = get(con(21,1),'Value');
        if nv == 1; return; end
        ind = setdiff(1:length(Obj),nv);
        tmp = get(con(21,1),'String');
        set(con(21,1),'String',tmp(ind), 'Value',length(ind));
        
        for ii = 1:length(hand)
            delete(hand{ii}(nv));
        end
        
        Obj = Obj(ind);
        hand{1} = hand{1}(ind);
        hand{2} = hand{2}(ind);
        hand{3} = hand{3}(ind);
        
        Count = length(ind)+1;
        switchObj;
    end

    function updateGraphics(HandInd,opt)
        if iscell(HandInd);
            HandInd2 = HandInd{2};
            HandInd = HandInd{1};
        else
            HandInd2 = 1:length(Obj);
        end
        
        
        for ii = HandInd;
            if Obj(1).point(ii)==Obj(1).lastpoint(ii) && opt~=1
                continue
            end
            for jj = HandInd2
                set(hand{ii}(jj),'CData',Obj(jj).frame{ii},'AlphaData',~isnan(Obj(jj).frame{ii}(:,:,1))*Obj(jj).Trans);
                Obj(jj).lastpoint(ii) = Obj(jj).point(ii);
            end
        end
        
    end

    function drawFresh(axx,opt,opt2,opt3)
        if nargin == 2;
            opt2 = 1:length(Obj);
        end
        if nargin < 4
            opt3 = 1;
        end
        axes(axx);
        
        for ii = opt2;
            if opt == 1;
                tmp = image(Obj(ii).pos{2}, Obj(ii).pos{3}, Obj(ii).frame{opt});
                if opt3; hand{opt}(ii) = tmp; end
                set(tmp,'AlphaData', ~isnan(Obj(ii).frame{opt}(:,:,1))*Obj(ii).Trans);
                set(axx, 'YDir','Normal'); set(gca, 'XDir','reverse'); axis equal;
            end
            if opt == 2;
                tmp = image(Obj(ii).pos{1}, Obj(ii).pos{3}, Obj(ii).frame{opt});
                if opt3; hand{opt}(ii) = tmp; end
                set(tmp,'AlphaData', ~isnan(Obj(ii).frame{opt}(:,:,1))*Obj(ii).Trans);
                set(axx, 'YDir','Normal'); set(gca, 'XDir','Normal');axis equal;
            end
            if opt == 3;
                tmp = image(Obj(ii).pos{1}, Obj(ii).pos{2}, Obj(ii).frame{opt});
                if opt3; hand{opt}(ii) = tmp; end
                set(tmp,'AlphaData', ~isnan(Obj(ii).frame{opt}(:,:,1))*Obj(ii).Trans);
                set(axx, 'YDir','Normal'); set(gca, 'XDir','Normal');axis equal;
            end
        end
        set(axx,'XTick',[],'YTick',[]);
        axis equal;
        axis tight;
    end

    function out = axLim(dat,h)

        t1 = [1 1 1 1; size(dat,1) 1 1 1]*h.mat';
        out{1} = t1(1,1):(t1(2,1)-t1(1,1))/(size(dat,1)-1):t1(2,1);
        t1 = [1 1 1 1; 1 size(dat,2) 1 1]*h.mat';
        out{2} = t1(1,2):(t1(2,2)-t1(1,2))/(size(dat,2)-1):t1(2,2);
        t1 = [1 1 1 1; 1 1 size(dat,3) 1]*h.mat';
        out{3} = t1(1,3):(t1(2,3)-t1(1,3))/(size(dat,3)-1):t1(2,3);
        
        out{1} = out{1}(end:-1:1);
    end

    function [h1 h2] = crossHairs(axx,loc)
        axes(axx)
        ax = axis;
        h1 = plot([loc(1) loc(1)], ax(3:4),'b');
        h2 = plot(ax(1:2),[loc(2) loc(2)],'b'); axis equal; axis tight;
    end

    function setupFrames(n,opt)
        for ii = n;
            ss = Obj(ii).axLims;            

            if Obj(ii).point(1)~=Obj(ii).lastpoint(1) || opt==1
                tmp = Obj(ii).I(Obj(ii).point(1),:,:);
                tmp(Obj(ii).mask(Obj(ii).point(1),:,:)==0)=NaN;
                if numel(Obj(ii).Thresh)==2
                    ind = find(tmp<Obj(ii).Thresh(1) | tmp>Obj(ii).Thresh(2));
                else
                    ind = find(tmp<Obj(ii).Thresh(1) | tmp>Obj(ii).Thresh(4) | (tmp>Obj(ii).Thresh(2) & tmp<Obj(ii).Thresh(3)));
                end
                tmp(ind) = NaN;
                tmp(tmp==0)=NaN;
                
                tmp = flipdim(rot90(squeeze(tmp),1),1);
                [cols cm cc] = cmap(tmp, Obj(ii).clim, cmaps{Obj(ii).col});
                Obj(ii).frame{1} = reshape(cols,[size(tmp) 3]);
            end
            
            if Obj(ii).point(2)~=Obj(ii).lastpoint(2) || opt==1
                tmp = Obj(ii).I(:,Obj(ii).point(2),:);
                tmp(Obj(ii).mask(:,Obj(ii).point(2),:)==0)=NaN;
                if numel(Obj(ii).Thresh)==2
                    ind = find(tmp<Obj(ii).Thresh(1) | tmp>Obj(ii).Thresh(2));
                else
                    ind = find(tmp<Obj(ii).Thresh(1) | tmp>Obj(ii).Thresh(4) | (tmp>Obj(ii).Thresh(2) & tmp<Obj(ii).Thresh(3)));
                end
                tmp(ind) = NaN;
                tmp(tmp==0)=NaN;
                
                tmp = flipdim(flipdim(rot90(squeeze(tmp),1),1),2);
                [cols cm cc] = cmap(tmp, Obj(ii).clim, cmaps{Obj(ii).col});
                Obj(ii).frame{2} = reshape(cols,[size(tmp) 3]);
            end
            
            if Obj(ii).point(3)~=Obj(ii).lastpoint(3) || opt==1
                tmp = Obj(ii).I(:,:,Obj(ii).point(3));
                tmp(Obj(ii).mask(:,:,Obj(ii).point(3))==0)=NaN;
                if numel(Obj(ii).Thresh)==2
                    ind = find(tmp<Obj(ii).Thresh(1) | tmp>Obj(ii).Thresh(2));
                else
                    ind = find(tmp<Obj(ii).Thresh(1) | tmp>Obj(ii).Thresh(4) | (tmp>Obj(ii).Thresh(2) & tmp<Obj(ii).Thresh(3)));
                end
                tmp(ind) = NaN;
                tmp(tmp==0)=NaN;
                
                tmp = flipdim(flipdim(rot90(squeeze(tmp),1),1),2);
                [cols cm cc] = cmap(tmp, Obj(ii).clim, cmaps{Obj(ii).col});
                Obj(ii).frame{3} = reshape(cols,[size(tmp) 3]);
            end
        end
    end

    function toggle(varargin)
        state = get(varargin{1},'Checked');
        if strcmpi(state,'on');
            set(varargin{1},'Checked','off');
        end
        
        if strcmpi(state,'off');
            set(varargin{1},'Checked','on');
        end
    end

    function toggleCrossHairs(varargin)
        state = get(varargin{1},'Checked');
        if strcmpi(state,'on');
            set(ch(:),'Visible','off')
            set(varargin{1},'Checked','off');
        end
        
        if strcmpi(state,'off');
            set(ch(:),'Visible','on')
            set(varargin{1},'Checked','on');
        end
    end

    function newFig(varargin)
        tf = gcf;
        vn = get(con(21,1),'Value');
        if vn == 1
            return;
        end
        
        NewObj = FIVE({Obj(vn).FullPath},ch);
        
        tmp1 = get(pane,'Position');
        tmp2 = tmp1;
        tmp2(1) = tmp1(1)+tmp1(3);
        set(gcf,'Position',tmp2);

%         links{end+1}=linkprop([ch(1,1),NewObj(1).ch(1,1)],'XData');
%         links{end+1}=linkprop([ch(1,2),NewObj(1).ch(1,2)],'YData');
%         links{end+1}=linkprop([ch(2,1),NewObj(1).ch(2,1)],'XData');
%         links{end+1}=linkprop([ch(2,2),NewObj(1).ch(2,2)],'YData');
%         links{end+1}=linkprop([ch(3,1),NewObj(1).ch(3,1)],'XData');
%         links{end+1}=linkprop([ch(3,2),NewObj(1).ch(3,2)],'YData');
%         
%         set(tf,'UserData',links);
        
        set(menu(4),'Enable','on','Checked','on');
        figure(tf);
        removeVolume
    end

    function syncViews(varargin)
        state = get(varargin{1},'Checked');
        
        %%% Remove any existing links this figure has to other figures
        for ii = 1:length(links);
            set(links{ii},'Enabled','on');
            removeprop(links{ii},'XData');
            removeprop(links{ii},'YData');
        end
        links = [];
        
        
        %%% Establish links between this figure and all other FIVE
        %%% Figures
        ind = findobj(0,'Type','Figure');
        try
            figs = ind(searchCellStr('FIVE',get(ind,'Name')));
        catch            
            set(varargin{1},'Checked','off');
            set(gcf,'WindowButtonMotionFcn',@buttonMotion);
            return
        end
        figs = setdiff(figs,gcf);
        for ii = 1:length(figs);
            hhh = findobj(figs(ii),'Color','b','type','line');
            if isempty(hhh)
                hhh = findobj(figs(ii),'Color','y','type','line');
            end
            links{end+1}=linkprop([ch(1,1),hhh(6)],'XData');
            links{end+1}=linkprop([ch(1,2),hhh(5)],'YData');
            links{end+1}=linkprop([ch(2,1),hhh(4)],'XData');
            links{end+1}=linkprop([ch(2,2),hhh(3)],'YData');
            links{end+1}=linkprop([ch(3,1),hhh(2)],'XData');
            links{end+1}=linkprop([ch(3,2),hhh(1)],'YData');
           
        end
        set(gcf,'UserData',links);
        
        
        %%% Set the state of the links
        if strcmpi(state,'on');
            for ii = 1:length(links);
                set(links{ii},'Enabled','off');
            end
            set(varargin{1},'Checked','off');
            set(gcf,'WindowButtonMotionFcn',@buttonMotion);
        end
        
        
        if strcmpi(state,'off');
            for ii = 1:length(links);
                set(links{ii},'Enabled','on');
            end
            set(varargin{1},'Checked','on');
            set(gcf,'WindowButtonMotionFcn',[]);
        end
        
        
        %%% Propogate changes to other Figure links if they exists.
        for ii = 1:length(figs)
            lnk = get(figs(ii),'UserData');
            hh = findobj(figs(ii),'Label','Sync Views');
            if strcmpi(state,'on');
                for jj = 1:length(lnk);
                    set(lnk{jj},'Enabled','off');
                end
                set(hh,'Checked','off');
                set(figs(ii),'WindowButtonMotionFcn',@buttonMotion);
            end
            if strcmpi(state,'off');
                for jj = 1:length(lnk);
                    set(lnk{jj},'Enabled','on');
                end
                set(hh,'Checked','on');
                set(figs(ii),'WindowButtonMotionFcn',[]);
            end
        end
    end

    function changeLayer(varargin)
        vn = get(con(21,1),'Value');
        if vn == 1
            return
        end
        if varargin{1} == con(12,1)
            for ii = 1:length(hand)
                uistack(hand{ii}(vn),'top');
                uistack(hand{ii}(vn),'down');
                uistack(hand{ii}(vn),'down');
            end
        elseif varargin{1} == con(13,1)
            for ii = 1:length(hand)
                uistack(hand{ii}(vn),'bottom');
                uistack(hand{ii}(vn),'up');
            end
        end
    end

    function axialView(varargin)
        if numel(varargin{1})==1 || isempty(varargin)
            out = popup('Which Slices (MNI coordinates)?');
        else
            out = varargin{1};
        end
        
        if isempty(out)
            n = 25;
            coords = min(Obj(2).pos{3}):6:max(Obj(2).pos{3});
            c = 0;
            while numel(coords)>n
                c = c+1;
                if mod(c,2)==1;
                    coords = coords(2:end);
                else
                    coords = coords(1:end-1);
                end
            end
        else
            n = numel(out);
            coords = out;
        end       
        
        IN = [];
        
        for jj = 1:numel(Obj);
           IN.IM{jj} = Obj(jj).I.*double(Obj(jj).mask);
           IN.H{jj} = Obj(jj).h;
           IN.TH{jj} = Obj(jj).Thresh;
           IN.LIMS{jj} = Obj(jj).clim;
           IN.TRANS{jj} = Obj(jj).Trans;
           IN.CM{jj} = cmaps{Obj(jj).col};
        end       
        IN.Coords = coords;
        IN.opt = 1;
        SliceView(IN);
    end

    function coronalView(varargin)
        if numel(varargin{1})==1 || isempty(varargin)
            out = popup('Which Slices (MNI coordinates)?');
        else
            out = varargin{1};
        end
                
        if isempty(out)
            n = 30;
            coords = min(Obj(2).pos{2}):6:max(Obj(2).pos{2});
            c = 0;
            while numel(coords)>n
                c = c+1;
                if mod(c,2)==1;
                    coords = coords(2:end);
                else
                    coords = coords(1:end-1);
                end
            end
        else
            n = numel(out);
            coords = out;
        end       
        
        IN = [];
        for ii = 1:numel(Obj);
           IN.IM{ii} = Obj(ii).I.*double(Obj(ii).mask);
           IN.H{ii} = Obj(ii).h;
           IN.TH{ii} = Obj(ii).Thresh;
           IN.LIMS{ii} = Obj(ii).clim;
           IN.TRANS{ii} = Obj(ii).Trans;
           IN.CM{ii} = cmaps{Obj(ii).col};
        end       
        IN.Coords = coords;
        IN.opt = 3;
        SliceView(IN);
    end

    function sagittalView(varargin)
        if numel(varargin{1})==1 || isempty(varargin)
            out = popup('Which Slices (MNI coordinates)?');
        else
            out = varargin{1};
        end
        
        if isempty(out)
            n = 25;
            coords = min(Obj(2).pos{1}):6:max(Obj(2).pos{1});
            c = 0;
            while numel(coords)>n
                c = c+1;
                if mod(c,2)==1;
                    coords = coords(2:end);
                else
                    coords = coords(1:end-1);
                end
            end
        else
            n = numel(out);
            coords = out;
        end       
        
        IN = [];
        for ii = 1:numel(Obj);
           IN.IM{ii} = Obj(ii).I.*double(Obj(ii).mask);
           IN.H{ii} = Obj(ii).h;
           IN.TH{ii} = Obj(ii).Thresh;
           IN.LIMS{ii} = Obj(ii).clim;
           IN.TRANS{ii} = Obj(ii).Trans;
           IN.CM{ii} = cmaps{Obj(ii).col};
        end       
        IN.Coords = coords;
        IN.opt = 2;
        SliceView(IN);
    end

    function allSliceView(varargin)
        axialView;
        coronalView;
        sagittalView;
    end

    function surfView(varargin)
        if numel(Obj)==1
            return
        end
        
        for jj = 2:numel(Obj)
            vn = jj;
            obj = [];
            
            fs = get(findobj(paramenu3(12),'Checked','On'),'Label');
            switch fs
                case '642'
                    obj.fsaverage = 'fsaverage3';
                case '2562'
                    obj.fsaverage = 'fsaverage4';
                case '10242'
                    obj.fsaverage = 'fsaverage5';
                case '40962'
                    obj.fsaverage = 'fsaverage6';
                case '163842'
                    obj.fsaverage = 'fsaverage';
                otherwise
            end
            
            obj.surface = lower(get(findobj(paramenu3(10),'Checked','On'),'Label'));

            if strcmpi(get(paramenu3(19),'Checked'),'off')
                obj.shading = lower(get(findobj(paramenu3(11),'Checked','On'),'Label'));
                obj.shadingrange = [str2num(get(findobj(paramenu3(16),'Checked','On'),'Label')) str2num(get(findobj(paramenu3(17),'Checked','On'),'Label'))];
            else
                switch obj.surface
                    case 'white'
                        obj.shading = 'curv';
                        obj.shadingrange = [-2 2];
                    case 'pi'
                        obj.shading = 'mixed';
                        obj.shadingrange = [-2 2];
                    case 'inflated'
                        obj.shading = 'logcurv';
                        obj.shadingrange = [-.75 .75];
                    case 'pial'
                        obj.shading = 'curv';
                        obj.shadingrange = [-2 3];
                    otherwise
                        obj.shading = 'curv';
                        obj.shadingrange = [-2 3];
                end
            end
            
            obj.input.m = Obj(vn).I.*double(Obj(vn).mask);
            obj.input.he = Obj(vn).h;
            
            if jj ==2
                obj.figno = 0;
                obj.newfig = 1;
            else
                obj.figno = gcf;
                obj.newfig = 0;
            end
            obj.colorlims = Obj(vn).clim;
            cm = get(con(1),'String');
            cm = cm{Obj(vn).col+1};
            obj.colomap = cm;
            
            thresh = Obj(vn).Thresh;
            if numel(thresh) == 2
                [t,i] = min(abs(thresh));
                sgn = sign(thresh(i));
                t = t*sgn;
                obj.overlaythresh = t;
                if sgn==1
                    obj.direction = '+';
                    obj.reverse = 0;
                else
                    obj.direction = '-';
                    obj.reverse = 0;
                end
                
                if sum(sign(thresh))==0
                    obj.direction='+';
                end
            else
                obj.overlaythresh = Obj(vn).Thresh(2:3);
                obj.reverse = 0;
                obj.direction = '+';
            end
            
            obj.mappingfile = [];
            
            nn = get(findobj(paramenu3(15),'Checked','On'),'Label');
            if strcmpi(nn,'Yes')
                obj.nearestneighbor=1;
            else
                obj.nearestneighbor=0;
            end
            
            %load(which('SperStandard_4p5prox_FS6_SM.mat'),'header')
            %compH = Obj(vn).h.mat==header.mat;
            %if mean(compH(:))==1
            %   obj.mappingfile = which('SperStandard_4p5prox_FS6_SM.mat');
            %end
            
            
            option = get(findobj(paramenu3(14),'Checked','On'),'Label');
            switch option
                case 'All'
                    obj.Nsurfs = 4;
                case 'Both'
                    obj.Nsurfs = 2;
                case 'Left 1'
                    obj.Nsurfs = -1;
                case 'Left 2'
                    obj.Nsurfs = 1.9;
                case 'Right 1'
                    obj.Nsurfs = 1;
                case 'Right 2'
                    obj.Nsurfs = 2.1;
                otherwise
            end
            
            
            [h1 hh1] = surfPlot(obj);
            Obj(vn).SurfInfo.h1 = h1;
            Obj(vn).SurfInfo.hh1 = hh1;
            Obj(vn).SurfInfo.par = obj;
            
            set(hh1,'FaceAlpha',Obj(vn).Trans);
        end
    end

    function changeSign(varargin)
        vn = get(con(21,1),'Value');
        if vn == 1
            return
        end
        
        Obj(vn).I = Obj(vn).I*-1;
        
        t1 = num2str(sort(str2num(get(con(4),'String'))*-1));
        for ii = 1:10; t1 = regexprep(t1,'  ',' '); end
        set(con(4),'String',t1);
        
        t2 = num2str(sort(str2num(get(con(5),'String'))*-1));
        for ii = 1:10; t2 = regexprep(t2,'  ',' '); end
        set(con(5),'String',t2)        
        
        Obj(vn).Range(1:2) = Obj(vn).Range([2 1])*-1;
        
        UpdateThreshold;
        
        %setupFrames(vn,1);
        %updateGraphics({[1 2 3] vn},1);
        
    end

    function setupParamMenu
        %%% Not Specified
        %           out: output prefix, default is to define using imagefile
        %           SPM: 0 or 1, see above for details
        %          mask: optional to mask your data
        %           df1: numerator degrees of freedom for T/F-test (if 0<thresh<1)
        %           df2: denominator degrees of freedom for F-test (if 0<thresh<1)
        %        thresh: T/F statistic or p-value to threshold the data or 0
        
        paramenu1(1) = uimenu(pane,'Label','Parameters');
        %paramenu2(1) = uimenu(paramenu1(1),'Label','Cluster Params');
        
        paramenu3(1) = uimenu(paramenu1(1),'Label','Sign');
        paramenu4(1) =     uimenu(paramenu3(1),'Label','Pos','Checked','on','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(1),'Label','Neg','CallBack',@groupCheck);
        
        paramenu3(2) = uimenu(paramenu1(1),'Label','Sphere Radius');
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','1mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','2mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','3mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','4mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','5mm','Checked','On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','6mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','7mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','8mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','9mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','10mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','11mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','12mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','13mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','14mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(2),'Label','15mm','CallBack',@groupCheck);
        
        paramenu3(3) = uimenu(paramenu1(1), 'Label','Peak Number Limit');
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','100','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','200','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','500','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','750','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','1000','Checked','on','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','1500','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','2000','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(3),'Label','3000','CallBack',@groupCheck);
        
        paramenu3(4) = uimenu(paramenu1(1),'Label','Peak Separation');
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','2mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','4mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','6mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','8mm','Checked', 'on','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','10mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','12mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','14mm','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(4),'Label','16mm','CallBack',@groupCheck);
        
        paramenu3(5) = uimenu(paramenu1(1),'Label','Neighbor Def.');
        paramenu4(end+1) = uimenu(paramenu3(5),'Label','6','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(5),'Label','18','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(5),'Label','26','Checked','on','CallBack',@groupCheck);
        
        
        paramenu3(7) = uimenu(paramenu1(1),'Label','Labeling');
        paramenu4(end+1) = uimenu(paramenu3(7),'Label','Use Nearest Label', 'Checked', 'On','CallBack',@groupCheck);
        
%         paramenu3(8) = uimenu(paramenu1(1),'Label','Label Map');
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','aal_MNI_V4', 'Checked', 'On','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','Nitschke_Lab','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','JHU_tracts','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','JHU_whitematter','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','Talairach','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','Thalamus','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','MNI','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','HarvardOxford_cortex','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','Cerebellum-flirt','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','Cerebellum-fnirt','CallBack',@changeLabelMap);
%         paramenu4(end+1) = uimenu(paramenu3(8),'Label','Juelich','CallBack',@changeLabelMap);
        
        paramenu3(8) = uimenu(paramenu1(1),'Label','Label Map');
        paramenu4(end+1) = uimenu(paramenu3(8),'Label','aal_MNI_V4', 'Checked', 'On','CallBack',@changeLabelMap);
        try
            [labellist]=getLabelMap();
            for ll=2:numel(labellist)
                paramenu4(end+1) = uimenu(paramenu3(8),'Label',labellist{ll},'CallBack',@changeLabelMap);
            end
        catch
            paramenu4(end+1) = uimenu(paramenu3(8),'Label','aal_MNI_V4', 'Checked', 'On','CallBack',@changeLabelMap);
            warning('Download Peak Nii to enable the use of other atlases');
        end
        paramenu3(18) = uimenu(paramenu1(1),'Label','Corrected Alpha');
        paramenu4(end+1) = uimenu(paramenu3(18),'Label','0.100', 'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(18),'Label','0.050', 'Checked', 'On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(18),'Label','0.025', 'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(18),'Label','0.010', 'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(18),'Label','0.001', 'Checked', 'Off','CallBack',@groupCheck);
        
        
        paramenu1(2) = uimenu(pane,'Label','Region Name');
        
        %menu(?) = uimenu(menu(1),'Label','Increase FontSize','CallBack',@increaseFont);
        %menu(?) = uimenu(menu(1),'Label','Decrease FontSize','CallBack',@decreaseFont);
        
        
        paramenu3(9) = uimenu(paramenu1(1),'Label','Surface Options');
        paramenu3(19) = uimenu(paramenu3(9),'Label','Use Surface Shading Defaults','Checked', 'On','CallBack',@groupCheck);

        paramenu3(10) = uimenu(paramenu3(9),'Label','Surface');
        paramenu4(end+1) = uimenu(paramenu3(10),'Label','Inflated', 'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(10),'Label','Pial',     'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(10),'Label','White',    'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(10),'Label','PI',       'Checked', 'On','CallBack',@groupCheck);
        
        paramenu3(11) = uimenu(paramenu3(9),'Label','Shading');
        paramenu4(end+1) = uimenu(paramenu3(11),'Label','Curv',     'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(11),'Label','Sulc',     'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(11),'Label','Thk',      'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(11),'Label','LogCurv',  'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(11),'Label','Mixed',    'Checked', 'On','CallBack',@groupCheck);
        
        paramenu3(12) = uimenu(paramenu3(9),'Label','N-verts');
        paramenu4(end+1) = uimenu(paramenu3(12),'Label','642',    'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(12),'Label','2562',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(12),'Label','10242',  'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(12),'Label','40962',  'Checked', 'On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(12),'Label','163842', 'Checked', 'Off','CallBack',@groupCheck);
        
        % paramenu3(13) = uimenu(paramenu3(9),'Label','Reverse Map');
        % paramenu4(end+1) = uimenu(paramenu3(13),'Label','No',    'Checked', 'On','CallBack',@groupCheck);
        % paramenu4(end+1) = uimenu(paramenu3(13),'Label','Yes',   'Checked', 'Off','CallBack',@groupCheck);

        paramenu3(14) = uimenu(paramenu3(9),'Label','N-Surfs');
        paramenu4(end+1) = uimenu(paramenu3(14),'Label','All',      'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(14),'Label','Both',     'Checked', 'On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(14),'Label','Left 1',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(14),'Label','Left 2',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(14),'Label','Right 1',  'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(14),'Label','Right 2',  'Checked', 'Off','CallBack',@groupCheck);
        
        paramenu3(15) = uimenu(paramenu3(9),'Label','Nearest Neighbor Only');
        paramenu4(end+1) = uimenu(paramenu3(15),'Label','No',    'Checked', 'On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(15),'Label','Yes',   'Checked', 'Off','CallBack',@groupCheck);
        
        
        paramenu3(16) = uimenu(paramenu3(9),'Label','Underlay Min');
        paramenu4(end+1) = uimenu(paramenu3(16),'Label',' 0.0',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(16),'Label','-0.5',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(16),'Label','-1.0',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(16),'Label','-1.5',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(16),'Label','-2.0',   'Checked', 'On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(16),'Label','-2.5',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(16),'Label','-3.0',   'Checked', 'Off','CallBack',@groupCheck);
        
        paramenu3(17) = uimenu(paramenu3(9),'Label','Underlay Max');
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+0.0',    'Checked','Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+0.5',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+1.0',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+1.5',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+2.0',   'Checked', 'On','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+2.5',   'Checked', 'Off','CallBack',@groupCheck);
        paramenu4(end+1) = uimenu(paramenu3(17),'Label','+3.0',   'Checked', 'Off','CallBack',@groupCheck);
        
    end

    function setupFigure
        
        tss = Obj(1).axLims;
        
        im = spm_imatrix(Obj(1).h.mat);
        %tmp = [1 1 1];
        tmp = abs(im(7:9));
        tss = tss.*tmp;
        
        height = tss(2)+tss(3);
        width =  tss(1)+tss(2);
        rat = height/width;
        
        ff = get(0,'ScreenSize');
        
        if ff(4)>800
            pro = 3.25;
        else
            pro = 2.5;
        end
        %[junk user] = UserTime;  if strcmp(user,'aschultz'); keyboard; end
        
        ss = get(0,'ScreenSize');
        if (ss(3)/ss(4))>2
            ss(3)=ss(3)/2;
        end
        op = floor([50 ss(4)-75-((ss(3)/pro)*rat)  ss(3)/pro (ss(3)/pro)*rat]);
        
        pane = figure('Visible','off');
        
        set(gcf, 'Position', op,'toolbar','none', 'Name', 'FIVE','Visible','off');
        set(gcf, 'WindowButtonUpFcn', @buttonUp);
        set(gcf, 'WindowButtonDownFcn', @buttonDown);
        set(gcf, 'WindowButtonMotionFcn', @buttonMotion);
        set(gcf, 'ResizeFcn', @resizeFig);
        
        hcmenu = uicontextmenu;
        set(hcmenu,'CallBack','movego = 0;');
        item = [];
        
        wid(1) = tss(2)/width;
        hei(1) = tss(3)/height;
        wid(2) = tss(1)/width;
        hei(2) = tss(3)/height;
        wid(3) = tss(1)/width;
        hei(3) = tss(2)/height;
        
        ax1 = axes; set(ax1,'Color','k','Position',[wid(2) hei(3) wid(1) hei(1)],'XTick',[],'YTick',[],'YColor','k','XColor','k'); hold on; colormap(gray(256));
        ax2 = axes; set(ax2,'Color','k','Position',[0      hei(3) wid(2) hei(2)],'XTick',[],'YTick',[],'YColor','k','XColor','k'); hold on; colormap(gray(256));
        ax3 = axes; set(ax3,'Color','k','Position',[0      0      wid(3) hei(3)],'XTick',[],'YTick',[],'YColor','k','XColor','k'); hold on; colormap(gray(256));
        ax4 = axes; set(ax4,'Color','w','Position',[wid(2) 0      .04     hei(3)*.99],'XTick',[],'YAxisLocation','right','YTick',[]); %,'YColor','k','XColor','k'
        
        
        set(gcf,'WindowKeyPressFcn', @keyMove);
        set(gcf,'WindowScrollWheelFcn',@scrollMove);
        set(gcf,'WindowKeyReleaseFcn',@keyHandler);
        
        %st = wid(2)+.01+.04+.03;
        st = wid(2)+.125;
        len1 = (1-st)-.01;
        len2 = len1/2;
        len3 = len1/3;
        len4 = len1/4;
        
        inc = (hei(3))/11;
        inc2 = .75*inc;
        
        % [wid(2)+.045      0  .025   hei(3)*.9]
        % [wid(2)+.070      0  .025   hei(3)*.9]
        % [wid(2)+.095      0  .025   hei(3)*.9]

        val = (Obj(1).Range(2)-Obj(1).Range(1))/Obj(1).Range(3);
        con(27,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[wid(2)+.070      0  .025   hei(3)*.9],'Value',val,'CallBack', @adjustUnderlay);
        val = (Obj(1).Range(1)-Obj(1).Range(1))/Obj(1).Range(3);
        con(28,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[wid(2)+.095      0  .025   hei(3)*.9],'Value',val,'CallBack', @adjustUnderlay);
        %con(29,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[wid(2)+.095      0  .025   hei(3)*.9],'Value',0,'CallBack', @adjustUnderlay);
        
        con(30,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[wid(2)+.070    hei(3)*.9   .05 .05]); % wid(2)+.0450    hei(3)*.9   .075 .05]
        
        % con(28,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[st-.1          (.01+(1*inc)+(0*inc2))  len1*.75   inc],'Value',1,'CallBack', @adjustTrans);
        
        con(1,1)   = uicontrol(pane,'style','popupmenu', 'Units','Normalized','Position',[st          .01+(0*inc)+(0*inc2)    len1   inc],'String',[{'Colormap'} cmaps(:)'],'CallBack', @changeColorMap); shg
        %con(1,1)   = uicontrol(pane,'style','popupmenu', 'Units','Normalized','Position',[st          .01+(0*inc)+(0*inc2)    len1   inc],'String',[{'Colormap' 'A' 'B' 'C' 'D'}],'CallBack', @changeColorMap); shg
        con(2,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[st          (.01+(1*inc)+(0*inc2))  len1*.75   inc],'Value',1,'CallBack', @adjustTrans);
        con(23,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+len1*.75 (.01+(1*inc)+(0*inc2))  len1*.25   inc],'String','All','CallBack', @applyToAll);
        con(3,1)   = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(2*inc)+(0*inc2)  len1   inc2],'String', 'Transparency','fontsize',12);
        
        con(4,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(2*inc)+(1*inc2)  len2   inc],'CallBack',@UpdateThreshold);
        con(5,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st+len2     .01+(2*inc)+(1*inc2)  len2   inc],'CallBack',@UpdateCLims);
        con(6,1)   = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(3*inc)+(1*inc2)  len2   inc2],'String', 'Thresh','fontsize',12);
        con(7,1)   = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len2     .01+(3*inc)+(1*inc2)  len2   inc2],'String', 'Color Limits','fontsize',12);
        con(8,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(3*inc)+(2*inc2)  len3   inc],'CallBack',@UpdatePVal);
        con(9,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st+len3     .01+(3*inc)+(2*inc2)  len3   inc],'CallBack',@UpdatePVal);
        con(23,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st+len3+len3     .01+(3*inc)+(2*inc2)  len3   inc],'CallBack',@ExtentThresh);
        con(10,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(4*inc)+(2*inc2)  len3   inc2],'String', 'DF','fontsize',12);
        con(11,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len3     .01+(4*inc)+(2*inc2)  len3   inc2],'String', 'P-Value','fontsize',12);
        con(24,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len3+len3     .01+(4*inc)+(2*inc2)  len3   inc2],'String', 'Extent','fontsize',12);
        
        con(15,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(4*inc)+(3*inc2)  len2   inc],'String',num2str(loc),'CallBack',@goTo);
        %con(15,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(4*inc)+(3*inc2)  len2   inc],'String',['0 0 0'],'CallBack',@goTo);
        con(25,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+len2     .01+(4*inc)+(3*inc2)  len2/2   inc],'String','FDR','CallBack', @correctThresh);
        con(26,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+(len2*1.5)  .01+(4*inc)+(3*inc2)  len2/2   inc],'String','FWE','CallBack', @correctThresh);
        con(17,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(5*inc)+(3*inc2)  len2   inc2],'String', 'MNI Coord',  'fontsize',12);
        con(18,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len2     .01+(5*inc)+(3*inc2)  len2   inc2],'String', 'MC Correct','fontsize',12);
        con(12,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st          .01+(5*inc)+(4*inc2) len2   inc],'String',{'Move Up'},'CallBack',@changeLayer);
        con(13,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+len2     .01+(5*inc)+(4*inc2) len2   inc],'String',{'Move Down'},'CallBack',@changeLayer);
        con(21,1)  = uicontrol(pane,'style','popupmenu', 'Units','Normalized','Position',[st          .01+(5*inc)+(5*inc2) len1   inc],'String',{'Overlays'},'CallBack',@switchObj);
        con(22,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[wid(2)-(len2/2) hei(3)-(inc*1.05) (len2/2) inc]);
        
        
        con(19,1)  = uicontrol(pane,'style','PushButton','Units','Normalized','Position',[st          .01+(6*inc)+(5*inc2) len2   inc],'String','Open Overlay','FontWeight','Bold','CallBack',@openOverlay);
        con(20,1)  = uicontrol(pane,'style','PushButton','Units','Normalized','Position',[st+len2     .01+(6*inc)+(5*inc2) len2   inc],'String','Remove Volume', 'FontWeight','Bold','CallBack',@removeVolume);
        
        
        
        menu(1) = uimenu(pane,'Label','Options');
        
        
        menu(2) = uimenu(menu(1),'Label','CrossHair Toggle','Checked','on','CallBack', @toggleCrossHairs);
        
        menu(43) = uimenu(menu(1),'Label','Reverse Image','CallBack', @changeSign);
        
        
        menu(15) = uimenu(menu(1),'Label','Change Underlay','CallBack', @changeUnderlay);
        menu(3) = uimenu(menu(1),'Label','Send Overlay To New Fig','CallBack',@newFig);
        if nargin<2
            menu(4) = uimenu(menu(1),'Label','Sync Views','Enable','on','CallBack',@syncViews);
        else
            menu(4) = uimenu(menu(1),'Label','Sync Views','Enable','on','Checked','on','CallBack',@syncViews);
        end
        menu(16) = uimenu(menu(1),'Label','SliceViews');
        menu(6) = uimenu(menu(16),'Label','Axial Slice View','CallBack',@axialView);
        menu(7) = uimenu(menu(16),'Label','Coronal Slice View','CallBack',@coronalView);
        menu(8) = uimenu(menu(16),'Label','Sagittal Slice View','CallBack',@sagittalView);
        menu(9) = uimenu(menu(16),'Label','All Slice View','CallBack',@allSliceView);
        
        
        
        menu(17) = uimenu(menu(1),'Label','Ploting');
        menu(10) = uimenu(menu(17),'Label','Load Plot Data','CallBack',@loadData);
        menu(13) = uimenu(menu(17),'Label','JointPlot','CallBack',@JointPlot);
        menu(12) = uimenu(menu(17),'Label','PaperPlot-Vert','CallBack',@PaperFigure_Vert);
        menu(33) = uimenu(menu(17),'Label','PaperPlot-Horz','CallBack',@PaperFigure_Horz);
        menu(23) = uimenu(menu(17),'Label','PaperPlot Overlay','Checked','off','CallBack', @toggle);
        menu(41) = uimenu(menu(17),'Label','Surface Render','Checked','off','CallBack', @surfView);
        
        menu(14) = uimenu(menu(1),'Label','Transparent Overlay','CallBack',@TestFunc);
        menu(11) = uimenu(menu(1),'Label','Get Peak Info','CallBack',@getPeakInfo);
        
        menu(25) = uimenu(menu(1),'Label','Resample');
        menu(26) = uimenu(menu(25),'Label','.5x.5x.5', 'CallBack', @resampleIm);
        menu(27) = uimenu(menu(25),'Label','1x1x1'   , 'CallBack', @resampleIm);
        menu(28) = uimenu(menu(25),'Label','2x2x2'   , 'CallBack', @resampleIm);
        menu(29) = uimenu(menu(25),'Label','3x3x3'   , 'CallBack', @resampleIm);
        menu(30) = uimenu(menu(25),'Label','4x4x4'   , 'CallBack', @resampleIm);
        menu(31) = uimenu(menu(25),'Label','5x5x5'   , 'CallBack', @resampleIm);
        menu(32) = uimenu(menu(25),'Label','6x6x6'   , 'CallBack', @resampleIm);
        
        menu(18) = uimenu(menu(1),'Label','Save Options');
        menu(19) = uimenu(menu(18),'Label','Save Thresholded Image','Enable','on','CallBack',@saveImg);
        menu(20) = uimenu(menu(18),'Label','Save Masked Image','Enable','on','CallBack',@saveImg);
        menu(21) = uimenu(menu(18),'Label','Save Cluster Image','Enable','on','CallBack',@saveImg);
        menu(22) = uimenu(menu(18),'Label','Save Cluster Mask','Enable','on','CallBack',@saveImg);
        
        menu(34) = uimenu(menu(1),'Label','Mask');
        menu(35) = uimenu(menu(34),'Label','Mask In','Enable','on','CallBack',@maskImage);
        menu(36) = uimenu(menu(34),'Label','Mask Out','Enable','on','CallBack',@maskImage);
        menu(37) = uimenu(menu(34),'Label','Un-Mask','Enable','on','CallBack',@maskImage);
        
        menu(38) = uimenu(menu(1),'Label','Movie Mode','Enable','on','CallBack',@movieMode);
        menu(39) = uimenu(menu(1),'Label','Conn Explore','Enable','on','CallBack',@initializeConnExplore);
        menu(40) = uimenu(menu(1),'Label','SS Connectivity','Enable','on','CallBack',@ssConn);
        
        item(1) = uimenu(hcmenu, 'Label', 'Go to local max',  'Callback', @gotoMinMax);
        item(2) = uimenu(hcmenu, 'Label', 'Go to local min',  'Callback', @gotoMinMax);
        item(3) = uimenu(hcmenu, 'Label', 'Go to global max',  'Callback', @gotoMinMax);
        item(4) = uimenu(hcmenu, 'Label', 'Go to global min',  'Callback', @gotoMinMax);
        item(5) = uimenu(hcmenu, 'Label', 'Plot Cluster', 'Callback', @plotVOI);
        item(6) = uimenu(hcmenu, 'Label', 'Plot Sphere', 'Callback', @plotVOI);
        item(7) = uimenu(hcmenu, 'Label', 'Plot Voxel',  'Callback', @plotVOI);
        item(8) = uimenu(hcmenu, 'Label', 'Plot Cached Cluster',  'Callback', @plotVOI);
        item(9) = uimenu(hcmenu, 'Label', 'Cache Cluster Index',  'Callback', @CachedPlot);
        %item(8) = uimenu(hcmenu, 'Label', 'RegionName',  'Callback', @regionName);
        
        menu(42) = uimenu(menu(1),'Label','Return Obj','Checked','off','CallBack', @returnInfo);
        %menu(43) = uimenu(menu(1),'Label','Return Obj as Global','Checked','off','CallBack', @returnInfoGlob);
        
        Obj(1).ax1 = ax1;
        Obj(1).ax2 = ax2;
        Obj(1).ax3 = ax3;
        Obj(1).ax4 = ax4;
        Obj(1).con = con;
        Obj(1).menu = menu;
    end

    function returnInfo(varargin)
        if isstruct(varargin{1})
            Obj = varargin{1};
        else
            assignin('base','Obj',Obj);
        end
    end
        
    function keyHandler(varargin)
        d = varargin{2};
        %disp([char(d.Modifier) ' - ' d.Key]);
        if strcmpi('shift',d.Modifier)
            switch d.Key
                
                case 'o'
                    maskImage(menu(36));
                case 'i'
                    maskImage(menu(35));
                case 'u'
                    maskImage(menu(37));
            end
        else
            switch d.Key
                
                case 'u'
                    updateConnMap;
                case 's'
                    ssUpdateConnMap;    
                otherwise
            end
        end
        
        if searchCellStr('alt',  {varargin{2}.Key})==1
            clickFlag = 0;
        end
    end

    function maskImage(varargin)
        vn = get(con(21,1),'Value');
        
        if varargin{1}==menu(37);
            Obj(vn).MaskInd = [];
            Obj(vn).mask(:) = 1;
            Obj(vn).mask(Obj(vn).Exclude)=0;
            
            setupFrames(vn,1);
            updateGraphics({1:3 vn},1);
            return
        end
        
        nnn = spm_select(inf,'image','Select a Mask Image:');
        th = spm_vol(nnn);
        for ii = 1:numel(th)
            m = resizeVol(th(ii),Obj(vn).h);

            if varargin{1}==menu(35)
                Obj(vn).MaskInd = [Obj(vn).MaskInd; find(isnan(m))];
            end
            
            if varargin{1}==menu(36)
                Obj(vn).MaskInd = [Obj(vn).MaskInd; find(~isnan(m))];
            end
        end
        Obj(vn).MaskInd = unique(Obj(vn).MaskInd);
        
        
        Obj(vn).mask(:) = 1;
        Obj(vn).mask(Obj(vn).MaskInd)=0;
        Obj(vn).mask(Obj(vn).Exclude)=0;
        
        setupFrames(vn,1);
        updateGraphics({1:3 vn},1); 
    end

    function CachedPlot(varargin)        
        if varargin{1} == item(9)
            CachedClusterLoc = [];            
            vn = get(con(21,1),'Value');
            
            CM = str2num(get(findobj(paramenu3(5),'Checked','On'),'Label'));
            lastItem = item(5);
            adir = [fileparts(Obj(vn).FullPath) filesep];
            
            %%% Using the diplayed image get the matrix vector indices
            %%% for the selected cluster
            mni = str2num(get(con(15),'string'));
            
            if ~isempty(Obj(vn).mask);
                tmpI = Obj(vn).I.*double(Obj(vn).mask);
            else
                tmpI = Obj(vn).I;
            end
            
            thresh = Obj(vn).Thresh;
            if numel(thresh)==2
                ind = find(tmpI<=thresh(2) & tmpI>=thresh(1));
                L = [];
                [L(:,1) L(:,2) L(:,3)] = ind2sub(Obj(vn).h.dim,ind);
                A = spm_clusters2(L(:,1:3)',CM);
            elseif numel(thresh)==4
                ind1 = find(tmpI<=thresh(2) & tmpI>=thresh(1));
                L1 = [];
                [L1(:,1) L1(:,2) L1(:,3)] = ind2sub(Obj(vn).h.dim,ind1);
                A1 = spm_clusters2(L1(:,1:3)',CM);
                
                ind2 = find(tmpI<=thresh(4) & tmpI>=thresh(3));
                L2 = [];
                [L2(:,1) L2(:,2) L2(:,3)] = ind2sub(Obj(vn).h.dim,ind2);
                A2 = spm_clusters2(L2(:,1:3)',CM);
                
                L = [L1; L2];
                A = [A1 A2+max(A1)];
                ind = [ind1; ind2];
            end
            
            L(:,4) = 1;
            L2 = L*Obj(vn).h.mat';
            dist = sqrt(sum((L2(:,1:3)-repmat(mni, size(L,1),1)).^2,2));
            if min(dist)>5
                disp('No Cluster is selected');
            end
            ind1 = find(dist==min(dist));
            ind2 = find(A==A(ind1(1)));
            ind3 = ind(ind2);
            
            CachedClusterLoc.index = ind3;
            CachedClusterLoc.mat = Obj(vn).h.mat;
        end
    end

    function PaperFigure_Vert(varargin)
        %%% Figure width is off
        ss2 = get(gcf,'position');
        tss = size(Obj(1).I);
        tmp = abs(matInfo(Obj(1).h.mat));
        tss = tss.*tmp;
        
        tot = tss(3)+tss(3)+tss(2);
        
        wid(1) = tss(2)/width;
        hei(1) = tss(3)/tot;
        wid(2) = tss(1)/width;
        hei(2) = tss(3)/tot;
        wid(3) = tss(1)/width;
        hei(3) = tss(2)/tot;
        
        y = get(ch(1,1),'XData');
        x = get(ch(2,1),'XData');
        z = get(ch(1,2),'YData');
        xyz = [x(1) y(1) z(1)];


        fg = figure; clf; colormap(gray);
        set(fg,'color','k')
        n = numel(Obj)-1;
        if n==0;
            return;
        end        
        if strcmpi(get(menu(23),'Checked'),'off')
            ww = floor((max(wid)*ss2(3)) * n);
        else
            ww = floor((max(wid)*ss2(3)) * 1);
            n = 1;
        end
        hh = floor((sum(tss([3 3 2])./height)*ss2(4))*1);
        dims = [10 10 ww+((1/8)*ww) hh];
        dims(3) = floor(dims(3));
        set(fg,'Position',dims); 
        
        dd = [.1/n (((dims(3)-dims(1))/n)*.1)/((dims(4)-dims(2)))];
        for zz = 2:length(Obj)
            if (zz==2) || strcmpi(get(menu(23),'Checked'),'off')
                aax1 = axes; set(aax1,'Color','k','Position',[(zz-2)/n sum(hei(1:2)) (1/n)-((1/8)/n) hei(3)],'Xcolor','k','Ycolor','k'); hold on;
                aax2 = axes; set(aax2,'Color','k','Position',[(zz-2)/n sum(hei(2))   (1/n)-((1/8)/n) hei(1)],'Xcolor','k','Ycolor','k'); hold on;
                aax3 = axes; set(aax3,'Color','k','Position',[(zz-2)/n 0             (1/n)-((1/8)/n) hei(2)],'Xcolor','k','Ycolor','k'); hold on;
                aax4 = axes; set(aax4,'Color','k','Position',[((zz-2)/n)+((7/8)/n) 0 (1/8)/n 1],'Xcolor','k','Ycolor','k'); hold on;
                tr(zz-1) = aax4;
                
                drawFresh(aax1,3,1,0);
                drawFresh(aax2,1,1,0);
                drawFresh(aax3,2,1,0);
            end
            h_in = (zz-2)/(numel(Obj)-1);
            
            drawFresh(aax1,3,zz,0);
            drawFresh(aax2,1,zz,0);
            drawFresh(aax3,2,zz,0);
            if strcmpi(get(menu(23),'Checked'),'off')
                an(zz-1,1) = annotation(gcf,'textbox',[h_in+.005 sum(hei(1:2))+(.88*hei(3)) dd],...
                    'Color','w', 'String',num2str(xyz(3)),'fontsize',18,'fontweight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
            
                an(zz-1,2) = annotation(gcf,'textbox',[h_in+.005 sum(hei(2))+(.88*hei(2)) dd],...
                    'Color','w', 'String',num2str(xyz(1)),'fontsize',18,'fontweight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
           
                an(zz-1,3) = annotation(gcf,'textbox',[h_in+.005 0+(.88*hei(2)) dd],...
                    'Color','w', 'String',num2str(xyz(2)),'fontsize',18,'fontweight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
                
                vn = zz;
                lim = [min(Obj(vn).clim) max(Obj(vn).clim)];
                yy = lim(1):(lim(2)-lim(1))/255:lim(2);
                axes(aax4); cla
                
                cm = colmap(cmaps{Obj(vn).col},256);
                
                imagesc(yy',1,reshape(cm,numel(yy),1,3)); axis tight
                set(aax4,'YDir','Normal','YAxisLocation','left','YTick', unique([1 get(aax4,'YTick')]));
                
                set(aax4,'YTickLabel',(round(((min(yy):spm_range(yy)/(numel(get(aax4,'YTick'))-1):max(yy)))*100)/100));
                
                
                th = moveYaxLabs(aax4,'y');
                set(th,'fontsize',get(con(2),'fontsize'),'fontweight','bold');
                set(th,'BackgroundColor','k'); shg
            else
                if zz == 2
                    delete(aax4);
                    tr = [];
                end
            end
        end
        for zz = 1:length(tr)
            uistack(tr(zz),'top');
        end
        
        if strcmpi(get(menu(23),'Checked'),'on')
            h_in = 0;
            an(1,1) = annotation(gcf,'textbox',[h_in+.005 sum(hei(1:2))+(.88*hei(3)) dd],...
                'Color','w', 'String',num2str(xyz(3)),'fontsize',18,'fontweight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
            
            an(1,2) = annotation(gcf,'textbox',[h_in+.005 sum(hei(2))+(.88*hei(2)) dd],...
                'Color','w', 'String',num2str(xyz(1)),'fontsize',18,'fontweight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
            
            an(1,3) = annotation(gcf,'textbox',[h_in+.005 0+(.88*hei(2)) dd],...
                'Color','w', 'String',num2str(xyz(2)),'fontsize',18,'fontweight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
        end
    end

    function PaperFigure_Horz(varargin)
        scrn = get(0, 'ScreenSize');
        tmp = get(Obj(1).ax1,'position');
        wid(1) = tmp(3); hei(1) = tmp(4);
        tmp = get(Obj(1).ax2,'position');
        wid(2) = tmp(3); hei(2) = tmp(4);
        tmp = get(Obj(1).ax3,'position');
        wid(3) = tmp(3); hei(3) = tmp(4);
        
        ss2 = get(gcf,'position'); 
        wid = wid*ss2(3);
        hei = hei*ss2(4);
        
        
        n = numel(Obj)-1;
        if n==0;
            return;
        end
        
        %fg = figure(400+gcf); clf; colormap(gray)
        fg = figure; clf; colormap(gray)
        set(fg,'color','k')
        
        tall = max(hei)./.9;
        
        if strcmpi(get(menu(23),'Checked'),'off')
            hh = tall*n;
        else
            hh = tall;
            n = 1;
        end

        dims = [10 scrn(3)-75 sum(wid) hh]; 
        dims = round(dims);
        set(fg,'Position',dims); 
        
        wid = wid./sum(wid);
                
        y = get(ch(1,1),'XData');
        x = get(ch(2,1),'XData');
        z = get(ch(1,2),'YData');
        xyz = [x(1) y(1) z(1)];

        dd = [.1/n .033];
        
        for zz = 2:length(Obj)    
            if (zz==2) || strcmpi(get(menu(23),'Checked'),'off')
                aax1 = axes; set(aax1,'Color','k','Position',[0                (( n-zz+1)/n)+(.1/n)   wid(3)  .9/n],'Xcolor','k','Ycolor','k'); hold on;
                aax2 = axes; set(aax2,'Color','k','Position',[wid(3)           (( n-zz+1)/n)+(.1/n)   wid(2)  .9/n],'Xcolor','k','Ycolor','k'); hold on;
                aax3 = axes; set(aax3,'Color','k','Position',[sum(wid([2 3]))  (( n-zz+1)/n)+(.1/n)   wid(1)  .9/n],'Xcolor','k','Ycolor','k'); hold on;
                aax4 = axes; set(aax4,'Color','k','Position',[0                (( n-zz+1)/n)          1       .1/n],'Xcolor','k','Ycolor','k'); hold on;
                tr(zz-1) = aax4;
                
                drawFresh(aax1,3,1,0);
                drawFresh(aax2,2,1,0);
                drawFresh(aax3,1,1,0);
            end
            
            h_in = (zz-2)/(numel(Obj)-1);
            
            drawFresh(aax1,3,zz,0);
            drawFresh(aax2,2,zz,0);
            drawFresh(aax3,1,zz,0);
            
            if strcmpi(get(menu(23),'Checked'),'off')
                an(zz-1,1) = annotation(gcf,'textbox',[0 ((n-zz+1)/n)+(.1/n)+(.8/n) .02 .033],...
                    'Color','w', 'String',num2str(xyz(3)),'fontsize',18,'fontweight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
                an(zz-1,2) = annotation(gcf,'textbox',[wid(2) ((n-zz+1)/n)+(.1/n)+(.8/n) .02 .033],...
                    'Color','w', 'String',num2str(xyz(2)),'fontsize',18,'fontweight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
                an(zz-1,3) = annotation(gcf,'textbox',[sum(wid([2 3]))+(.04) ((n-zz+1)/n)+(.1/n)+(.8/n) .02 .033],...
                    'Color','w', 'String',num2str(xyz(1)),'fontsize',18,'fontweight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
                
                
                vn = zz;
                lim = [min(Obj(vn).clim) max(Obj(vn).clim)];
                yy = lim(1):(lim(2)-lim(1))/255:lim(2);
                axes(aax4); cla
                
                cm = colmap(cmaps{Obj(vn).col},256);
                imagesc(1,yy,reshape(cm,1,numel(yy),3)); axis tight
                %imagesc(1,yy',reshape(cmap{Obj(vn).col,1},1,size(cmap{Obj(vn).col,1},1),3)); axis tight
                
                set(aax4,'XDir','Normal','XAxisLocation','top','XTick', unique([1 get(aax4,'XTick')]));
                
                set(aax4,'XTickLabel',(round(((min(yy):spm_range(yy)/(numel(get(aax4,'XTick'))-1):max(yy)))*100)/100));
                
                set(aax4,'xcolor','w','fontsize',12,'fontweight','bold')
                
                th = moveYaxLabs(aax4,'x');
                %set(th,'fontsize',get(con(2),'fontsize'),'fontweight','bold');
                set(th,'BackgroundColor','k'); shg
            else
                if zz == 2
                    delete(aax4);
                    tr = [];
                end
            end
        end
        
        
        for zz = 1:length(tr)
            uistack(tr(zz),'top');
        end
        
        if strcmpi(get(menu(23),'Checked'),'on')
            h_in = 0;
            an(1,1) = annotation(gcf,'textbox',[0                     (.1/n)+(.8/n) .02 .033],...
                'Color','w', 'String',num2str(xyz(3)),'fontsize',18,'fontweight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
            an(1,2) = annotation(gcf,'textbox',[wid(2)                (.1/n)+(.8/n) .02 .033],...
                'Color','w', 'String',num2str(xyz(2)),'fontsize',18,'fontweight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
            an(1,3) = annotation(gcf,'textbox',[sum(wid([2 3]))+(.04) (.1/n)+(.8/n) .02 .033],...
                'Color','w', 'String',num2str(xyz(1)),'fontsize',18,'fontweight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end

    function JointPlot(varargin)
        if length(Obj)>3
            disp('This option only works if there are two and only two overlays');
            return
        end
        
%         Obj(2).col = 16;
%         set(con(21,1),'Value',2);
%         switchObj({2});
%         UpdateThreshold;
%         
%         Obj(3).col = 14;
%         set(con(21,1),'Value',3);
%         switchObj({3});
%         UpdateThreshold;
        %%%
        tmp = Obj(2).I.*double(Obj(2).mask);
        if numel(Obj(2).Thresh)==2
            ind = find(tmp<Obj(2).Thresh(1) | tmp>Obj(2).Thresh(2));
        else
            ind = find(tmp<Obj(2).Thresh(1) | tmp>Obj(2).Thresh(4) | (tmp>Obj(2).Thresh(2) & tmp<Obj(2).Thresh(3)));
        end
        tmp(ind)=NaN;
        ind1 = find(~isnan(tmp));
        %%%
        tmp = Obj(3).I.*double(Obj(3).mask);
        if numel(Obj(3).Thresh)==2
            ind = find(tmp<Obj(3).Thresh(1) | tmp>Obj(3).Thresh(2));
        else
            ind = find(tmp<Obj(3).Thresh(1) | tmp>Obj(3).Thresh(4) | (tmp>Obj(3).Thresh(2) & tmp<Obj(3).Thresh(3)));
        end
        tmp(ind)=NaN;
        ind2 = find(~isnan(tmp));
        both = intersect(ind1,ind2);
        
        Obj(4) = Obj(3);
        Obj(4).Name = 'Overlap.nii';
        Obj(4).I(:) = 0;
        Obj(4).I(both)=1;
        Obj(4).Thresh = [.5 1];
        Obj(4).clim = [.5 1.5];
        Obj(4).col = 15;
        Obj(4).mask = ones(size(Obj(4).I),'uint8');
        
        set(con(21,1),'String', [get(con(21,1),'String'); Obj(4).Name],'Value',4);

        setupFrames(4,1)
        
        drawFresh(ax1,1,4);
        drawFresh(ax2,2,4);
        drawFresh(ax3,3,4);
        
        uistack(ch(1,1),'top'); uistack(ch(1,2),'top');
        uistack(ch(2,1),'top'); uistack(ch(2,2),'top');
        uistack(ch(3,1),'top'); uistack(ch(3,2),'top');
        
        for ii = 1:length(hand);
            set(hand{ii}, 'uicontextmenu',hcmenu);
        end
        
        switchObj;
    end

    function gotoMinMax(varargin)
        movego = 0;
        vn = get(con(21,1),'Value');
        if vn==1
            return
        end
        
        if ~isempty(Obj(vn).mask);
            tmpI = Obj(vn).I.*double(Obj(vn).mask);
        else
            tmpI = Obj(vn).I;
        end
        
        if varargin{1} == item(1)
            currloc = Obj(vn).point(1:3);
            
            neigh = neighbors(currloc,Obj(vn).h.dim);
            vi = sub2ind(Obj(vn).h.dim,neigh(:,1),neigh(:,2),neigh(:,3));
            whvi = find(tmpI(vi)==max(tmpI(vi)));
            whvi = vi(whvi(1));
            currvi = vi(15);
            
            while currvi ~= whvi
                currvi = whvi;
                [currloc(1) currloc(2) currloc(3)] = ind2sub(Obj(vn).h.dim,whvi);
                neigh = neighbors(currloc,Obj(vn).h.dim);
                vi = sub2ind(Obj(vn).h.dim,neigh(:,1),neigh(:,2),neigh(:,3));
                whvi = find(tmpI(vi)==max(tmpI(vi)));
                whvi = vi(whvi(1));
            end
            
            [currloc(1) currloc(2) currloc(3)] = ind2sub(Obj(vn).h.dim,whvi);
            newloc = round([currloc 1]*Obj(vn).h.mat');
            set(con(15),'String', num2str(newloc(1:3)));
            goTo(con(15,1)); shg
        elseif varargin{1} == item(2)
            currloc = Obj(vn).point(1:3);
            
            neigh = neighbors(currloc,Obj(vn).h.dim);
            vi = sub2ind(Obj(vn).h.dim,neigh(:,1),neigh(:,2),neigh(:,3));
            whvi = find(tmpI(vi)==min(tmpI(vi)));
            whvi = vi(whvi(1));
            currvi = vi(15);
            
            while currvi ~= whvi
                currvi = whvi;
                [currloc(1) currloc(2) currloc(3)] = ind2sub(Obj(vn).h.dim,whvi);
                neigh = neighbors(currloc,Obj(vn).h.dim);
                vi = sub2ind(Obj(vn).h.dim,neigh(:,1),neigh(:,2),neigh(:,3));
                whvi = find(tmpI(vi)==min(tmpI(vi)));
                whvi = vi(whvi(1));
            end
            
            [currloc(1) currloc(2) currloc(3)] = ind2sub(Obj(vn).h.dim,whvi);
           
            newloc = round([currloc 1]*Obj(vn).h.mat');
            set(con(15),'String', num2str(newloc(1:3)));
            goTo(con(15,1)); shg
        elseif varargin{1} == item(3)
            %keyboard;
            i1 = find(tmpI==max(tmpI(:)));
            [x y z] = ind2sub(Obj(vn).h.dim,i1);
            if numel(x)>1
                disp('There are multiple maxima!');
                x = x(1); y = y(1); z = z(1);
            end
            
            newloc = round([x y z 1]*Obj(vn).h.mat');
            set(con(15),'String', num2str(newloc(1:3)));
            goTo(con(15,1)); shg
        elseif varargin{1} == item(4)
            i1 = find(tmpI==min(tmpI(:)));
            [x y z] = ind2sub(Obj(vn).h.dim,i1);
            if numel(x)>1
                disp('There are multiple minima!');
                x = x(1); y = y(1); z = z(1);
            end
            
            newloc = round([x y z 1]*Obj(vn).h.mat');
            set(con(15),'String', num2str(newloc(1:3)));
            goTo(con(15,1)); shg
        end
        movego = 0;
    end

    function scrollMove(varargin)
        if isempty(gco);
            return
        end
        
        switch gca%get(gco,'Parent')
            case Obj(1).ax1
                switch (varargin{2}.VerticalScrollCount>0)
                    case 1
                        cl = str2num(get(con(15,1),'String'));
                        cl(1) = cl(1)+1;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    otherwise
                        cl = str2num(get(con(15,1),'String'));
                        cl(1) = cl(1)-1;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                end
            case Obj(1).ax2
                switch (varargin{2}.VerticalScrollCount>0)
                    case 1
                        cl = str2num(get(con(15,1),'String'));
                        cl(2) = cl(2)+1;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    otherwise
                        cl = str2num(get(con(15,1),'String'));
                        cl(2) = cl(2)-1;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                end
            case Obj(1).ax3
                switch (varargin{2}.VerticalScrollCount>0)
                    case 1
                        cl = str2num(get(con(15,1),'String'));
                        cl(3) = cl(3)+1;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    otherwise
                        cl = str2num(get(con(15,1),'String'));
                        cl(3) = cl(3)-1;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                end
            otherwise
                return
        end 
    end

    function keyMove(varargin)
        if searchCellStr('alt',  {varargin{2}.Key})==1
            clickFlag = 1;
        end
        
        if strcmpi(varargin{2}.Modifier,'shift');
            mult = 5;
        else
            mult = 1;
        end
        
        if char(get(gcf,'currentcharacter')) == 't'
            vn = get(con(21,1),'Value');
            state = get(hand{1}(vn),'visible');
            
            if strcmpi(state,'on')
                set(hand{1}(vn),'visible','off');
                set(hand{2}(vn),'visible','off');
                set(hand{3}(vn),'visible','off');
            else
                set(hand{1}(vn),'visible','on');
                set(hand{2}(vn),'visible','on');
                set(hand{3}(vn),'visible','on');
            end
            return
        end

        if ~isempty(gco)
            return;
        end

        switch gca 
            case Obj(1).ax1
                switch char(get(gcf,'currentcharacter'))
                    case char(28)
                        cl = str2num(get(con(15,1),'String'));
                        cl(2) = cl(2)+mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(29)
                        cl = str2num(get(con(15,1),'String'));
                        cl(2) = cl(2)-mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(30)
                        cl = str2num(get(con(15,1),'String'));
                        cl(3) = cl(3)+mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(31)
                        cl = str2num(get(con(15,1),'String'));
                        cl(3) = cl(3)-mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    otherwise
                        return
                end
            case Obj(1).ax2
                switch char(get(gcf,'currentcharacter'))
                    case char(28)
                        cl = str2num(get(con(15,1),'String'));
                        cl(1) = cl(1)-mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(29)
                        cl = str2num(get(con(15,1),'String'));
                        cl(1) = cl(1)+mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(30)
                        cl = str2num(get(con(15,1),'String'));
                        cl(3) = cl(3)+mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(31)
                        cl = str2num(get(con(15,1),'String'));
                        cl(3) = cl(3)-mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    otherwise
                        return
                end
            case Obj(1).ax3
                switch char(get(gcf,'currentcharacter'))
                    case char(28)
                        cl = str2num(get(con(15,1),'String'));
                        cl(1) = cl(1)-mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(29)
                        cl = str2num(get(con(15,1),'String'));
                        cl(1) = cl(1)+mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(30)
                        cl = str2num(get(con(15,1),'String'));
                        cl(2) = cl(2)+mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    case char(31)
                        cl = str2num(get(con(15,1),'String'));
                        cl(2) = cl(2)-mult;
                        set(con(15,1),'string',num2str(cl))
                        goTo(con(15,1))
                    otherwise
                        return
                end
            otherwise
                return
        end
    end

    function groupCheck(varargin)        
        if varargin{1}==paramenu3(19)
            state = get(varargin{1},'Checked');
            if strcmpi(state,'on')
                set(varargin{1},'Checked','off');
            else
                set(varargin{1},'Checked','on');
            end
            return
        end
        
        swh1 = findobj(get(varargin{1},'Parent'),'Checked','on');
        swh2 = findobj(get(varargin{1},'Parent'),'Checked','off');
        
        if swh1 == varargin{1}
            if strcmpi(get(varargin{1},'Checked'),'on')
                set(varargin{1},'Checked','off');
            else
                set(varargin{1},'Checked','on');
            end
        else
            set(swh1,'Checked','off');
            set(varargin{1},'Checked','on');
        end
        
        %if get(varargin{1},'Parent') == paramenu3(8)
        %    get(varargin{1},'Label')
        %    keyboard;
        %    RNH = spm_vol([which('aal_MNI_V4.img')]);
        %    [RNI Rxyz] = spm_read_vols(RNH);
        %    [RNI RNH] = openIMG([which('aal_MNI_V4.img')]);
        %    RNames = load('aal_MNI_V4_List.mat');
        %end
        
        
    end

    function plotVOI(varargin)
        movego = 0;
        vn = get(con(21,1),'Value');
        
        if isempty(contrasts) || isempty(DataHeaders);
            try
                loadData
            catch
                error('No go, no support file found.');
            end
        else
            try
                if ~strcmpi(fileparts(Obj(vn).FullPath),Des.OutputDir)
                    try
                        loadData
                    catch
                        error('No go, no support file found.');
                    end
                end
            catch
                if ~strcmpi(fileparts(Obj(vn).FullPath),Des.swd)
                    try
                        loadData
                    catch
                        error('No go, no support file found.');
                    end
                end
            end
        end
        
        if ~isempty(Obj(vn).mask);
            tmpI = Obj(vn).I.*double(Obj(vn).mask);
        else
            tmpI = Obj(vn).I;
        end
        
        if strcmpi(modelType,'SPM')
            tmp = Obj(vn).DispName;
            i1 = find(tmp=='_');
            i2 = find(tmp=='.');
            tmp = tmp(i1(end)+1:i2(end)-1);
            spmCon = contrasts(str2num(tmp));
        end

        switch varargin{1}
            case item(5)  %% Cluster
                CM = str2num(get(findobj(paramenu3(5),'Checked','On'),'Label'));
                lastItem = item(5);
                adir = [fileparts(Obj(vn).FullPath) filesep];
                
                %%% Using the diplayed image get the matrix vector indices
                %%% for the selected cluster
                mni = str2num(get(con(15),'string'));
                
                thresh = Obj(vn).Thresh;
                if numel(thresh)==2
                    ind = find(tmpI<=thresh(2) & tmpI>=thresh(1));
                    L = [];
                    [L(:,1) L(:,2) L(:,3)] = ind2sub(Obj(vn).h.dim,ind);
                    A = spm_clusters2(L(:,1:3)',CM);
                elseif numel(thresh)==4
                    ind1 = find(tmpI<=thresh(2) & tmpI>=thresh(1));
                    L1 = [];
                    [L1(:,1) L1(:,2) L1(:,3)] = ind2sub(Obj(vn).h.dim,ind1);
                    A1 = spm_clusters2(L1(:,1:3)',CM);
                    
                    ind2 = find(tmpI<=thresh(4) & tmpI>=thresh(3));
                    L2 = [];
                    [L2(:,1) L2(:,2) L2(:,3)] = ind2sub(Obj(vn).h.dim,ind2);
                    A2 = spm_clusters2(L2(:,1:3)',CM);
                    
                    L = [L1; L2];
                    A = [A1 A2+max(A1)];
                    ind = [ind1; ind2];
                end
                
                L(:,4) = 1;
                L2 = L*Obj(vn).h.mat';
                dist = sqrt(sum((L2(:,1:3)-repmat(mni, size(L,1),1)).^2,2));
                if min(dist)>5
                    disp('No Cluster is selected');
                end
                ind1 = find(dist==min(dist));
                ind2 = find(A==A(ind1(1)));
                ind3 = ind(ind2);
                
                %%% Read in an orignal volume
                try
                    th = spm_vol(strtrim(Flist(1,:))); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                catch
                    th = spm_vol(Obj(vn).FullPath); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                end
                th = th(1);
                
                %try
                matLoc = [mni 1]*inv(th(1).mat'); matLoc = matLoc(1:3);
                %catch
                %    keyboard; 
                %end
                
                %%% Covert the display index into the native volume index
                Q = [];
                [Q(:,1) Q(:,2) Q(:,3)] = ind2sub(Obj(vn).h.dim,ind3);
                Q(:,4) = 1; 
                Q = (Q*Obj(vn).h.mat')*inv(th.mat');
                
                Q =  unique(round(Q),'rows');
                i1 = unique([find(Q(:,1)<=0 | Q(:,1)>th.dim(1)); find(Q(:,2)<=0 | Q(:,2)>th.dim(2)); find(Q(:,3)<=0 | Q(:,3)>th.dim(3))]);
                i2 = setdiff(1:size(Q,1),i1);
                Q = Q(i2,:);
                voxInd = sub2ind(th.dim,Q(:,1),Q(:,2),Q(:,3));

                ind = find(Obj(vn).Name == '.');
                indfilesep = find(Obj(vn).Name == filesep); %DGM Added
                wind=find(ind>indfilesep); %DGM Added
                ind=ind(wind);%DGM Added
                wh = str2num(Obj(vn).Name(ind-4:ind-1));
                
                if isempty(wh);
                    ind = find(Obj(vn).Name==filesep);
                    if isempty(ind); ind=0; end
                    wh = str2num(Obj(vn).Name(ind+1:ind+4));
                end
                
                
                VOI = [];
                VOI.mniLoc = mni;
                VOI.matLoc = matLoc;
                VOI.index = voxInd;
                VOI.isCluster = 1;
                
                try
                    [x y z] = ind2sub(th.dim,voxInd);
                    dat = zeros(numel(DataHeaders),numel(x));
                    for zz = 1:numel(DataHeaders);
                        dat(zz,:) = spm_sample_vol(DataHeaders(zz),x,y,z,0);
                    end
                    
                    if isfield(Des,'ZeroDrop')
                        if Des.ZeroDrop == 1
                            dat(dat==0)=NaN;
                        end
                    end
                    
                    if isfield(Des,'OL')
                        ti1 = intersect(voxInd,Des.OL{2}(:,2));
                        for zz = 1:numel(ti1);
                            ti2 = find(Des.OL{2}(:,2)==ti1(zz));
                            ti3 = find(voxInd == ti1(zz));
                            dat(ti3,Des.OL{2}(ti2,1))=NaN;
                        end 
                    end
                    
                    VOI.allData = dat;
                    VOI.data = nanmean(dat,2);
                catch
                    VOI.allData = [];
                end
            case item(6)  %% Sphere
                lastItem = item(6);
                adir = Obj(vn).FullPath;
                ind = find(adir==filesep);
                if isempty(ind);
                    adir = [];
                else
                    adir = adir(1:ind(end));
                end
                
                try
                    th = spm_vol(strtrim(Flist(1,:))); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                catch
                    th = spm_vol(Obj(vn).FullPath); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                end
                
                
                mniLoc = str2num(get(con(15),'string'));
                matLoc = [mniLoc 1]*inv(th(1).mat'); matLoc = matLoc(1:3);
                
                rad = get(findobj(paramenu3(2),'Checked','on'),'Label');
                rad = str2num(rad(1:end-2));
                
                voxInd = find(sqrt(sum((XYZ'-repmat(mniLoc,size(XYZ,2),1)).^2,2))<rad);
                if isempty(ind);
                    voxInd = find(sqrt(sum((XYZ'-repmat(mniLoc,size(XYZ,2),1)).^2,2))==min(sqrt(sum((XYZ'-repmat(mniLoc,size(XYZ,2),1)).^2,2))));
                end
                
                ind = find(Obj(vn).Name == '.');
                indfilesep = find(Obj(vn).Name == filesep); %DGM Added
                wind=find(ind>indfilesep); %DGM Added
                ind=ind(wind);%DGM Added
                wh = str2num(Obj(vn).Name(ind-4:ind-1));
                if isempty(wh);
                    ind = find(Obj(vn).Name==filesep);
                    if isempty(ind); ind=0; end
                    wh = str2num(Obj(vn).Name(ind+1:ind+4));
                end
                
                VOI = [];
                VOI.mniLoc = mniLoc;
                VOI.matLoc = matLoc;
                VOI.index = voxInd;
                VOI.rad = rad;
                VOI.isCluster = 0;
                
                try
                    [x y z] = ind2sub(th.dim,voxInd);
                    dat = zeros(numel(DataHeaders),numel(x));
                    for zz = 1:numel(DataHeaders);
                        dat(zz,:) = spm_sample_vol(DataHeaders(zz),x,y,z,0);
                    end
                    
                    if isfield(Des,'ZeroDrop')
                        if Des.ZeroDrop == 1
                            dat(dat==0)=NaN;
                        end
                    end
                    
                    if isfield(Des,'OL')
                        ti1 = intersect(voxInd,Des.OL{2}(:,2));
                        for zz = 1:numel(ti1);
                            ti2 = find(Des.OL{2}(:,2)==ti1(zz));
                            ti3 = find(voxInd == ti1(zz));
                            dat(ti3,Des.OL{2}(ti2,1))=NaN;
                        end 
                    end
                    
                    VOI.allData = dat;
                    VOI.data = nanmean(dat,2);
                catch
                    VOI.allData = [];
                end
            case item(7) %% Voxel
                lastItem = item(7);
                adir = Obj(vn).FullPath;
                ind = find(adir==filesep);
                if isempty(ind);
                    adir = [];
                else
                    adir = adir(1:ind(end));
                end
                
                try
                    th = spm_vol(strtrim(Flist(1,:))); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                catch
                    th = spm_vol(Obj(vn).FullPath); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                end
                
                mniLoc = str2num(get(con(15),'string'));
                matLoc = [mniLoc 1]*inv(th(1).mat'); matLoc = matLoc(1:3);
                voxInd = find(sqrt(sum((XYZ'-repmat(mniLoc,size(XYZ,2),1)).^2,2))==min(sqrt(sum((XYZ'-repmat(mniLoc,size(XYZ,2),1)).^2,2))));
                
                ind = find(Obj(vn).Name == '.');
                indfilesep = find(Obj(vn).Name == filesep); %DGM Added
                wind=find(ind>indfilesep); %DGM Added
                ind=ind(wind);%DGM Added
                wh = str2num(Obj(vn).Name(ind-4:ind-1));
                if isempty(wh);
                    ind = find(Obj(vn).Name==filesep);
                    if isempty(ind); ind=0; end
                    wh = str2num(Obj(vn).Name(ind+1:ind+4));
                end
               
                matLoc = round(matLoc);
                VOI = [];
                VOI.mniLoc = mniLoc;
                VOI.matLoc = matLoc;
                VOI.index = voxInd;
                VOI.rad = 0;
                VOI.isCluster = 0;
                 
                try
                    [x y z] = ind2sub(th.dim,voxInd);
                    dat = zeros(numel(DataHeaders),numel(x));
                    for zz = 1:numel(DataHeaders);
                        dat(zz,:) = spm_sample_vol(DataHeaders(zz),x,y,z,0);
                    end
                    
                    if isfield(Des,'ZeroDrop')
                        if Des.ZeroDrop == 1
                            dat(dat==0)=NaN;
                        end
                    end
                    
                    if isfield(Des,'OL')
                        ti1 = intersect(voxInd,Des.OL{2}(:,2));
                        for zz = 1:numel(ti1);
                            ti2 = find(Des.OL{2}(:,2)==ti1(zz));
                            ti3 = find(voxInd == ti1(zz));
                            dat(ti3,Des.OL{2}(ti2,1))=NaN;
                        end 
                    end
                    
                    VOI.allData = dat;
                    VOI.data = nanmean(dat,2);
                catch
                    VOI.allData = [];
                end
            case item(8)  %% Cached Location
                
                if ~isfield(CachedClusterLoc,'index');
                    disp('No Cluster has been Cached');
                    return
                end
                if isempty(CachedClusterLoc.index)
                    disp('No Cluster has been Cached');
                    return
                end
                tmp = CachedClusterLoc.mat-Obj(vn).h.mat;
                if sum(abs(tmp))>1e-10
                    disp('The Current Image is a different size than the one used to Cache the cluster.  One image must be resized to the other for this to work');
                    return; 
                end
                    
                try
                    th = spm_vol(strtrim(Flist(1,:))); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                catch
                    th = spm_vol(Obj(vn).FullPath); [trash,XYZ] = spm_read_vols(th(1)); clear trash
                end
                
                ind3 = CachedClusterLoc.index;
                %mni = [];
                %matLoc = []    ;
                
                %%% Covert the display index into the native volume index
                Q = [];
                [Q(:,1) Q(:,2) Q(:,3)] = ind2sub(Obj(vn).h.dim,ind3);
                Q(:,4) = 1; 
                Q = (Q*Obj(vn).h.mat')*inv(th.mat');
                
                Q =  unique(round(Q),'rows');
                i1 = unique([find(Q(:,1)<=0 | Q(:,1)>th.dim(1)); find(Q(:,2)<=0 | Q(:,2)>th.dim(2)); find(Q(:,3)<=0 | Q(:,3)>th.dim(3))]);
                i2 = setdiff(1:size(Q,1),i1);
                Q = Q(i2,:);
                voxInd = sub2ind(th.dim,Q(:,1),Q(:,2),Q(:,3));
                

                ind = find(Obj(vn).Name == '.');
                indfilesep = find(Obj(vn).Name == filesep); %DGM Added
                wind=find(ind>indfilesep); %DGM Added
                ind=ind(wind);%DGM Added
                wh = str2num(Obj(vn).Name(ind-4:ind-1));
                
                if isempty(wh);
                    ind = find(Obj(vn).Name==filesep);
                    if isempty(ind); ind=0; end
                    wh = str2num(Obj(vn).Name(ind+1:ind+4));
                end
                
                VOI = [];
                VOI.mniLoc = 'From Cached Cluster';
                VOI.matLoc = [];
                VOI.index = voxInd;
                VOI.isCluster = 1;
                
                try
                    [x y z] = ind2sub(th.dim,voxInd);
                    dat = zeros(numel(DataHeaders),numel(x));
                    for zz = 1:numel(DataHeaders);
                        dat(zz,:) = spm_sample_vol(DataHeaders(zz),x,y,z,0);
                    end
                    
                    if isfield(Des,'ZeroDrop')
                        if Des.ZeroDrop == 1
                            dat(dat==0)=NaN;
                        end
                    end
                    
                    if isfield(Des,'OL')
                        ti1 = intersect(voxInd,Des.OL{2}(:,2));
                        for zz = 1:numel(ti1);
                            ti2 = find(Des.OL{2}(:,2)==ti1(zz));
                            ti3 = find(voxInd == ti1(zz));
                            dat(ti3,Des.OL{2}(ti2,1))=NaN;
                        end 
                    end
                    
                    VOI.allData = dat;
                    VOI.data = nanmean(dat,2);
                catch
                    VOI.allData = [];
                end
            otherwise
        end
        
        if ~plotGo
            assignin('base','VOI',VOI);
            warning('Plotting options are not available for this contrast');
            return;
        end
        
        switch modelType
            case 'spm'
                c = spmCon.c;
                VOI.con = c;
                VOI.DM = DM;
                VOI.SPM = Des;
                assignin('base','VOI',VOI);
                GLM_Plot_spm(VOI,777+gcf);
            case 'regular'
                if isfield(contrasts,'c') && ~isempty(contrasts(wh).c)
                    c = contrasts(wh).c;
                else
                    try
                        c = getContrastMat(contrasts(wh));
                    catch
                        assignin('base','VOI',VOI);
                        return
                    end
                end
                 VOI.con = c;
                 VOI.DM = DM;
                 VOI.F = Des.F;
                 VOI.ConSpec = Des.Cons(wh);
                 assignin('base','VOI',VOI);
                 [yres xres] = GLM_Plot(VOI,777+gcf);
                 VOI.yres = yres;
                 VOI.xres = xres;
                 assignin('base','VOI',VOI);
            case 'fast'
                [pr1 pr2 pr3] = fileparts(Obj(vn).FullPath);
                ind = find(pr2=='_');
                effect = pr2(ind(2)+1:end);
                %figure(555+gcf); clf;
                figure; clf;
                [yres xres part] = GLM_Plot_Fast(VOI.data,DM,effect);
                VOI.effect = effect;
                VOI.mod = DM;
                VOI.part = part;
                VOI.yres = yres;
                VOI.xres = xres;
                assignin('base','VOI',VOI);
            otherwise
        end
  
        movego = 0;
    end

    function loadData(varargin)
        vn = get(con(21,1),'Value');
        
        adir = Obj(vn).FullPath;
        ind = find(adir==filesep);
        if isempty(ind);
            adir = [];
        else
            adir = adir(1:ind(end));
        end
        
        if exist([adir 'I.mat'])>0;
            %if ~isempty(searchCellStr('aschultz',{UserTime})) keyboard; end
            HH = load([adir 'I.mat']);
            if exist([adir 'FinalDataSet.nii'])
                DataHeaders = spm_vol([adir 'FinalDataSet.nii']);
                disp('Loading FinalDataSet.nii');
                tmp = [(DataHeaders.n)]; 
                Flist = [char(DataHeaders.fname) repmat(',',numel(DataHeaders),1) char(num2str(tmp(1:2:end)'))];                
                plotGo = 1;
            elseif isfield(HH.I,'Scans');
                disp('Reading in original input files');
                Flist = char(HH.I.Scans);
                DataHeaders = spm_vol(Flist);
                plotGo = 1;
            else
                plotGo = 0;
            end
            if isfield(HH.I,'OL');
                Outliers = HH.I.OL;
            end
            if isfield(HH.I,'Cons')
                contrasts = HH.I.Cons;
                DM = HH.I.F.XX;
                modelType = 'regular';
            else
                modelType = 'fast';
                DM = HH.I.MOD;
            end
            Des = HH.I;
            nvox = prod(HH.I.v.dim);
            ns = size(DM,1);
        elseif exist([adir 'SPM.mat'])>0;
            modelType = 'spm';
            HH = load([adir 'SPM.mat']);
            contrasts = HH.SPM.xCon;
            tmp = spm_read_vols(HH.SPM.xY.VY);
            ss = size(tmp);
            Des = HH.SPM;
            origDat = double(reshape(tmp, prod(ss(1:3)),ss(4))');
            DM = HH.SPM.xX.X;
            Flist = char(HH.SPM.xY.P);
            DataHeaders = HH.SPM.xY.VY;
            plotGo = 1;
        else
            disp('No stat support files were found');
            return;
        end
    end

    function applyToAll(varargin)
        vn = get(con(21,1),'Value');
        sett = setdiff(2:length(Obj),vn);
        for ii = sett
            %Obj(ii).Thresh = Obj(vn).Thresh;
            Obj(ii).clim = Obj(vn).clim;
            Obj(ii).PVal = Obj(vn).PVal;
            %Obj(ii).DF = Obj(vn).DF;
            Obj(ii).Trans = Obj(vn).Trans;
            
            set(con(21,1),'Value',ii);
            %AutoUpdate;
            switchObj;
            %UpdateThreshold;
        end
        set(con(21,1),'Value',vn);
        switchObj;
    end

    function getClusterParams
        vn = get(con(21,1),'Value');
        if vn<2
            return
        end
        
        clear S;
        
        S.mask = [];
        S.SPM = 0;
        
        if strcmpi(get(paramenu4(1),'Checked'),'on')
            S.sign = 'pos';
        elseif strcmpi(get(paramenu4(2),'Checked'),'on')
            S.sign = 'neg';
        end;
        
        if numel(Obj(vn).DF)==1;
            if any(isnan(Obj(vn).DF))
                S.type = 'none';
                S.df1 = [];
            elseif any(isinf(Obj(vn).DF))
                S.type = 'Z';
                S.df1 = inf;
            else
                S.type = 'T';
                S.df1 = Obj(vn).DF;
            end
            
        elseif numel(Obj(vn).DF)==2;
            S.type = 'F';
            S.df1 = Obj(vn).DF(1);
            S.df2 = Obj(vn).DF(2);
        end

        if strcmpi(S.type,'none');
            if numel(Obj(vn).Thresh)==4;
                error('Threshold can only contain two terms (one-sided)');
            end
            i1 = find(abs(Obj(vn).Thresh)==min(abs(Obj(vn).Thresh)));
            S.thresh = Obj(vn).Thresh(i1);
            if strcmpi(S.sign,'pos'); l = 1; else l = -1;end
            if sign(S.thresh)~=l
                error('Direction of peaks (from Parameters-->Sign) does not match sign of the specified threshold.');
            end
        else
            S.thresh = Obj(vn).PVal;
        end
        
        S.voxlimit = str2num(get(findobj(paramenu3(3),'Checked','on'),'Label'));
           
        S.separation = str2num(strtok(get(findobj(paramenu3(4),'Checked','on'),'Label'),'m'));
        S.conn = str2num(get(findobj(paramenu3(5),'Checked','on'),'Label'));

        tmp = str2num(get(con(23),'String'));
        if isempty(tmp); tmp = 0; end
        S.cluster = tmp;
        
        if strcmpi(get(findobj(paramenu3(7),'Label','Use Nearest Label'),'Checked'),'on')
            S.nearest = 1;
        else
            S.nearest = 0;
        end
        
        S.label = get(findobj(paramenu3(8),'Checked','on'),'Label');
        
        %%% New Params
        % S.UID = [];
        % S.out = [];
        % S.df2 = ?  
        
        th = Obj(vn).h;
        [a b c] = fileparts(th.fname);
        
        %%% Create a mask based on the current masking field.
        % th.pinfo = [1 0 352]';
        % th.dt = [2 0];
        % th.fname = [a filesep 'tmpPeakMask.nii'];
        % spm_write_vol(th,Obj(vn).mask);
        % S.mask = th.fname;
        
        S.exact = 0;
        % S.sphere = 0;
        % S.clustersphere=0;
        S.SV=0;

        A = load([a filesep 'I.mat']);
        % S.RESELS = [];
        
        %%% Will need to double check this with a repeated measures design.
        tmp = regexp(b,'_','split');
        effect = tmp{3};
        et = [];
        for ii = 1:numel(A.I.MOD.RFMs)
            for jj = 1:numel(A.I.MOD.RFMs(ii).Effect)
                if strcmpi(A.I.MOD.RFMs(ii).Effect(jj).name,effect)
                    et = ii;
                    break
                end
            end
        end
        
        try
            if ~isempty(et)
                S.FWHM = A.I.FWHM{et};
                % S.RESELS = spm_resels_vol(Obj(vn).h,S.FWHM)';
            end
        catch
        end
        
        S.threshc = str2num(get(findobj(paramenu3(18),'Checked','on'),'Label'));
        % S.exactvoxel = [];
        % S.FIVE = [];`
        % S.savecorrected.do = 1;
        % S.savecorrected.type = {'cFWE' 'vFWE' 'vFDR' 'cFDR'};
     end

    function getPeakInfo(varargin)
        vn = get(con(21,1),'Value');
        S = [];

        getClusterParams;
        peak = [];   
        % save S.mat S;
        [peak.voxels, peak.voxelstats, peak.clusterstats, peak.sigthresh, peak.regions, peak.mapparameters, peak.UID] = peak_nii(Obj(vn).FullPath(1:end-2),S);

        assignin('base','peak',peak);

        figure(666); clf; %reset(666);
        
        %%% Add in some sorting options for the table
        [a b] = sortrows([cellstr(char(num2str(peak.voxels{1}(:,end)))) peak.regions(:,2)]);
        tabHand = uitable('Parent',666,...
            'ColumnName',{'Cluster Size' 'T/F-Stat' 'X' 'Y' 'Z' 'N_peaks' 'Cluster Num' 'Region Num' 'Region Name'},...
            'data', [mat2cell(peak.voxels{1}(b,:),ones(size(peak.voxels{1},1),1), ones(size(peak.voxels{1},2),1)) peak.regions(b,:)],...
            'Units','Normalized','Position', [0 0 1 1],...
            'ColumnWidth', 'auto', 'RearrangeableColumns','on','CellSelectionCallback',@goToCluster);
        %set(tabHand, 'uicontextmenu',tableMenu);
        
        Out = cell(size(peak.voxels{1},1)+1 ,9);
        Out(1,:) = {'Cluster Size' 'T/F-Stat' 'X' 'Y' 'Z' 'N_peaks' 'Cluster Num'  'Region Num' 'Region Name'};
        % 'Other1' 'Other2' 'Other3' 'Other4'
        
        Out(2:end,1:7) = num2cell(peak.voxels{1}(:,1:7));
        
        Out(2:end,8:9) = peak.regions;
        
        [a b c] = fileparts(Obj(vn).FullPath);
        WriteDataToText(Out,[a filesep 'peakinfo.csv'],'w',',');
        
        %delete(th.fname);
    end

    function goToCluster(varargin)

       data = get(tabHand,'Data');
       row = varargin{2}.Indices(1);
       mni = [data{row,3:5}];
       
       set(con(15,1),'String',num2str(mni));
       goTo(con(15,1));
    end

    function out = initializeUnderlay(M,HH)
        out.Name = 'Structural Underlay';
        
        out.I = double(M);
        out.h = HH;
        out.Thresh = [min(M(:)) max(M(:))];
        if ~isempty(searchCellStr('defaultUnderlay.nii',{HH.fname}))
            out.clim = [max(M(:))*.15 max(M(:))*.8];            
        else
            out.clim = out.Thresh;
        end
        out.PVal = [];
        out.col = 1;
        out.pos = axLim(M,HH);
        out.Exclude = [];
        out.MaskInd = [];
        out.mask = ones(size(M),'uint8');
        out.Trans = 1;
        out.mask = ones(size(M),'uint8');
        out.Thresh = [min(M(:)) max(M(:))];
        out.clim = [min(M(:)) max(M(:))];
        out.Range = [out.Thresh(1) out.Thresh(2) diff(out.Thresh) ];
           
        out.axLims = size(M);
    
        mniLims = [[min(out.pos{1}) min(out.pos{2}) min(out.pos{3})]; ...
            [max(out.pos{1}) max(out.pos{2}) max(out.pos{3})]];
    end

    function changeUnderlay(varargin)
        
        ufn = spm_select(1,'image','Choose the new underlay');
        
        
        h = spm_vol(ufn);

        [I mmm] = SliceAndDice3(h,MH,[],[],[0 NaN],[]);
        h.dim = size(I);
        h.mat = mmm;
        
        nv = 1;
        for ii = 1:length(hand)
            delete(hand{ii}(nv));
        end
        
        tmp = initializeUnderlay(I,h);
        flds = fields(tmp);
        for ii = 1:length(flds)
            Obj(1).(flds{ii}) = tmp.(flds{ii});
        end
        Obj(1).point =  round([loc 1] * inv(Obj(1).h.mat)');
        setupFrames(1,1);
        
        tss = Obj(1).axLims;
        
        tmp = [1 1 1];
        tss = tss.*tmp;
        
        height = tss(2)+tss(3);
        width =  tss(1)+tss(2);
        rat = height/width;
        
        ff = get(0,'ScreenSize');
        if ff(4)>800
            pro = 3.25;
        else
            pro = 2.5;
        end
        

        ss = get(0,'ScreenSize');
        if (ss(3)/ss(4))>2
            ss(3)=ss(3)/2;
        end
        op = floor([50 ss(4)-75-((ss(3)/pro)*rat)  ss(3)/pro (ss(3)/pro)*rat]);
        
        
        set(gcf, 'Position', op);
        
        wid(1) = tss(2)/width;
        hei(1) = tss(3)/height;
        wid(2) = tss(1)/width;
        hei(2) = tss(3)/height;
        wid(3) = tss(1)/width;
        hei(3) = tss(2)/height;
          

        set(ax1,'Color','k','Position',[wid(2) hei(3) wid(1) hei(1)]); hold on;
        set(ax2,'Color','k','Position',[0      hei(3) wid(2) hei(2)]); hold on;
        set(ax3,'Color','k','Position',[0      0      wid(3) hei(3)]); hold on;
        set(ax4,'Color','w','Position',[wid(2) 0      .04     hei(3)]); axis tight;
        
        %st = wid(2)+.01+.04+.03;
        st = wid(2)+.125;
        len1 = (1-st)-.01;
        len2 = len1/2;
        len3 = len1/3;
        
        inc = (hei(3))/11;
        inc2 = .75*inc;        
        
        val = (Obj(1).Range(2)-Obj(1).Range(1))/Obj(1).Range(3);
        con(27,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[wid(2)+.070      0  .025   hei(3)*.9],'Value',val,'CallBack', @adjustUnderlay);
        val = (Obj(1).Range(1)-Obj(1).Range(1))/Obj(1).Range(3);
        con(28,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[wid(2)+.095      0  .025   hei(3)*.9],'Value',val,'CallBack', @adjustUnderlay);
        %con(29,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[wid(2)+.095      0  .025   hei(3)*.9],'Value',0,'CallBack', @adjustUnderlay);
        
        con(30,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[wid(2)+.070    hei(3)*.9   .05 .05]); % wid(2)+.0450    hei(3)*.9   .075 .05]
        
        % con(28,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[st-.1          (.01+(1*inc)+(0*inc2))  len1*.75   inc],'Value',1,'CallBack', @adjustTrans);
        
        con(1,1)   = uicontrol(pane,'style','popupmenu', 'Units','Normalized','Position',[st          .01+(0*inc)+(0*inc2)    len1   inc],'String',[{'Colormap'} cmaps(:)'],'CallBack', @changeColorMap); shg
        %con(1,1)   = uicontrol(pane,'style','popupmenu', 'Units','Normalized','Position',[st          .01+(0*inc)+(0*inc2)    len1   inc],'String',[{'Colormap' 'A' 'B' 'C' 'D'}],'CallBack', @changeColorMap); shg
        con(2,1)   = uicontrol(pane,'style','slider',    'Units','Normalized','Position',[st          (.01+(1*inc)+(0*inc2))  len1*.75   inc],'Value',1,'CallBack', @adjustTrans);
        con(23,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+len1*.75 (.01+(1*inc)+(0*inc2))  len1*.25   inc],'String','All','CallBack', @applyToAll);
        con(3,1)   = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(2*inc)+(0*inc2)  len1   inc2],'String', 'Transparency','fontsize',12);
        
        con(4,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(2*inc)+(1*inc2)  len2   inc],'CallBack',@UpdateThreshold);
        con(5,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st+len2     .01+(2*inc)+(1*inc2)  len2   inc],'CallBack',@UpdateCLims);
        con(6,1)   = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(3*inc)+(1*inc2)  len2   inc2],'String', 'Thresh','fontsize',12);
        con(7,1)   = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len2     .01+(3*inc)+(1*inc2)  len2   inc2],'String', 'Color Limits','fontsize',12);
        con(8,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(3*inc)+(2*inc2)  len3   inc],'CallBack',@UpdatePVal);
        con(9,1)   = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st+len3     .01+(3*inc)+(2*inc2)  len3   inc],'CallBack',@UpdatePVal);
        con(23,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st+len3+len3     .01+(3*inc)+(2*inc2)  len3   inc],'CallBack',@ExtentThresh);
        con(10,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(4*inc)+(2*inc2)  len3   inc2],'String', 'DF','fontsize',12);
        con(11,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len3     .01+(4*inc)+(2*inc2)  len3   inc2],'String', 'P-Value','fontsize',12);
        con(24,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len3+len3     .01+(4*inc)+(2*inc2)  len3   inc2],'String', 'Extent','fontsize',12);
        
        con(15,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(4*inc)+(3*inc2)  len2   inc],'String',num2str(loc),'CallBack',@goTo);
        %con(15,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[st          .01+(4*inc)+(3*inc2)  len2   inc],'String',['0 0 0'],'CallBack',@goTo);
        con(25,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+len2     .01+(4*inc)+(3*inc2)  len2/2   inc],'String','FDR','CallBack', @correctThresh);
        con(26,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+(len2*1.5)  .01+(4*inc)+(3*inc2)  len2/2   inc],'String','FWE','CallBack', @correctThresh);
        con(17,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st          .01+(5*inc)+(3*inc2)  len2   inc2],'String', 'MNI Coord',  'fontsize',12);
        con(18,1)  = uicontrol(pane,'style','text',      'Units','Normalized','Position',[st+len2     .01+(5*inc)+(3*inc2)  len2   inc2],'String', 'MC Correct','fontsize',12);
        con(12,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st          .01+(5*inc)+(4*inc2) len2   inc],'String',{'Move Up'},'CallBack',@changeLayer);
        con(13,1)  = uicontrol(pane,'style','pushbutton','Units','Normalized','Position',[st+len2     .01+(5*inc)+(4*inc2) len2   inc],'String',{'Move Down'},'CallBack',@changeLayer);
        con(21,1)  = uicontrol(pane,'style','popupmenu', 'Units','Normalized','Position',[st          .01+(5*inc)+(5*inc2) len1   inc],'String',{'Overlays'},'CallBack',@switchObj);
        con(22,1)  = uicontrol(pane,'style','edit',      'Units','Normalized','Position',[wid(2)-(len2/2) hei(3)-(inc*1.05) (len2/2) inc]);
        
        
        con(19,1)  = uicontrol(pane,'style','PushButton','Units','Normalized','Position',[st          .01+(6*inc)+(5*inc2) len2   inc],'String','Open Overlay','FontWeight','Bold','CallBack',@openOverlay);
        con(20,1)  = uicontrol(pane,'style','PushButton','Units','Normalized','Position',[st+len2     .01+(6*inc)+(5*inc2) len2   inc],'String','Remove Volume', 'FontWeight','Bold','CallBack',@removeVolume);

        
        drawFresh(ax1,1,1); axis equal;
        drawFresh(ax2,2,1); axis equal;
        drawFresh(ax3,3,1); axis equal;
         
        set(ax1,'color','k')
        set(ax2,'color','k')
        set(ax3,'color','k')

        axes(ax1); ax = axis;
        set(ch(1,1),'YData',ax(3:4)); set(ch(1,2),'XData',ax(1:2));
        axes(ax2); ax = axis;
        set(ch(2,1),'YData',ax(3:4)); set(ch(2,2),'XData',ax(1:2));
        axes(ax3); ax = axis;
        set(ch(3,1),'YData',ax(3:4)); set(ch(3,2),'XData',ax(1:2));
        
        uistack(hand{1}(1),'bottom');
        uistack(hand{2}(1),'bottom');
        uistack(hand{3}(1),'bottom');
        
        set(con(27),'value',1);
        set(con(28),'value',0);
        
        if get(con(21,1),'Value') == 1
            set(con(5),'String', sprintf('%0.3f %0.3f',Obj(1).Thresh(1),Obj(1).Thresh(2)));
        end
    end

    function plotConfig
       fh = figure(777); clf; reset(777);
       
       pop(1) = uicontrol(fh,'style','text','Units','Normalized','position',[.05 .95 .9 .05],'String','Choose the plot Type','Fontsize',12);
       bg = uibuttongroup(fh,'Units','Normalized','position',[.05 .9 .9 .05]);

       pop(2) = uicontrol(bg,'style','radiobutton','Units','Normalized','position',[ 00 0 .33 1],'String','BoxPlot','Fontsize',12);
       pop(3) = uicontrol(bg,'style','radiobutton','Units','Normalized','position',[.33 0 .33 1],'String','InteractionPlot','Fontsize',12);
       pop(4) = uicontrol(bg,'style','radiobutton','Units','Normalized','position',[.66 0 .33 1], 'String', 'ScatterPlot','Fontsize',12);
       
       pop(5) = uicontrol(fh,'style','text','Units','Normalized','position',[.05 .83 .9 .05],'String','Set the Data Groups','Fontsize',12);
       pop(6) = uicontrol(fh,'style','edit','Units','Normalized','position',[.05 .75 .9 .08],'String','{{1:10} {11:20} {etc..}}','Fontsize',12);
       
       pop(7) = uicontrol(fh,'style','text','Units','Normalized','position',[.05 .68 .9 .05],'String','Specify Design Matirx Columns','Fontsize',12);
       pop(8) = uicontrol(fh,'style','edit','Units','Normalized','position',[.05 .60 .9 .08],'String','{{1:10} {11:20} {etc..}}','Fontsize',12);
    end

    function ExtentThresh(varargin)
        %keyboard;
        %CM = 18;%
        CM = str2num(get(findobj(paramenu3(5),'Checked','On'),'Label'));

        vn = get(con(21,1),'Value');
        ClusterExtent = str2num(get(con(23),'String'));
        if isempty(ClusterExtent) || ClusterExtent==0;
            Obj(vn).mask = ones(size(Obj(vn).I),'uint8');
            Obj(vn).Exclude = [];
            Obj(vn).mask(Obj(vn).MaskInd)=0;
            return
        end
        
        
        Obj(vn).ClusterThresh = ClusterExtent;
        
         
        if numel(Obj(vn).Thresh)==2
            % ind = find(Obj(vn).I<Obj(vn).Thresh(1) | Obj(vn).I>Obj(vn).Thresh(2));
            ind = find(Obj(vn).I>Obj(vn).Thresh(1) & Obj(vn).I<Obj(vn).Thresh(2));
        else
            % ind = find(Obj(vn).I<Obj(vn).Thresh(1) | Obj(vn).I>Obj(vn).Thresh(4) | (Obj(vn).I>Obj(vn).Thresh(2) & Obj(vn).I<Obj(vn).Thresh(3)));
            ind = find( (Obj(vn).I>Obj(vn).Thresh(1) & Obj(vn).I<Obj(vn).Thresh(2)) | (Obj(vn).I>Obj(vn).Thresh(3) & Obj(vn).I<Obj(vn).Thresh(4)));
        end
        % ind = setdiff(1:prod(Obj(vn).axLims),ind);
        ind = ind(:)';
        
        L = [];
        [L(:,1) L(:,2) L(:,3)] = ind2sub(Obj(vn).h.dim,ind);
        A = spm_clusters2(L(:,1:3)',CM);
        
        Exclude = [];
        list = unique(A);
        for ii = list
           i1 = find(A==ii);
           if numel(i1)<ClusterExtent
               Exclude = [Exclude ind(i1)];
           end
        end

        Obj(vn).mask(:) = 1;
        Obj(vn).Exclude = Exclude;
        Obj(vn).mask(Exclude)=0;
        Obj(vn).mask(Obj(vn).MaskInd)=0;

        setupFrames(vn,1);
        updateGraphics({1:3 vn},1); 
    end

    function out = popup(Message)
        %%% Fix this up later so that multiple messages mean multiple
        %%% inputs
        fh = figure(777); clf; %reset(777);
                
        pop(1) = uicontrol(fh,'style','text','Units','Normalized','position',[.05 .83 .9 .05],'String',Message,'Fontsize',12);
        pop(2) = uicontrol(fh,'style','edit','Units','Normalized','position',[.05 .75 .9 .08],'String','','Fontsize',12,'Callback','uiresume');
        
        uiwait
        eval(['out = [' get(pop(2),'String') '];']);
        close(gcf);

    end

    function TestFunc(varargin)
        
        vn = get(con(21,1),'Value');
        wh = Count;
        Obj(wh) = Obj(1);

        tmpVol = SliceAndDice3(Obj(vn).FullPath,Obj(1).h,Obj(1).h,Obj(1).h,[1 NaN],[]);
        
        ind = find(tmpVol>Obj(vn).Thresh(1) & tmpVol<Obj(vn).Thresh(2));
        Obj(wh).I(setdiff(1:numel(Obj(wh).I),ind)) = NaN;
        tmp = Obj(1).CM(:,:,:,1); tmp = round((tmp./max(tmp(:)))*(size(cmap{1,1},1)-1));
        nn = zeros(numel(tmp),3)*NaN;
        ind2 = find(~isnan(tmp));
        nn(ind2,:) = cmap{2,1}(tmp(ind2)+1,:);
        
        a = zeros(size(Obj(1).CM))*NaN;
        b = zeros(size(Obj(1).CM))*NaN;
        c = zeros(size(Obj(1).CM))*NaN;
        a(ind) = nn(ind,1);
        b(ind) = nn(ind,2);
        c(ind) = nn(ind,3);
        Obj(wh).CM = cat(4,a,b,c);
        
        m = Obj(wh).I;
        n2 = 'TmpImg';
        
        set(con(21,1),'String', [get(con(21,1),'String'); n2],'Value',Count);
        
        set(con(4),'String',[num2str(floor(min(m(:))*1000)/1000) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
        set(con(5),'String',[num2str(floor(min(m(:))*1000)/1000) ', ' num2str(ceil(max(m(:))*1000)/1000)]);
        set(con(8),'String','NaN');
        set(con(9),'String','NaN');
        
        Obj(wh).DF = NaN;
        Obj(wh).PVal = NaN;
        Obj(wh).Thresh = [floor(min(m(:))*1000)/1000 ceil(max(m(:))*1000)/1000];
        Obj(wh).clim = [floor(min(m(:))*1000)/1000 ceil(max(m(:))*1000)/1000];
        Obj(wh).col = 2;
        Obj(wh).Trans = 1;
        
        
        y = get(ch(1,1),'XData');
        x = get(ch(2,1),'XData');
        z = get(ch(1,2),'YData');
        loc = [x(1) y(1) z(1)];
        
        Obj(wh).point = round([loc 1] * inv(Obj(wh).h.mat)');
        
        
        setupFrames(wh,0);
        Obj(wh).pos = axLim(Obj(wh).I,Obj(wh).h);
        
        tmp = size(Obj(wh).CM);
        Obj(wh).axLims = tmp(1:3);
                
        set(con(1,1),'Value',wh+1); 
        
        drawFresh(ax1,1,wh);
        drawFresh(ax2,2,wh);
        drawFresh(ax3,3,wh);
        
        uistack(ch(1,1),'top'); uistack(ch(1,2),'top');
        uistack(ch(2,1),'top'); uistack(ch(2,2),'top');
        uistack(ch(3,1),'top'); uistack(ch(3,2),'top');
        
        for ii = 1:length(hand);
            set(hand{ii}, 'uicontextmenu',hcmenu);
        end
        
        Count = Count+1;
    end

    function c = getContrastMat(Con)
        tx = Des.F.XX;
        
        if ~isfield(Con,'c') || isempty(Con.c)
            
            if iscell(Con.Groups)
                tc = [];
                for zz = 1:length(Con.Groups)
                    if size(Con.Groups{zz},1)>1
                        tmpX = Con.Groups{zz};
                    else
                        tmpX = tx(:,Con.Groups{zz});
                    end
                    cc = MakePreCons(tmpX,tx,Des.F.CovarCols);
                    tc(:,zz) = mean(cc,2);
                end
                levs = Con.Levs(end:-1:1);
                if levs == 0
                    c = tc;
                else
                    tmp = reshape(tc,[size(tc,1) levs]);
                    tmp = squeeze(tmp);
                    if size(tmp,2)==1
                        c = tmp;
                    else
                        c = differencer(tmp);
                        if numel(levs)==1
                            c = c*-1;
                        end
                    end
                end
            else
                [r,c,cc] = MakeContrastMatrix('a',Con.Groups,Con.Levs,tx);
            end
            
        else
            c = Con.c;
            c0 = eye(size(txx,2))-(c*pinv(c));
            x0 = txx*c0;
            r = eye(size(x0,1))-(x0*pinv(x0));
        end
        
    end

    function saveImg(varargin)
        %%% Write something in the descip field to describe how the image
        %%% was made
        from = varargin{1};
        vn = get(con(21,1),'Value');
        CM = str2num(get(findobj(paramenu3(5),'Checked','On'),'Label'));
        
        switch from
            case menu(19)
                disp('Save Thresholded Image');
                dotI = find(Obj(vn).FullPath=='.');
                fn = [Obj(vn).FullPath(1:dotI(end)-1) '_thresh.nii'];
                th = spm_vol(Obj(vn).FullPath);
                hh = Obj(vn).h;
                hh.fname = fn;

                if ~isempty(Obj(vn).mask);
                    tmp = Obj(vn).I.*double(Obj(vn).mask);
                else
                    tmp = Obj(vn).I;
                end
                
                if numel(Obj(vn).Thresh)==2
                    ind = find(tmp<Obj(vn).Thresh(1) | tmp>Obj(vn).Thresh(2));
                else
                    ind = find(tmp<Obj(vn).Thresh(1) | tmp>Obj(vn).Thresh(4) | (tmp>Obj(vn).Thresh(2) & tmp<Obj(vn).Thresh(3)));
                end
                tmp(ind)=NaN;

                spm_write_vol(hh,tmp);
                
                hh2 = spm_vol(fn);
                N = resizeVol2(hh2,th);
                th.fname = fn;
                th.descrip = [];
                spm_write_vol(th,N);                    
            case menu(20)
                disp('Save Masked Image');
                dotI = find(Obj(vn).FullPath=='.');
                fn = [Obj(vn).FullPath(1:dotI(end)-1) '_mask.nii'];
                th = spm_vol(Obj(vn).FullPath);
                hh = Obj(vn).h;
                hh.fname = fn;
                
                hh.dt = [2 0];
                
                tmp = Obj(vn).I;
                if numel(Obj(vn).Thresh)==2
                    ind = find(tmp>Obj(vn).Thresh(1) & tmp<Obj(vn).Thresh(2));
                else
                    ind = find((tmp>Obj(vn).Thresh(1) & tmp<Obj(vn).Thresh(2)) | (tmp>Obj(vn).Thresh(3) & tmp<Obj(vn).Thresh(4)));
                end
                tmp(:)=0;
                tmp(ind)=1;
                
                spm_write_vol(hh,tmp);
                
                hh2 = spm_vol(fn);
                N = resizeVol(hh2,th);
                th.fname = fn;
                th.descrip = [];
                th.dt = [2 0];
                spm_write_vol(th,N);
                
            case menu(21)
                disp('Save Cluster Image');

                mni = str2num(get(con(15),'string'));
                
                
                tmp = Obj(vn).I;
                if numel(Obj(vn).Thresh)==2
                    ind = find(tmp<Obj(vn).Thresh(1) | tmp>Obj(vn).Thresh(2));
                else
                    ind = find(tmp<Obj(vn).Thresh(1) | tmp>Obj(vn).Thresh(4) | (tmp>Obj(vn).Thresh(2) & tmp<Obj(vn).Thresh(3)));
                end
                tmp(ind)=NaN;
                ind = find(~isnan(tmp));
                
                [L(:,1) L(:,2) L(:,3)] = ind2sub(Obj(vn).h.dim,ind);
                A = spm_clusters2(L(:,1:3)',18);
                
                L(:,4) = 1;
                L2 = L*Obj(vn).h.mat';
                dist = sqrt(sum((L2(:,1:3)-repmat(mni, size(L,1),1)).^2,2));
                if min(dist)>5
                    disp('No Cluster is selected');
                end
                ind1 = find(dist==min(dist));
                ind2 = find(A==A(ind1(1)));
                ind3 = ind(ind2);
                
                tmp(:) =  NaN;
                tmp(ind3) = Obj(vn).I(ind3);
                
                dotI = find(Obj(vn).FullPath=='.');
                fn = [Obj(vn).FullPath(1:dotI(end)-1) '_cluster.nii'];
                th = spm_vol(Obj(vn).FullPath);
                hh = Obj(vn).h;
                hh.fname = fn;
                spm_write_vol(hh,tmp);
                
                hh2 = spm_vol(fn);
                N = resizeVol2(hh2,th);
                th.fname = fn;
                th.descrip = [];
                spm_write_vol(th,N);
                
            case menu(22)
                disp('Save Cluster Mask');
                
                mni = str2num(get(con(15),'string'));
                
                tmp = Obj(vn).I;
                if numel(Obj(vn).Thresh)==2
                    ind = find(tmp<Obj(vn).Thresh(1) | tmp>Obj(vn).Thresh(2));
                else
                    ind = find(tmp<Obj(vn).Thresh(1) | tmp>Obj(vn).Thresh(4) | (tmp>Obj(vn).Thresh(2) & tmp<Obj(vn).Thresh(3)));
                end
                tmp(ind)=NaN;
                ind = find(~isnan(tmp));

                [L(:,1) L(:,2) L(:,3)] = ind2sub(Obj(vn).h.dim,ind);
                A = spm_clusters2(L(:,1:3)',18);
                
                L(:,4) = 1;
                L2 = L*Obj(vn).h.mat';
                dist = sqrt(sum((L2(:,1:3)-repmat(mni, size(L,1),1)).^2,2));
                if min(dist)>5
                    disp('No Cluster is selected');
                end
                ind1 = find(dist==min(dist));
                ind2 = find(A==A(ind1(1)));
                ind3 = ind(ind2);
                
                tmp(:) =  0;
                tmp(ind3) = 1;
                
                dotI = find(Obj(vn).FullPath=='.');
                fn = [Obj(vn).FullPath(1:dotI(end)-1) '_clusterMask.nii'];
                th = spm_vol(Obj(vn).FullPath);
                hh = Obj(vn).h;
                hh.fname = fn;
                spm_write_vol(hh,tmp);
                
                hh2 = spm_vol(fn);
                N = resizeVol(hh2,th);
                th.fname = fn;
                th.descrip = [];
                th.dt = [2 0];
                spm_write_vol(th,N);
                
            otherwise
        end
    end

    function resampleIm(varargin)
        vn = get(con(21,1),'Value');
        %keyboard;
        voxdim = str2num(char(regexp(get(varargin{1},'Label'),'x','split')))';

        nnn = Obj(vn).FullPath;
        hh = spm_vol(nnn);
        [mm mat] = SliceAndDice3(hh, MH, voxdim, [],[1 0],[]);
        hh.dim = size(mm);
        hh.mat = mat;
        
        Obj(vn).h = hh;
        Obj(vn).I = double(mm);
        
        
        Obj(vn).mask = ones(size(mm),'uint8');
        Obj(vn).pos = axLim(Obj(vn).I,Obj(vn).h);
        
        Obj(vn).axLims = Obj(vn).h.dim;

        loc = str2num(get(con(15),'String'));
        Obj(vn).point =  round([loc 1] * inv(Obj(vn).h.mat)');
        
        setupFrames(vn,1);
        
        for ii = 1:length(hand)
            delete(hand{ii}(vn));
        end
        
        drawFresh(ax1,1,vn);
        drawFresh(ax2,2,vn);
        drawFresh(ax3,3,vn);
        
        uistack(ch(1,1),'top'); uistack(ch(1,2),'top');
        uistack(ch(2,1),'top'); uistack(ch(2,2),'top');
        uistack(ch(3,1),'top'); uistack(ch(3,2),'top');
        
        for ii = 1:length(hand);
            set(hand{ii}, 'uicontextmenu',hcmenu);
        end        
    end

    function correctThresh(varargin)
        
        vn = get(con(21,1),'Value');
        if numel(Obj(vn).DF)==1
            stat = 'T';
            df = [1 Obj(vn).DF];
        elseif numel(Obj(vn).DF)==2
            stat = 'F';
            df = Obj(vn).DF;
        end
        
        if ~isempty(Obj(vn).mask);
            tmpI = Obj(vn).I.*double(Obj(vn).mask);
        else
            tmpI = Obj(vn).I;
        end
        
        ind = find(Obj(vn).Name == '.');
            indfilesep = find(Obj(vn).Name == filesep); %DGM Added
            wind=find(ind>indfilesep); %DGM Added
            ind=ind(wind);%DGM Added
            wh = str2num(Obj(vn).Name(ind-4:ind-1));
            if isempty(wh);
                ind = find(Obj(vn).Name==filesep);
                if isempty(ind); ind=0; end
                wh = str2num(Obj(vn).Name(ind+1:ind+4));
            end
        
        if varargin{1}==con(25)  %FDR Correction
            if numel(Obj(vn).Thresh) == 2
                ind = find(abs(tmpI)>0);
                div = 1;
                if sum(sign(Obj(vn).Thresh))<0
                    %ind = find(Obj(vn).I<0);                   
                    t = abs(tmpI(ind));
                    
                elseif sum(sign(Obj(vn).Thresh))>0
                    %ind = find(Obj(vn).I>0);
                    t = abs(tmpI(ind));
                end
                
            elseif numel(Obj(vn).Thresh) == 4
                ind = find(abs(tmpI)>0);
                t = abs(tmpI(ind));
                div = 2;
            end
            
            
            if strcmpi(stat,'f');
                pp = 1-spm_Fcdf(t,df);
            else
                pp = (1-spm_Tcdf(t,df(2))).*div;
            end
            
            pp = sort(pp);
            %keyboard;
            alpha = str2num(get(findobj(paramenu3(18),'Checked','On'),'Label'));
            p = pp';
            rh = ((1:numel(p))./numel(p)).*alpha;
            below = p<=rh;
            in = (find(below==1));
            if isempty(in)
                warning('FDR Correction is not possible for this alpha');
                return
            end
            %;
            t = sort(t,'descend');
            thresh = t(in(end));
            pthresh = p(in(end));
            %[FDR_alpha FDR_p FDR_t] = computeFDR('0003_T_TFO_Pos_Corr.nii',114,.05,0)
        end
        
        if varargin{1}==con(26)  %% FWE correction using RFT
            [a b c] = fileparts(Obj(vn).FullPath);
            try
                HH = load([a filesep 'I.mat']);
                %R = HH.I.ReslInfo{HH.I.Cons(wh).ET};
                if isfield(HH.I,'Cons')
                    FWHM = HH.I.FWHM{HH.I.Cons(wh).ET};
                end
                if isfield(HH.I,'MOD')
                    ind = find(b=='_');
                    effect = b(ind(2)+1:end);
                    for ii = 1:numel(HH.I.MOD.RFMs)
                        for jj = 1:numel(HH.I.MOD.RFMs(ii).Effect)
                            if strcmpi(effect,HH.I.MOD.RFMs(ii).Effect(jj).name)
                                ET = ii;
                            end
                        end
                    end
                    FWHM = HH.I.FWHM{ET};
                end
            catch
                HH = load([a filesep 'SPM.mat']);
                %R = HH.SPM.xVol.R;
                FWHM = HH.SPM.xVol.FWHM;
            end
            
            p = str2num(get(findobj(paramenu3(18),'Checked','On'),'Label'));
            
            try
                [msk mh] = openIMG([a filesep 'mask.img']);
            catch
                try
                    [msk mh] = openIMG([a filesep 'mask.nii']);
                catch
                    [msk mh] = openIMG([a filesep 'NN.nii']);
                    mh.fname = 'mask.nii';
                    spm_write_vol(mh,msk>0);
                    [msk mh] = openIMG([a filesep 'mask.nii']);
                end
            end
            
            if ~isempty(Obj(vn).mask);
                msk = Obj(vn).I.*double(Obj(vn).mask);
            end
            
            reselinfo = spm_resels_vol(mh,FWHM)';
            
            if strcmpi(stat,'f');
                thresh = spm_uc(p,df,stat,reselinfo,1,numel(find(msk>0)));
                pthresh = 1-spm_Fcdf(thresh,df);
            else
                if numel(Obj(vn).Thresh) == 4
                    thresh = spm_uc(p/2,df,stat,reselinfo,1,numel(find(msk>0)));
                    pthresh = (1-spm_Tcdf(thresh,df(2)));
                else
                    thresh = spm_uc(p,df,stat,reselinfo,1,numel(find(msk>0)));
                    pthresh = (1-spm_Tcdf(thresh,df(2)));
                end
            end
        end
        
        thresh = round(thresh*1000)/1000;
        if numel(Obj(vn).Thresh) == 2
            if sum(sign(Obj(vn).Thresh))<0
                if thresh<Obj(vn).Thresh(1)
                    warning(['No voxels beyond the corrected threshold of ' num2str(-thresh)]);
                    return
                end
                Obj(vn).Thresh = [Obj(vn).Thresh(1) -thresh];
                Obj(vn).PVal = pthresh;
            else
                if thresh>Obj(vn).Thresh(2)
                    warning(['No voxels beyond the corrected threshold of ' num2str(thresh)]);
                    return
                end
                Obj(vn).Thresh = [thresh Obj(vn).Thresh(2)];
                Obj(vn).PVal = pthresh;
            end
        elseif numel(Obj(vn).Thresh) == 4
            up = ~(thresh>Obj(vn).Thresh(4));
            down = ~(thresh<Obj(vn).Thresh(1));
            
            if ~up && ~down
                warning(['Nothing Left at this Threshold. Exiting Correction.']);
                return;
            end
            if up && down
                Obj(vn).Thresh = [Obj(vn).Thresh(1) -thresh thresh Obj(vn).Thresh(4)];
                Obj(vn).PVal = pthresh;
            end
            if up && ~down
                warning(['No voxels below the corrected threshold of ' num2str(-thresh)]);
                Obj(vn).Thresh = [thresh Obj(vn).Thresh(4)];
                Obj(vn).PVal = pthresh;
            end
            if down && ~up
                warning(['No voxels above the corrected threshold of ' num2str(thresh)]);
                Obj(vn).Thresh = [Obj(vn).Thresh(1) -thresh];
                Obj(vn).PVal = pthresh;
            end            
        end
        
        ll = []; for zz = 1:numel(Obj(vn).Thresh); ll = [ll num2str(Obj(vn).Thresh(zz)) ' ']; end;
        ll = ll(1:end-1);
        set(con(4,1),'String',num2str(ll));

        UpdateThreshold(con(4,1));
    end

    function changeLabelMap(varargin)
        groupCheck(varargin{1})
        
%         fn1 = which([get(varargin{1},'Label') '.img']);
%         fn2 = which([get(varargin{1},'Label') '.mat']);

        try 
            [fn1,fn2]=getLabelMap(get(varargin{1},'Label'));
        catch
            fn1=[];fn2=[];
        end
        
        if isempty(fn1) && isempty(fn2)
            warning(['Label Map ' get(varargin{1},'Label')  ' cannot be found. Reverting to aal_MNI_V4']);
            set(varargin{1},'Checked','off');
            set(findobj(gcf,'Label','aal_MNI_V4'),'Checked','on');
            return
        end
        
        RNH = spm_vol(fn1);
        [RNI Rxyz] = spm_read_vols(RNH);
        RNames = load(fn2);
    end

    function h = moveYaxLabs(currAx,wh)
        if nargin==0 || isempty(currAx)
            currAx = gca;
        end
        
        if wh == 'y'
            ax = axis(currAx);
            y = get(currAx,'YTick');
            lab = cellstr(get(currAx,'YTickLabel'));
            for ii = 1:numel(lab)
                if lab{ii}(1) ~= '-'
                    lab{ii} = ['+' lab{ii}];
                end
            end
            
            adj = range(y)*.015;
            adj2 = range(ax(1:2))*.5;
             
            for ii = 1:numel(lab)
                if ii==1
                    h(ii) = text(ax(1)+adj2,y(ii)+adj,lab{ii},'FontName', get(currAx,'FontName'), 'FontSize',get(currAx,'FontSize'),'color','w','HorizontalAlignment','Center'); shg
                elseif ii==numel(lab)
                    h(ii) = text(ax(1)+adj2,y(ii)-adj,lab{ii},'FontName', get(currAx,'FontName'), 'FontSize',get(currAx,'FontSize'),'color','w','HorizontalAlignment','Center'); shg
                else
                    h(ii) = text(ax(1)+adj2,y(ii),    lab{ii},'FontName', get(currAx,'FontName'), 'FontSize',get(currAx,'FontSize'),'color','w','fontweight','bold','HorizontalAlignment','Center'); shg
                    %p = get(h(ii),'Position'); p(3) = 15; set(h(ii),'Position',p);
                end
            end
            
            set(gca,'XTick', [],'YTick',[]);
        elseif wh == 'x'
            ax = axis(currAx);
            x = get(currAx,'XTick');
            lab = cellstr(get(currAx,'XTickLabel'));
            for ii = 1:numel(lab)
                if lab{ii}(1) ~= '-'
                    lab{ii} = ['+' lab{ii}];
                end
            end
            
            adj2 = range(ax(3:4))*.5;
            
            for ii = 1:numel(lab)
                if ii==1
                    h(ii) = text(ax(1)+(.005*adj2), ax(3)+adj2, lab{ii},'FontName', get(gca,'FontName'), 'FontSize',get(gca,'FontSize'),'color','w','HorizontalAlignment','Left'); shg 
                elseif ii==numel(lab)
                    h(ii) = text(ax(2)-(.005*adj2), ax(3)+adj2, lab{ii},'FontName', get(gca,'FontName'), 'FontSize',get(gca,'FontSize'),'color','w','HorizontalAlignment','Right'); shg
                else
                    h(ii) = text(x(ii)-(.005*adj2), ax(3)+adj2, lab{ii},'FontName', get(gca,'FontName'), 'FontSize',get(gca,'FontSize'),'color','w','fontweight','bold','HorizontalAlignment','Center'); shg
                end
            end
            
            set(gca,'XTick', [],'YTick',[]);
        end
    end

    function initializeConnExplore(varargin)
        vn = get(con(21,1),'Value');
        
        if ConExp==1
            switchObj(ConLayer)
            return
        end
        
        fn = ['/autofs/space/schopenhauer_004/users/ConnectivityAtlas/Maps/blank.nii'];
        
        openOverlay(fn);
        
        ConExp = 1;
        ConHeader = spm_vol(fn);
        vn = get(con(21,1),'Value');
        ConLayer = vn;
        %Obj(2).h.mat
        
        Obj(vn).Thresh = [-1 1];
        Obj(vn).clim = [-1 1];
        
        set(con(4),'String','-1 1');
        set(con(5),'String','-1 1');
    end

    function updateConnMap(varargin)
        if ConExp==0
            return
        end
        vn = get(con(21,1),'Value');
        if  vn ~= ConLayer;
            return
        end
            
        MNI = Obj(1).point*Obj(1).h.mat';
        matLoc = round(MNI*inv(ConHeader.mat)');

        thisLoc = sub2ind(ConHeader.dim,matLoc(1),matLoc(2),matLoc(3));
        fn = ['/autofs/space/schopenhauer_004/users/ConnectivityAtlas/Maps/Vox_' sprintf('%0.6d',thisLoc) '.nii'];
        
        if exist(fn)>0
            Obj(vn).I = openIMG(fn);
            %Obj(vn).I = resizeVol2(spm_vol(fn),Obj(vn).h);
        else
            Obj(vn).I(:) = NaN;
        end
        UpdateThreshold;
    end

    function ssConn(varargin)
        vn = get(con(21,1),'Value');
        
        if ssConExp==1
            switchObj(ssConLayer)
            return
        end
        
        ff = spm_select(inf,'image');
        
        fn = [];
        for zz = 1:size(ff,1)
            fn{zz} = (ff(zz,1:end-2));
        end
        
        
        openOverlay(ff(1,:));
                
        ssConExp = 1;
        vn = get(con(21,1),'Value');
        ssConHeader = Obj(vn).h;
        ssConLayer = vn;
        

        ssData = [];
        for zz = 1:numel(fn);
            pth = fileparts(fn{zz});
%             R = [];
%             fl = [pth '/ExtraRegressors.mat'];
%             if exist(fl,'file')>0
%                 load(fl);
%                 [rows,cols] = find(R==1);
%             else
%                 rows = [];
%             end
            
            th = spm_vol(fn{zz});
            
%             wh = setdiff(1:numel(th), rows);
%             th = th(wh);
            dat = zeros(numel(th),prod(Obj(vn).h.dim));
            for qq = 1:numel(th)
                dd = resizeVol2(th(qq),Obj(vn).h);
                dat(qq,:) = dd(:);
            end
            
            ssData = [ssData; zscore(dat)];
        end
        
        ssData = zscore(ssData);
        
        Obj(vn).I(:) = 0;
        
        
        Obj(vn).Thresh = [-1 1];
        Obj(vn).clim = [-1 1];
        set(con(4),'String','-1 1');
        set(con(5),'String','-1 1');
        UpdateThreshold;
    end

    function ssUpdateConnMap(varargin)
        if ssConExp==0
            return
        end
        vn = get(con(21,1),'Value');
        if  vn ~= ssConLayer;
            return
        end
        
        MNI = Obj(1).point*Obj(1).h.mat';
        matLoc = round(MNI*inv(ssConHeader.mat)');
        
        thisLoc = sub2ind(ssConHeader.dim,matLoc(1),matLoc(2),matLoc(3));
        
        
        rad = get(findobj(paramenu3(2),'Checked','on'),'Label');
        rad = str2num(rad(1:end-2));
        [ml vi] = getMatCoord(ssConHeader,MNI(1:3),rad*2);
        
        seed = zscore(nanmean(ssData(:,vi),2));
        
        beta = pinv(seed)*ssData;        
%         keyboard;
        Obj(vn).I(:) = beta;
        UpdateThreshold;
    end

    function movieMode(varargin)
        nnn = spm_select(inf,'image');
        if size(nnn,1)==1
           ind = find(nnn==',');
           nnn = nnn(1:ind-1);
        end
        dh = spm_vol(nnn);
        
        for ii = 1:numel(dh)
            [t1 t2] = SliceAndDice3(dh(ii),MH,[],Obj(1).h,[0 NaN],[]);
            if ii == 1
                dd = nan([size(t1) numel(dh)]);
            end
            dd(:,:,:,ii) = t1;
        end
        rang = [min(dd(:)) min(dd(:))];
        
        xyz = [20 20 20]; jj = 1;
        tmp = [];
        tmp{1} = flipdim(rot90(squeeze(dd(xyz(1),:,:,jj)),1),1);
        tmp{2} = flipdim(flipdim(rot90(squeeze(dd(:,xyz(2),:,jj)),1),1),2);
        tmp{3} = flipdim(flipdim(rot90(squeeze(dd(:,:,xyz(3),jj)),1),1),2);
        
%         tmp = pcolor(Obj(ii).pos{2}, Obj(ii).pos{3}, Obj(ii).frame{opt});
%         if opt3; hand{opt}(ii) = tmp; end
%         colormap(gray(256)); shading interp; hold on;
        keyboard;
    end

    function adjustUnderlay(varargin)
        if varargin{1}==con(29)
            return
        end
        
        vn = get(con(21,1),'Value');
        
        if isempty(Obj(vn).Thresh)
           Obj(vn).Thresh = [min(Obj(vn).I(:)) max(Obj(vn).I(:))];
           Obj(vn).clim = [min(Obj(vn).I(:)) max(Obj(vn).I(:))];
        end
        
        
        if ~isfield(Obj(vn),'Range') || isempty(Obj(vn).Range)
           Obj(vn).Range = [min(Obj(vn).I(:)) max(Obj(vn).I(:)) max(Obj(vn).I(:))-min(Obj(vn).I(:)) ];
        end
        
        
        val = get(varargin{1},'Value');
        
        
        if varargin{1}==con(27)
            Obj(vn).clim(2) = (val*Obj(vn).Range(3))+Obj(vn).Range(1);
            %Obj(vn).clim(2) = val*Obj(vn).Thresh(2);
        end
        
        if varargin{1}==con(28)
            Obj(vn).clim(1) = (val*Obj(vn).Range(3))+Obj(vn).Range(1);
            %Obj(vn).clim(1) = val*Obj(vn).Thresh(2);
        end
        
        set(con(5),'String',[sprintf('%0.3f %0.3f',Obj(vn).clim(1),Obj(vn).clim(2))]);
        
        setupFrames(vn,1);
        updateGraphics([1 2 3],1);
        
         
        %%%%%
        yy = Obj(vn).clim(1):(Obj(vn).clim(2)-Obj(vn).clim(1))/255:Obj(vn).clim(2);
        axes(ax4); cla
        imagesc(yy,1,reshape(colmap(cmaps{Obj(vn).col},256),256,1,3));
        axis tight;
        set(ax4,'YDir','Normal','YAxisLocation','right','YTick',  unique([1 get(ax4,'YTick')]));
        set(ax4,'YTickLabel',round(min(yy):spm_range(yy)/(numel(get(ax4,'YTick'))-1):max(yy)));        
        set(ax4,'fontsize',6);

    end

    function UpdateCLims(varargin)
        vn = get(con(21,1),'Value');
        
        b = get(con(5),'String');
        
        if searchCellStr(',',{b})
            ind = find(b==',');
            b = [str2num(strtrim(b(1:ind-1))) str2num(strtrim(b(ind+1:end)))];
        else
            b = str2num(b);
        end
        
        if b(1)==-inf
            b(1) = Obj(vn).Range(1);
        end
        
         if b(2)==inf
            b(2) = Obj(vn).Range(2);
        end
        
        Obj(vn).clim = b;
        
        tmp1 = (Obj(vn).clim(2)-Obj(vn).Range(1))/Obj(vn).Range(3); if tmp1>1; tmp1=1; end; if tmp1<0; tmp1=0; end;
        set(con(27),'Value', tmp1);
        tmp1 = (Obj(vn).clim(1)-Obj(vn).Range(1))/Obj(vn).Range(3); if tmp1>1; tmp1=1; end; if tmp1<0; tmp1=0; end;
        set(con(28),'Value', tmp1);
        
        
        setupFrames(vn,1);
        updateGraphics([1 2 3],1);
        
        
        set(con(5),'String', sprintf('%0.3f %0.3f',b(1),b(2)));
        
        %%%%%%%%%%%%%
        yy = Obj(vn).clim(1):(Obj(vn).clim(2)-Obj(vn).clim(1))/255:Obj(vn).clim(2);
        axes(ax4); cla
        imagesc(yy,1,reshape(colmap(cmaps{Obj(vn).col},256),256,1,3));
        axis tight;
        set(ax4,'YDir','Normal','YAxisLocation','right','YTick',  unique([1 get(ax4,'YTick')]));
        set(ax4,'YTickLabel',round(min(yy):spm_range(yy)/(numel(get(ax4,'YTick'))-1):max(yy)));        
        set(ax4,'fontsize',6);

    end
end

function out = differencer(tmp,count)

if nargin == 1;
    count = numel(size(tmp));
end

out = diff(tmp,1,count);

if count ~= 2;
    out = differencer(out,count-1);
else
    if numel(size(out))>2   
        ss = size(out);
        out = reshape(out,size(out,1),prod(ss(2:end)));
    end
    return
end
end