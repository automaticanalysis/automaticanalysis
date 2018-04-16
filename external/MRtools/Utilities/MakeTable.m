function hands = MakeTable(P,savename)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
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

if iscell(P);
    tmp = P;
    clear P;
    P.Input = tmp;
end

P2 = struct('Input',[], ...
           'Caption',[], ...
           'Footer',[], ...
           'FontSize',10, ...
           'Font','Helvetica',...
           'FontWeight','normal',...
           'BackgroundColor',[1 1 1],...
           'RowColors',[1 1 1; 1 1 1],...
           'DataLine',2); 
       

if ~isfield(P,'Input') || isempty(P.Input)
    error('No input table was provided.'); 
else
    P2.Input = P.Input;
end

if isfield(P,'Caption') && ~isempty(P.Caption)
    P2.Caption = P.Caption;
end

if isfield(P,'FontSize') && ~isempty(P.FontSize)
    P2.FontSize = P.FontSize;
end

if isfield(P,'Font') && ~isempty(P.Font)
    P2.Font = P.Font;
end

if isfield(P,'BackgroundColor') && ~isempty(P.BackgroundColor)
    P2.BackgroundColor = P.BackgroundColor;
end

if isfield(P,'RowColors') && ~isempty(P.RowColors)
    P2.RowColors = P.RowColors;
end

if isfield(P,'DataLine') && ~isempty(P.DataLine)
    P2.DataLine = P.DataLine;
end

if isfield(P,'Footer') && ~isempty(P.Footer)
    P2.Footer = P.Footer;
end

if isfield(P,'FontWeight') && ~isempty(P.FontWeight)
    P2.FontWeight = P.FontWeight;
end
       
fsize = P2.FontSize;
ftype = P2.Font;
cols = P2.RowColors;
input = P2.Input;
Caption = P2.Caption;
Footer = P2.Footer;
FontWeight = P2.FontWeight;
% keyboard
tmp = input; tmp2 = input;
for ii = 1:size(tmp,1)
    for jj = 1:size(tmp,2)
        if searchCellStr('[<<<|>>>]',tmp(1,2))
            tmp{ii,jj} = regexprep(tmp{ii,jj},'<<<',' ');
            tmp{ii,jj} = regexprep(tmp{ii,jj},'>>>',' ');
            tmp2{ii,jj} = ' ';
        end
    end
end


%%
wid = [];
for ii = 1:size(tmp2,2);
    wid(ii) = size(char(tmp2(:,ii)),2);
end
wid = wid+4;

for ii = 1:size(input,1)
    for jj = 1:size(input,2)
        curJ = jj;
        if ~isempty(searchCellStr('>>>',input(ii,jj))) && ~isempty(regexprep(input{ii,jj},'>>>',''));
            str = regexprep(input{ii,jj},'>>>','');
            whcols = jj;
            while searchCellStr('>>>',input(ii,jj))
                jj = jj+1;
                if jj>size(input,2)
                    break
                end
                whcols = [whcols jj];
                
            end
            cl = sum(wid(whcols));
            if cl < (size(str,2)+(numel(whcols)*4));
                adj = ceil(((size(str,2)+(numel(whcols)*4))-cl)/numel(whcols));
                wid(whcols) = wid(whcols)+adj;
            end
        end
    end
end

n1 = size(input,2);
n2 = size(input,1);


figure; clf;
set(gcf,'Visible','off');

set(gcf,'color',P2.BackgroundColor); %'toolbar','none','menubar','none',
ss(1) = 20;
ss(2) = 400;
% ss(4) = (fsize/.5)*n2;
% ss(4) = (fsize/(8/fsize))*n2;
ss(3) = (sum(wid)*(fsize/(sqrt(fsize)/2)));
MsgHandle = uicontrol('style', 'text', 'position', [0 0 ss(3) 20],'visible','off','fontsize',fsize,'fontname',ftype,'fontweight',FontWeight);
[mess vpos] = textwrap(MsgHandle,{'AbX9g'});
ssize = get(0,'ScreenSize');
ss(4) = vpos(4)*n2;
ss(2) = ssize(4)-25-ss(4);

% keyboard;
if ~isempty(P2.Caption)
    MsgHandle = uicontrol('style', 'text', 'position', [0 0 ss(3) 20],'visible','off','fontsize',fsize-0,'fontname',ftype,'fontweight',FontWeight);
    [NewCap CapPos]=textwrap(MsgHandle,Caption);
    delete(MsgHandle);
    
    ss(4) = ss(4)+CapPos(4)+10;
    capadj = ((CapPos(4)*.8)+10)/ss(4);
else
    CapPos = [0 0 0 0];
    capadj = 0;
end


if ~isempty(P2.Footer)
    MsgHandle = uicontrol('style', 'text', 'position', [0 0 ss(3) 20],'visible','off','fontsize',fsize-0,'fontname',ftype,'fontweight',FontWeight);
    [NewFoot FootPos]=textwrap(MsgHandle,P2.Footer);
    delete(MsgHandle);
    
    ss(4) = ss(4)+FootPos(4)+10;
    footadj = ((FootPos(4)*.8)+10)/ss(4);
else
    FootPos = [0 0 0 0];
    footadj = 0;
end

set(gcf,'Position',ss);
% keyboard;
axes; %axis([0 0 1 1]);
set(gca,'units','normalized','position',[0 0 1 1],'visible','off');



hst = [0 cumsum(wid)]./sum(wid);

i1 = 1-(5/ss(4))-capadj;
i2 = (5/ss(4))+footadj;
inc = (i1-i2)/n2;
vst = (i1-inc):-inc:i2;
%%
wid = [0 wid./sum(wid) 1];
a = cumsum(wid);
wid = [a(1:end-2)' diff(a(1:end-1))'];

%%
ll = [];
hhh = []; th = [];
hh = []; ww = 0; str = []; adj = 0; go = 0; tt = 0;
for jj = 1:size(input,2)
    for ii = 1:size(input,1)
        if isempty(setdiff(input{ii,jj}, ' '))
            continue
        end
        if isempty(setdiff(input{ii,jj}, '<'))
            continue
        end
        if isempty(setdiff(input{ii,jj}, '>'))
            continue
        end
        
        if ~isempty(regexp(input{ii,jj},'>>>$'))
             h1 = hst(jj);
             v1 = vst(ii);
             tt = inc;
            for kk = jj:size(input,2);
                ww = ww+(wid(kk)./sum(wid));
                str = [str input{ii,kk}(1:end-3)];
                go = 1;
                if isempty(regexp(input{ii,kk},'>>>$'))
                    break
                end
            end
        elseif ~isempty(regexp(input{ii,jj},'<<<$'))
            h1 = hst(jj);
            v1 = vst(ii);
            ww = (wid(jj)./sum(wid));
            for kk = ii:size(input,1);                
                str = [str input{kk,jj}(1:end-3)];
                tt = tt+inc;
                v1 = vst(kk);
                go = 2;
                if isempty(regexp(input{kk,jj},'<<<$'))
                    break
                end
            end
        else
             ww = ww+(wid(jj)./sum(wid));
             str = [str input{ii,jj}];
              h1 = hst(jj);
              v1 = vst(ii);
              tt = inc;
        end
        
                
%         hh(end+1,1) = annotation(gcf,'textbox',[h1 v1 ww tt],'units','normalized','string',str, ...
%             'fontsize',fsize,'fontweight',FontWeight,'fontname',ftype,...
%             'HorizontalAlignment', 'Center','VerticalAlignment','Middle', ...
%             'BackgroundColor', cols(mod(jj+1,2)+1,:),'LineStyle','none');
        
%         hh(end+1,1) = uicontrol(gcf,'style','text','units','normalized','Position',[h1 v1 ww tt],'string',str, ...
%             'fontsize',fsize,'fontweight',FontWeight,'fontname',ftype,...
%             'HorizontalAlignment', 'Center', ... %'VerticalAlignment','Middle',
%             'BackgroundColor', cols(mod(jj+1,2)+1,:));%,'LineStyle','none'

%         if jj==1
%              tmp = text((h1+h1+ww)/2, (v1+v1+tt)/2,str, ...
%             'fontsize',fsize,'fontweight',FontWeight,'fontname',ftype,...
%             'HorizontalAlignment', 'Center','VerticalAlignment','Middle', ...
%             'BackgroundColor', 'none');%,'LineStyle','none
%         
%             b = get(tmp,'extent');
%             th(ii,1) = annotation(gcf,'rectangle',[0 b(2) 1 b(4)],'FaceColor',cols(mod(ii+1,2)+1,:),'FaceAlpha',1,'linestyle','none');
%             delete(tmp)
%         end
        

        hh(end+1,1) = annotation(gcf,'textbox',[wid(jj,1) v1*1 wid(jj,2) tt*1],'string',str, ...
            'fontsize',fsize,'fontweight',FontWeight,'fontname',ftype,...
            'HorizontalAlignment', 'Center','VerticalAlignment','Middle', ...
            'BackgroundColor', cols(mod(ii+1,2)+1,:),'linestyle','none');

%         hh(end+1,1) = text((h1+h1+ww)/2, (v1+v1+tt)/2,str, ...
%             'fontsize',fsize,'fontweight',FontWeight,'fontname',ftype,...
%             'HorizontalAlignment', 'Center','VerticalAlignment','Middle', ...
%             'BackgroundColor', 'none');%,'LineStyle','none'
%         % 'BackgroundColor', cols(mod(ii+1,2)+1,:)

        if go==1
            x = diff([h1 h1+ww])*.05;
            ll(end+1,1) = annotation(gcf,'line',[h1+x h1+ww-x],[v1 v1],'units','normalized','color','k',...
                'linewidth',1.5,'LineStyle','-');
        end
        
        if go==2
            ll(end+1,1) = annotation(gcf,'line',[h1+(ww/10) .99],[v1 v1],'units','normalized','color','k',...
                'linewidth',1.5,'LineStyle','-');
        end
        hhh(ii,jj) = hh(end);
        ww = 0;
        str = [];
        adj = 0;
        go = 0;
        tt = 0;
    end
%     keyboard;
%     %%
%     ex = [];
%     for kk = 1:numel(hh);
%         ex(kk,:) = get(hh(kk),'Position');
%     end
%     ex = max(ex);
%     set(hh,'Position',ex);
    %ex = [min(ex(:,1)) max(ex(:,2)) min(ex(:,3)) max(ex(:,4))];
%     set(hh,'Extent',ex);
end

%%

ll(end+1,1) = annotation(gcf,'line',[0 1],[vst(1)+inc vst(1)+inc],'units','normalized','color','k',...
                'linewidth',2,'LineStyle','-');
ll(end+1,1) = annotation(gcf,'line',[0 1],[vst(end) vst(end)],'units','normalized','color','k',...
                'linewidth',2,'LineStyle','-'); 

if ~isempty(P2.DataLine) && ~isnan(P2.DataLine(1));
    %keyboard;
    tadj = 1/ss(4);
    if P2.DataLine(1)-1 > 0
        ll(end+1,1) = annotation(gcf,'line',[0 1],[vst(P2.DataLine(1)-1)-tadj vst(P2.DataLine(1)-1)-tadj],'units','normalized','color','k',...
            'linewidth',1,'LineStyle','-');
        ll(end+1,1) = annotation(gcf,'line',[0 1],[vst(P2.DataLine(1)-1)+tadj vst(P2.DataLine(1)-1)+tadj],'units','normalized','color','k',...
            'linewidth',1,'LineStyle','-');
    end
end
            

if ~isempty(P2.Caption)
%     for ii = 2:numel(NewCap);
%         NewCap{ii} = ['  ' NewCap{ii}];
%     end
    cap = annotation('textbox',[0 vst(1)+(inc) 1 capadj], 'units','normalized','string',NewCap, ...
        'fontsize',fsize-0,'fontweight',FontWeight,'fontname',ftype,...
        'HorizontalAlignment', 'Left','VerticalAlignment','Middle', ...
        'BackgroundColor', 'w','LineStyle','none'); %,'FitBoxToText','on'
else
    cap = [];
end

if ~isempty(Footer)
%     for ii = 2:numel(Footer);
%         Footer{ii} = ['  ' Footer{ii}];
%     end
    
    foot = annotation('textbox',[0 (5/ss(4))  1 footadj], 'units','normalized','string',Footer, ...
        'fontsize',fsize-0,'fontweight',FontWeight,'fontname',ftype,...
        'HorizontalAlignment', 'Left','VerticalAlignment','Middle', ...
        'BackgroundColor', 'w','LineStyle','none'); %,'FitBoxToText','on'
else
    foot = [];
end


if ~isempty(P2.DataLine);
    for ii = 2:numel(P2.DataLine);
        if isnan(P2.DataLine)
            continue
        end
        tadj = 1/ss(4);
        ll(end+1,1) = annotation(gcf,'line',[0 1],[vst(P2.DataLine(ii)-1)-tadj vst(P2.DataLine(ii)-1)-tadj],'units','normalized','color','k',...
            'linewidth',1,'LineStyle','-');
        ll(end+1,1) = annotation(gcf,'line',[0 1],[vst(P2.DataLine(ii)-1)+tadj vst(P2.DataLine(ii)-1)+tadj],'units','normalized','color','k',...
            'linewidth',1,'LineStyle','-');
    end
end

uistack(ll,'top');

hands.hh = hh;
hands.ll = ll;
hands.cap = cap;
hands.foot = foot;

% set(hands.hh,'FontUnits','Normalized')
set(gcf,'Visible','on');

if nargin==2
    ss = get(gcf,'Position');
    set(gcf,'PaperUnits','points');

    adj = 5;
    set(gcf,'PaperSize',[ss(3) ss(4)+adj],'PaperPosition',[0 0 ss(3) ss(4)+adj]);
        
    %saveas(gcf,[savename '.fig'],'fig');
    %saveas(gcf,[savename '.pdf'],'pdf');    
    %system(['evince ' savename '.pdf &'])
    
    saveas(gcf,[savename '.fig'],'fig');
    saveas(gcf,[savename '.eps'],'eps');    
    system(['epstopdf ' savename '.eps']);
    delete([savename '.eps'])
%     system(['evince ' savename '.pdf &']) 
end

