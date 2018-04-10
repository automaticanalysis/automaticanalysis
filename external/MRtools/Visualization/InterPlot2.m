function InterPlot2(dat,Factors,altX)
%%% Displays in an interaction line plot with 1*SE error bars.  This
%%% function only works for two-way interactions.
%%%
%%% Inputs:
%%% dat = a cell array where each cell contains the data for a particular
%%% confition.
%%%
%%% Factors =  a double cell array of condition names (strings).
%%%
%%% FogNo = the figure number to display the plot.
%%%
%%% Example:  InterPlot({randn(30,1)+1,randn(30,1)+2,randn(30,1)+3,randn(30,1)+2,randn(30,1)+4,randn(30,1)+6}, {{'PiB+' 'PiB-'} {'Young' 'Normal' 'MCI'}});
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
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
FigNo = gcf;

if size(Factors,2)~=2;
    error('This function is for Two-Way interactions only.');
end

j1 = unique(Factors(:,1),'stable');
j2 = unique(Factors(:,2),'stable');
% j1 = unique(Factors(:,1));
% j2 = unique(Factors(:,2));

%double check the ordering for altX
if nargin<3
    altX = (1:numel(j1))';
    altX = repmat(altX,1,numel(j2));
else
    if any(size(altX)==1)
        altX = altX(:);
        altX = repmat(altX,1,numel(j2));
    end
end
    
for ii = 1:numel(j1);
    for jj = 1:numel(j2);
        if isobject(Factors)
            warning('not yet configured');
            return
        end
        if iscell(Factors)
            i1 = strmatch(j1{ii},Factors(:,1),'exact');
            i2 = strmatch(j2{jj},Factors(:,2),'exact');
            
            i3 = intersect(i1,i2);
            
            mm(ii,jj) = mean(dat(i3));
            %se(ii,jj) = 2*(std(dat(i3))./sqrt(numel(dat(i3))));% Standard Error
            se(ii,jj) = (std(dat(i3))./sqrt(numel(dat(i3))))*1.96;%95% CI
            er{ii,jj} = [mm(ii,jj)-se(ii,jj) mm(ii,jj)+se(ii,jj)];
            %er{ii,jj} = [mean(dat(i3))-(1*std(dat(i3))) mean(dat(i3))+(1*std(dat(i3)))];
        end
    end
end

if size(mm,2) == 1
    cols = [0 0 1];
else
    cols = zeros(size(mm,2),3);
    cols(:,3) = 1:-1/(size(mm,2)-1):0;
    cols(:,1) = 0:1/(size(mm,2)-1):1;
end

figure(FigNo); clf


adj = 0;
for ii = 1:size(mm,2);
    h(ii) = plot(altX(:,ii)+adj,mm(:,ii), '.-','color',cols(ii,:),'linewidth',2,'markersize',20); hold on;

    for jj = 1:size(mm,1)
        loc = altX(jj,ii);
        
        plot([loc+adj loc+adj],[er{jj,ii}(1) er{jj,ii}(2)],'-','color',cols(ii,:),'linewidth',2);
        plot([loc-.025+adj loc+.025+adj],[er{jj,ii}(1) er{jj,ii}(1)],'-','color',cols(ii,:),'linewidth',2);
        plot([loc-.025+adj loc+.025+adj],[er{jj,ii}(2) er{jj,ii}(2)],'-','color',cols(ii,:),'linewidth',2);
    end
    adj = adj+.05;
end

for ii = 1:length(h);
    uistack(h(ii),'bottom');
end

% keyboard;
legend(h,j2, 'Fontsize',18,'location','Best');
if numel(j1)<15
    set(gca,'XTick', mean(altX,2),'XTickLabel',j1, 'Fontsize',16);
else
    set(gca,'XTick',mean(altX(1:2:end,:),2),'XTickLabel',j1(1:2:end), 'Fontsize',16);
end

ax = axis;
ax(1:2) = [altX(1)-.2 altX(end)+.2];

axis(ax);



