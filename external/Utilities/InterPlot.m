function InterPlot(dat,Factors,FigNo)
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

for ii = 1:length(Factors);
    for jj = 1:length(Factors{ii});
        tmp = Factors{ii}{jj};
        ind = find(tmp=='_');
        tmp(ind)= ' ';
        Factors{ii}{jj} = tmp;
    end
end

for ii = 1:length(Factors);
    levs(ii) = numel(Factors{ii});
end

dat2 = reshape(dat,levs(2),levs(1));

if numel(Factors)~=2;
    error('This function is for Two-Factor interactions only.');
end

for ii = 1:size(dat2,1);
    for jj = 1:size(dat2,2);
        mm(ii,jj) = mean(dat2{ii,jj});
        %se(ii,jj) = 2*(std(dat2{ii,jj})./sqrt(numel(dat2{ii,jj})));% Standard Error
        se(ii,jj) = (std(dat2{ii,jj})./sqrt(numel(dat2{ii,jj})))*1.96;%95% CI
        er{ii,jj} = [mm(ii,jj)-se(ii,jj) mm(ii,jj)+se(ii,jj)];
%         er{ii,jj} = [mean(dat2{ii,jj})-(1*std(dat2{ii,jj})) mean(dat2{ii,jj})+(1*std(dat2{ii,jj}))];
    end
end

cols = {'b.-' 'r.-' 'k.-' 'm.-' 'c.-'};
cols2 = {'b-' 'r-' 'k-' 'm-' 'c-'};
if nargin == 3;
    figure(FigNo); clf
else
    figure; clf
end

adj = 0;
for ii = 1:size(mm,2);
    h(ii) = plot((1:numel(mm(:,ii)))+adj,mm(:,ii), cols{ii},'linewidth',2,'markersize',20);
    hold on;
    
    for jj = 1:size(mm,1)
       plot([jj+adj jj+adj],[er{jj,ii}(1) er{jj,ii}(2)],cols2{ii},'linewidth',2);
       plot([jj-.025+adj jj+.025+adj],[er{jj,ii}(1) er{jj,ii}(1)],cols2{ii},'linewidth',2);
       plot([jj-.025+adj jj+.025+adj],[er{jj,ii}(2) er{jj,ii}(2)],cols2{ii},'linewidth',2);
    end
    adj = adj+.05;
end

for ii = 1:length(h);
    uistack(h(ii),'bottom');
end

legend(h,Factors{1}, 'Fontsize',18,'location','Best');
if numel(Factors{2})<8
    set(gca,'XTick', 1:numel(Factors{2}),'XTickLabel',Factors{2}, 'Fontsize',16);
else
    set(gca,'XTick', 1:3:numel(Factors{2}),'XTickLabel',Factors{2}(1:3:end), 'Fontsize',16);
end

ax = axis;
ax(1:2) = [.8 numel(Factors{2})+.2];
axis(ax);

% return

% a = randn(30,1)+3;
% b = randn(30,1)-1;
% c = randn(30,1)-2;
% d = randn(30,1)+4;
% 
% Factors = {{'Normal' 'MCI'} {'PiB+' 'PiB-'}};
% 
% dat = {a b c d};


