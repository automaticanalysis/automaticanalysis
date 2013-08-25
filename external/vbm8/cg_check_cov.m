function cg_check_cov(vargin)
%cg_check_cov	to check covriance across sample
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. normalized images)
% An output image will be save with SD at each voxel.
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_check_cov.m 425 2011-08-22 14:40:10Z gaser $

global fname jY h1 h2 YpY slice_array
rev = '$Rev: 425 $';

if nargin == 1
  P = char(vargin.data);
  norm = vargin.scale;
  if isempty(vargin.nuisance)
    nuisance = [];
  else
    nuisance = vargin.nuisance.c;
  end
  slice_mm = vargin.slice;
  gap = vargin.gap;
end

if nargin < 1
  P = spm_select(Inf,'image','Select images');
end

V = spm_vol(deblank(P));
n = size(P,1);

if length(V)>1 & any(any(diff(cat(1,V.dim),1,1),1))
  error('images don''t all have same dimensions')
end
if max(max(max(abs(diff(cat(3,V.mat),1,3))))) > 1e-8
  error('images don''t all have same orientation & voxel size')
end

if nargin < 1
  norm = spm_input('Prop. scaling (e.g. for T1- or modulated images)?',1,'yes|no',[1 0],2);
  def_nuis = spm_input('Variable to covariate out (nuisance parameter)?','+1','yes|no',[1 0],2);
  if def_nuis
    nuisance = spm_input('Nuisance parameter:','+1','r',[],n);
  else
    nuisance = [];
  end
	slice_mm = spm_input('Slice [mm]?','+1','e',0,1);
	gap = spm_input('Gap for slices to speed up','+1','e',0,1);
end

if ~isempty(nuisance)
  if size(nuisance,2) ~= n
    nuisance = nuisance';
  end
end

% voxelsize and origin
vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
Orig = V(1).mat\[0 0 0 1]';

% range
range = ([1 V(1).dim(3)] - Orig(3))*vx(3);

% calculate slice from mm to voxel
sl = round(slice_mm/vx(3)+Orig(3));
while (sl < 1) | (sl > V(1).dim(3))
	slice_mm = spm_input(['Slice (in mm) [' num2str(range(1)) '...' num2str(range(2)) ']'],1,'e',0);
	sl = round(slice_mm/vx(3)+Orig(3));
end

% global scaling
if norm
  gm=zeros(size(V,1),1);
  disp('Calculating globals...');
  for i=1:size(V,1), gm(i) = spm_global(V(i)); end
  gm_all = mean(gm);
  for i=1:n
    V(i).pinfo(1:2,:) = gm_all*V(i).pinfo(1:2,:)/gm(i);
  end
end

%-Start progress plot
%-----------------------------------------------------------------------
vol = zeros(prod(V(1).dim(1:2)), n);
YpY = zeros(n);

spm_progress_bar('Init',V(1).dim(3),'Check covariance','planes completed')
slice_array = zeros([V(1).dim(1:2) n]);

% consider gap for slices to speed up the process
slices = 1:gap:V(1).dim(3);
[mn, ind] = min(abs(slices-sl));
slices = slices - slices(ind) + sl;
if slices(1) < 1
  slices = slices(2:end);
end
if slices(end) > V(1).dim(3)
  slices = slices(1:end-1);
end
if slices(end) < V(1).dim(3) - gap
  slices = [slices slices(end)+gap]; 
end

for j=slices

  M  = spm_matrix([0 0 j]);

  for i = 1:n
    img = spm_slice_vol(V(i),M,V(1).dim(1:2),[1 0]);
    vol(:,i) = img(:);
  end

	% get slice data
  if j == sl
  	slice_array = reshape(vol,[V(1).dim(1:2) n]);
  end

  mean_slice = mean(reshape(vol,[V(i).dim(1:2) n]),3);
  mask = find(mean_slice ~= 0 & ~isnan(mean_slice));
  % remove nuisance and calculate again mean
  if ~isempty(nuisance) 
    vol(mask,:) = vol(mask,:) - vol(mask,:)*pinv(nuisance)*nuisance;
  end

  if ~isempty(mask)
    % make sure data is zero mean
    tmp_vol = vol(mask,:);
    tmp_vol = tmp_vol - repmat(mean(tmp_vol,1), [length(mask) 1]);
	  YpY = YpY + (tmp_vol'*tmp_vol)/n;
  end 
  spm_progress_bar('Set',j);  
end

spm_progress_bar('Clear');

% normalize YpY
d = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd = d*d';
YpY = YpY./(dd+eps);
t = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n+1:end) = sign(diag(YpY));

YpYsum = sum(YpY,1);
[iY, jY] = sort(YpYsum, 2, 'descend');
YpYsorted = YpY(jY,jY);
Nsorted = P(jY,:);

% extract mean covariance
mean_cov = zeros(n,1);
for i=1:n
	% extract row for each subject
	cov0 = YpY(i,:);
	% remove cov with its own
	cov0(i) = [];
	mean_cov(i) = mean(cov0);
end

threshold_cov = mean(mean_cov) - 2*std(mean_cov);

[tmp fname] = spm_str_manip(char(V.fname),'C');
fprintf('Compressed filenames: %s  \n',tmp);

% print suspecious files with cov>0.9
YpY_tmp = YpY - tril(YpY);
[indx, indy] = find(YpY_tmp>0.9);
if ~isempty(indx) & (sqrt(length(indx)) < 0.5*n)
  fprintf('\nUnusual large covariances (check that subjects are not identical):\n');
  for i=1:length(indx)
    % exclude diagonal
    if indx(i) ~= indy(i)
      % report file with lower mean covariance first
      if mean_cov(indx(i)) < mean_cov(indy(i))
        fprintf('%s and %s: %3.3f\n',fname.m{indx(i)},fname.m{indy(i)},YpY(indx(i),indy(i)));
      else
        fprintf('%s and %s: %3.3f\n',fname.m{indy(i)},fname.m{indx(i)},YpY(indy(i),indx(i)));
      end
    end
  end
end

% sort files
fprintf('\nMean covariance for data below 2 standard deviations:\n');
[mean_cov_sorted, ind] = sort(mean_cov,'descend');
n_thresholded = min(find(mean_cov_sorted < threshold_cov));

for i=n_thresholded:n
  fprintf('%s: %3.3f\n',V(ind(i)).fname,mean_cov_sorted(i));
end

Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
FS    = spm('FontSizes');

xpos = 2*(0:n-1)/(n-1);
for i=1:n
  text(xpos(i),mean_cov(i),fname.m{i},'FontSize',FS(8),'HorizontalAlignment','center')
end

hold on
cg_boxplot({mean_cov});
set(gca,'XTick',[],'XLim',[-.25 2.25]);
if max(mean_cov) > min(mean_cov)
  set(gca,'YLim',[0.9*min(mean_cov) 1.1*max(mean_cov)]);
end
title(sprintf('Boxplot: mean covariance  \nCommon filename: %s*%s',fname.s,fname.e),'FontSize',FS(12),'FontWeight','Bold');
ylabel('<----- low (poor quality) --- mean covariance --- large (good quality)------>  ','FontSize',FS(10),'FontWeight','Bold');
xlabel('<----- first --- file order --- last ------>  ','FontSize',FS(10),'FontWeight','Bold');
hold off

% covariance
f = figure(4);
ws = spm('Winsize','Graphics');

set(f,'Name','Click in image to get file names','NumberTitle','off');
h = datacursormode(f);
set(h,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
set(f,'MenuBar','none','Position',[10 10 ws(3) ws(3)]);

cmap = [gray(64); hot(64)];

% scale YpY to 0..1
mn = min(YpY(:));
mx = max(YpY(:));
YpY_scaled = (YpY - mn)/(mx - mn);
YpYsorted_scaled = (YpYsorted - mn)/(mx - mn);

% show upper right triangle in gray
ind_tril = find(tril(ones(size(YpY))));
ima = YpY_scaled;
ima(ind_tril) = 0;
ima(ind_tril) = 1 + 1/64 + YpY_scaled(ind_tril);
image(64*ima)
a = gca;
set(a,'XTickLabel','','YTickLabel','');
axis image
xlabel('<----- first --- file order --- last ------>  ','FontSize',10,'FontWeight','Bold');
ylabel('<----- last --- file order --- first ------>  ','FontSize',10,'FontWeight','Bold');
title('Covariance','FontSize',12,'FontWeight','Bold');
colormap(cmap)

% ordered covariance
f = figure(5);
set(f,'Name','Click in image to get file names','NumberTitle','off');
h = datacursormode(f);
set(h,'UpdateFcn',@myupdatefcn_ordered,'SnapToDataVertex','on','Enable','on');
set(f,'MenuBar','none','Position',[11+ws(3) 10 ws(3) ws(3)]);

% show upper right triangle in gray
ind_tril = find(tril(ones(size(YpY))));
ima = YpYsorted_scaled;
ima(ind_tril) = 0;
ima(ind_tril) = 1 + 1/64 + YpYsorted_scaled(ind_tril);
image(64*ima)
if n_thresholded <= n
	hold on
	line([n_thresholded-0.5, n_thresholded-0.5], [0.5,n_thresholded-0.5])
	line([0.5,n_thresholded-0.5],[n_thresholded-0.5, n_thresholded-0.5])
	hold off
end
a = gca;
set(a,'XTickLabel','','YTickLabel','');
axis image
xlabel('<----- high --- mean covariance --- low ------>  ','FontSize',10,'FontWeight','Bold');
ylabel('<----- low --- mean covariance --- high ------>  ','FontSize',10,'FontWeight','Bold');
title({'Sorted Covariance','Blue line indicates 2-SD threshold'},'FontSize',12,'FontWeight','Bold');
colormap(cmap)

% slice preview
f = figure(6);
set(f,'MenuBar','none','Name','Slice preview','NumberTitle','off','Position',[12+2*ws(3) 10 2*V(1).dim(2) 4*V(1).dim(1)]);

% Close button
hCloseButton = uicontrol(f,...
    		'position',[V(1).dim(2)-40 4*V(1).dim(1)-25 80 20],...
    		'style','Pushbutton',...
    		'string','Close',...
    		'callback','try close(6); end; try close(5); end;try close(4);end;',...
    		'ToolTipString','Close windows',...
        'Interruptible','on','Enable','on');

% range 0..64
slice_array = 64*slice_array/max(slice_array(:));

% check for replicates
for i=1:n
  for j=1:n
  if (i>j) & (mean_cov(i) == mean_cov(j))
    [s,differ] = unix(['diff ' V(i).fname ' ' V(j).fname]);
    if (s==0), fprintf(['\nWarning: ' V(i).fname ' and ' V(j).fname ' are same files?\n']); end
  end
  end
end

%-End
%-----------------------------------------------------------------------
spm_progress_bar('Clear')

show = spm_input('Show files with poorest cov?','+1','yes|no',[1 0],2);
if show
  number = min([n 15]);
  number = spm_input('How many files ?','+1','e',number);
  
  list = str2mat(V(ind(n:-1:1)).fname);
  list2 = list(1:number,:);
  spm_check_registration(list2)
end
return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global fname jY h1 h2 YpY slice_array
pos = get(event_obj, 'Position');
h = gca;

x = pos(1);
y = pos(2);

txt = {sprintf('Covariance: %3.3f',YpY(x,y)),fname.m{x},fname.m{y}};

f = figure(6);
img = [slice_array(:,:,x); slice_array(:,:,y)];
image(img)
p = get(f,'Position');
p(3:4) = 2*size(img');
set(f,'Position',p);
set(f,'MenuBar','none','Colormap',gray);
set(gca,'XTickLabel','','YTickLabel','');
h = xlabel({['Top: ',fname.m{x}],['Bottom: ',fname.m{y}]});
set(h,'Interpreter','none');
axis image

return

%-----------------------------------------------------------------------
function txt = myupdatefcn_ordered(obj, event_obj)
%-----------------------------------------------------------------------
global fname jY h1 h2 YpY slice_array
pos = get(event_obj, 'Position');
h = gca;

x = jY(pos(1));
y = jY(pos(2));

txt = {sprintf('Covariance: %3.3f',YpY(x,y)),fname.m{x},fname.m{y}};

f = figure(6);
img = [slice_array(:,:,x); slice_array(:,:,y)];
image(img)
p = get(f,'Position');
p(3:4) = 2*size(img');
set(f,'Position',p);
set(f,'MenuBar','none','Colormap',gray);
set(gca,'XTickLabel','','YTickLabel','');
h = xlabel({['Top: ',fname.m{x}],['Bottom: ',fname.m{y}]});
set(h,'Interpreter','none');
axis image

return

%-----------------------------------------------------------------------
function s = cg_boxplot (data,notched,symbol,vertical,maxwhisker)
%-----------------------------------------------------------------------
% usage: s = cg_boxplot (data,notched,symbol,vertical,maxwhisker);
%
% The box plot is a graphical display that simultaneously describes several 
% important features of a data set, such as center, spread, departure from 
% symmetry, and identification of observations that lie unusually far from
% the bulk of the data.
%
% data is a matrix with one column for each dataset, or data is a cell
% vector with one cell for each dataset.
% notched = 1 produces a notched-box plot. Notches represent a robust 
% estimate of the uncertainty about the median.
% notched = 0 (default) produces a rectangular box plot. 
% notched in (0,1) produces a notch of the specified depth.
% notched values outside [0,1] are amusing if not exactly practical.
% symbol sets the symbol for the outlier values, default symbol for
% points that lie outside 3 times the interquartile range is 'o',
% default symbol for points between 1.5 and 3 times the interquartile
% range is '+'. 
% Examples
% symbol = '.' points between 1.5 and 3 times the IQR is marked with
% '.' and points outside 3 times IQR with 'o'.
% symbol = ['x','*'] points between 1.5 and 3 times the IQR is marked with
% 'x' and points outside 3 times IQR with '*'.
% vertical = 0 makes the boxes horizontal, by default vertical = 1.
% maxwhisker defines the length of the whiskers as a function of the IQR
% (default = 1.5). If maxwhisker = 0 then boxplot displays all data  
% values outside the box using the plotting symbol for points that lie
% outside 3 times the IQR.   
%
% The returned matrix s has one column for each dataset as follows:
%
%    1  minimum
%    2  1st quartile
%    3  2nd quartile (median)
%    4  3rd quartile
%    5  maximum
%    6  lower confidence limit for median
%    7  upper confidence limit for median
%
% Example
%
%   title("Grade 3 heights");
%   tics("x",1:2,["girls";"boys"]);
%   axis([0,3]);
%   boxplot({randn(10,1)*5+140, randn(13,1)*8+135});
%

% Author: Alberto Terruzzi <t-albert@libero.it>
% Version: 1.4
% Created: 6 January 2002
% Copyright (C) 2002 Alberto Terruzzi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% modified by Christian Gaser (christian.gaser@uni-jena.de)
% original version was written for octave by Alberto Terruzzi

% assign parameter defaults
if nargin < 1 || nargin > 5
   error('s = boxplot (data,notch,symbol,vertical,maxwhisker)')
end
if nargin < 5, maxwhisker = 1.5; end
if nargin < 4, vertical = 1; end
if nargin < 3, symbol = ['+','o']; end
if nargin < 2, notched = 0; end

if length(symbol)==1, symbol(2)=symbol(1); end

if notched==1, notched=0.25; end
a = 1-notched;

% figure out how many data sets we have
if iscell(data), 
  nc = length(data);
else
  if isvector(data), data = data(:); end
  nc = columns(data);
end

% compute statistics
% s will contain
%    1,5    min and max
%    2,3,4  1st, 2nd and 3rd quartile
%    6,7    lower and upper confidence intervals for median
s = zeros(7,nc);
box = zeros(1,nc);
whisker_x = ones(2,1)*[1:nc,1:nc];
whisker_y = zeros(2,2*nc);
outliers_x = [];
outliers_y = [];
outliers2_x = [];
outliers2_y = [];

for i=1:nc
  % Get the next data set from the array or cell array
  if iscell(data)
    col = data{i}(:);
  else
    col = data(:,i);
  end
  % Skip missing data
  col(isnan(col)) = [];
  % Remember the data length
  nd = length(col);
  box(i) = nd;
  if (nd > 1)
    % min,max and quartiles
    s(1:5,i) = [min(col) prctile(col,[25 50 75]) max(col)]';
    % confidence interval for the median
    est = 1.57*(s(4,i)-s(2,i))/sqrt(nd);
    s(6,i) = max([s(3,i)-est, s(2,i)]);
    s(7,i) = min([s(3,i)+est, s(4,i)]);
    % whiskers out to the last point within the desired inter-quartile range
    IQR = maxwhisker*(s(4,i)-s(2,i));
    whisker_y(:,i) = [min(col(col >= s(2,i)-IQR)); s(2,i)];
    whisker_y(:,nc+i) = [max(col(col <= s(4,i)+IQR)); s(4,i)];
    % outliers beyond 1 and 2 inter-quartile ranges
    outliers = col((col < s(2,i)-IQR & col >= s(2,i)-2*IQR) | (col > s(4,i)+IQR & col <= s(4,i)+2*IQR));
    outliers2 = col(col < s(2,i)-2*IQR | col > s(4,i)+2*IQR);
    outliers_x = [outliers_x; i*ones(size(outliers))];
    outliers_y = [outliers_y; outliers];
    outliers2_x = [outliers2_x; i*ones(size(outliers2))];
    outliers2_y = [outliers2_y; outliers2];
  elseif (nd == 1)
    % all statistics collapse to the value of the point
    s(:,i) = col;
    % single point data sets are plotted as outliers.
    outliers_x = [outliers_x; i];
    outliers_y = [outliers_y; col];
  else
    % no statistics if no points
    s(:,i) = NaN;
  end
end

% Note which boxes don't have enough stats
chop = find(box <= 1);
    
% Draw a box around the quartiles, with width proportional to the number of
% items in the box. Draw notches if desired.
box = box*0.3/max(box);
quartile_x = ones(11,1)*[1:nc] + [-a;-1;-1;1;1;a;1;1;-1;-1;-a]*box;
quartile_y = s([3,7,4,4,7,3,6,2,2,6,3],:);

% Draw a line through the median
median_x = ones(2,1)*[1:nc] + [-a;+a]*box;
median_y = s([3,3],:);

% Chop all boxes which don't have enough stats
quartile_x(:,chop) = [];
quartile_y(:,chop) = [];
whisker_x(:,[chop,chop+nc]) = [];
whisker_y(:,[chop,chop+nc]) = [];
median_x(:,chop) = [];
median_y(:,chop) = [];

% Add caps to the remaining whiskers
cap_x = whisker_x;
cap_x(1,:) = cap_x(1,:) - 0.05;
cap_x(2,:) = cap_x(2,:) + 0.05;
cap_y = whisker_y([1,1],:);

% Do the plot
if vertical
	plot(quartile_x, quartile_y, 'b-')
	hold on
	plot(whisker_x, whisker_y, 'b-')
	plot(cap_x, cap_y, 'b-')
	plot(median_x, median_y, 'r-')
	plot(outliers_x, outliers_y, [symbol(1),'r'])
        plot(outliers2_x, outliers2_y, [symbol(2),'r']);
else
	plot(quartile_y, quartile_x, 'b-')
	hold on
	plot(whisker_y, whisker_x, 'b-')
	plot(cap_y, cap_x, 'b-')
	plot(median_y, median_x, 'r-')
	plot(outliers_y, outliers_x, [symbol(1),'r'])
        plot(outliers2_y, outliers2_x, [symbol(2),'r']);
end

hold off
return

%-----------------------------------------------------------------------
function y = prctile(x,p);
%-----------------------------------------------------------------------
%PRCTILE gives the percentiles of the sample in X.
%   Y = PRCTILE(X,P) returns a value that is greater than P percent
%   of the values in X. For example, if P = 50  Y is the median of X. 
%
%   P may be either a scalar or a vector. For scalar P, Y is a row   
%   vector containing Pth percentile of each column of X. For vector P,
%   the ith row of Y is the P(i) percentile of each column of X.

%   Copyright 1993-2002 The MathWorks, Inc. 

[prows pcols] = size(p);
if prows ~= 1 & pcols ~= 1
    error('P must be a scalar or a vector.');
end
if any(p > 100) | any(p < 0)
    error('P must take values between 0 and 100');
end

if (~any(isnan(x)))
   y = prctilecol(x,p);
else                    % if there are NaNs, process each column
   if (size(x,1) == 1)
      x = x';
   end
   c = size(x,2);
   np = length(p);
   y = zeros(np,c);
   for j=1:c
      xx = x(:,j);
      xx = xx(~isnan(xx));
      y(:,j) = prctilecol(xx,p)';
   end
   if (min(size(x)) == 1)
      y = y';
   end
end

return
      
%-----------------------------------------------------------------------
function y = prctilecol(x,p);
%-----------------------------------------------------------------------
xx = sort(x);
[m,n] = size(x);

if m==1 | n==1
    m = max(m,n);
	if m == 1,
	   y = x*ones(length(p),1);
	   return;
	end
    n = 1;
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx(:); max(x)];
else
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx; max(x)];
end

q = [0 q 100];
y = interp1(q,xx,p);

return
