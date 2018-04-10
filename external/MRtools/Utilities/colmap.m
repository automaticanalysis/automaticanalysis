function [cm list] = colmap(map,depth,custom)
%%% Written by Aaron P. Schultz - aschultz@martinos.org
%%%
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

list = {'gray' 'hot' 'cool' 'jet' 'hsv' 'pink' 'spring' 'summer' 'autumn' 'winter' 'bone' 'copper' 'red' 'red2' 'green' 'blue' 'yellow' 'yellow2' 'purple' 'purple2' 'cyan2' 'mag2' 'ocean' 'orange' 'my' 'ym' 'mby' 'cy' 'yc' 'mc' 'cm' ...
    'jet1' 'jet2' 'cool2' 'rb1' 'rb2' 'hc1' 'hc2' 'sunrise' ...
    'cbd_rd_yl_gn' 'cbd_rd_yl_bu' 'cbd_rd_gy' 'cbd_rd_bu' 'cbd_pu_or' 'cbd_pu_gn' 'cbd_pi_yl_gn' 'cbd_br_bg' ...
    'cbs_yl_or_rd' 'cbs_yl_or_br' 'cbs_yl_gn_bu' 'cbs_yl_gn' 'cbs_rd_pu' 'cbs_pu' 'cbs_pu_rd' 'cbs_pu_bu_gn' 'cbs_pu_bu' 'cbs_or_rd' 'cbs_or' 'cbs_gr' 'cbs_gn' 'cbs_gn_bu' 'cbs_bu_pu' 'cbs_bn_gn' 'cbs_bu' ...
    'cbq_set1' 'cbq_paired' 'cbq_dark2' 'cbq_accent' 'cbq_pastel1' 'cbq_pastel2' 'cbq_set2' 'cbq_set3' 'parula'  ...
    'lut' 'test1' 'test2' 'lightgray' ...
    'custom'};

if isempty(map);
    map = ' ';
end

if nargin<3
    custom = [];
end

switch lower(map)
    case('rb1')
        dat = [0  0  1
               1  1  1
               1  0  0];
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx =0:1/(depth-1):1;
            ny = interp1([0 .5 1],dat(:,ii),nx,'linear');
            cm(:,ii) = ny;
        end
        cm = reshape(cm,256,1,3);
    case('rb2')
        dat = [0  0  1
            0  0  0
            1  0  0];
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx =0:1/(depth-1):1;
            ny = interp1([0 .5 1],dat(:,ii),nx,'linear');
            cm(:,ii) = ny;
        end
        cm = reshape(cm,256,1,3);
    case('hc1')
        dat = [ 0  0  1;
                0 .5  1;
                0  1  1;
                1  1  1;
                1  1  1;
                1  1  0;
                1 .5  0;
                1  0  0];
        nx = (0:255);
        x = linspace(0,255,size(dat,1))';
        cm = [];
        cm(:,1) = interp1(x,dat(:,1),nx,'linear');
        cm(:,2) = interp1(x,dat(:,2),nx,'linear');
        cm(:,3) = interp1(x,dat(:,3),nx,'linear');
        cm = cm.^.75;
    case('hc2');
        dat = [ 1  0  0;
                1 .5  0;
                1  1  0;
                1  1  1;
                0  1  1;
                0 .5  1;
                0  0  1];

        nx = (0:255);
        x = linspace(0,255,size(dat,1))';
        cm = [];
        cm(:,1) = interp1(x,dat(:,1),nx,'linear');
        cm(:,2) = interp1(x,dat(:,2),nx,'linear');
        cm(:,3) = interp1(x,dat(:,3),nx,'linear');
        cm = cm.^.75;
    case('gray')
        cm = gray(depth);
    case('hot')
        cm = hot(depth);
    case('hot3')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(1,1,depth)';
        cm(:,2) = linspace(0,1,depth)';
        cm(:,3) = linspace(0,0,depth)';  
    case('cool')
        cm = cool(depth);
    case('cool2')
        cm = cool(depth);    
        cm = cm(end:-1:1,:);
     case('cool3')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(0,0,depth)';
        cm(:,2) = linspace(0,1,depth)';
        cm(:,3) = linspace(1,1,depth)';  
    case('jet')
        cm = [0 0 0; jet(depth-1)];
    case('jet1')
        cm = jet(depth*2);
        cm = cm(1:depth,:);
    case('jet2')
        cm = jet(depth*2);    
        cm = cm(end:-1:end-(depth-1),:);
    case('hsv')
        cm = hsv(depth);
    case('pink')
        cm = pink(depth);
    case('spring')
        cm = spring(depth);
    case('summer')
        cm = summer(depth);
    case('autumn')
        cm = autumn(depth);
    case('winter')
        cm = winter(depth);
    case('bone')
        cm = bone(depth);
    case('copper')
        cm = copper(depth);
    case('parula')
        cm = parula(depth);
    case('red')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,1) = x(end:-1:1);
    case('red2')
        cm = zeros(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,1) = x(1:end);
    case('green')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,2) = x(end:-1:1);
    case('blue')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,3) = x(1:end);
    case('lightblue')
        cm = reshape([linspace(.25,.8,depth)' linspace(.25,.8,depth)' ones(depth,1)],depth,1,3);
    case('yellow')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,1) = x(end:-1:1);
        cm(:,2) = x(end:-1:1);
    case('yellow2')
        cm = zeros(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,1) = x(1:end);
        cm(:,2) = x(1:end);
    case('purple')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,1) = x(end:-1:1);
        cm(:,3) = x(end:-1:1);
    case('purple2')
        cm = zeros(depth,3);
        x = 0:(1/(depth-1)):1;
        cm(end:-1:1,1) = x(end:-1:1);
        cm(end:-1:1,3) = x(end:-1:1);
    case('cyan2')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(0,0,depth)';
        cm(:,2) = linspace(0,1,depth)';
        cm(:,3) = linspace(0,1,depth)';
    case('mag2')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(0,1,depth)';
        cm(:,2) = linspace(0,0,depth)';
        cm(:,3) = linspace(0,1,depth)';
    case('ocean')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,2) = x(end:-1:1);
        cm(:,3) = x(end:-1:1);
    case('orange')
        cm = ones(depth,3)*0;
        x = 0:(1/(depth-1)):1;
        cm(:,1) = x(end:-1:1);
        cm(:,2) = x(end:-1:1)*.5;
    case('my')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(1,1,depth)';
        cm(:,2) = linspace(0,1,depth)';
        cm(:,3) = linspace(1,0,depth)';
    case('ym')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(1,1,depth)';
        cm(:,2) = linspace(1,0,depth)';
        cm(:,3) = linspace(0,1,depth)';
    case 'mby'
        n = depth/2;
        cm = zeros(depth,3)*0;
        
        cm(1:n,1) = linspace(1,0,n)';
        cm(1:n,2) = linspace(0,0,n)';
        cm(1:n,3) = linspace(1,0,n)';
        
        cm(n+1:end,1) = linspace(0,1,n)';
        cm(n+1:end,2) = linspace(0,1,n)';
        cm(n+1:end,3) = linspace(0,0,n)';
    case 'ybm'
        n = depth/2;
        cm = zeros(depth,3)*0;
        
        cm(1:n,1) = linspace(1,0,n)';
        cm(1:n,2) = linspace(0,0,n)';
        cm(1:n,3) = linspace(1,0,n)';
        
        cm(n+1:end,1) = linspace(0,1,n)';
        cm(n+1:end,2) = linspace(0,1,n)';
        cm(n+1:end,3) = linspace(0,0,n)';
        
        cm = cm(end:-1:1,:);
            
    case('cy')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(0,1,depth)';
        cm(:,2) = linspace(1,1,depth)';
        cm(:,3) = linspace(1,0,depth)';
    case('yc')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(1,0,depth)';
        cm(:,2) = linspace(1,1,depth)';
        cm(:,3) = linspace(0,1,depth)';      
    case('mc')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(1,0,depth)';
        cm(:,2) = linspace(0,1,depth)';
        cm(:,3) = linspace(1,1,depth)';
    case('cm')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(0,1,depth)';
        cm(:,2) = linspace(1,0,depth)';
        cm(:,3) = linspace(1,1,depth)';
    case('lightgray')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(.2,.8,depth)';
        cm(:,2) = linspace(.2,.8,depth)';
        cm(:,3) = linspace(.2,.8,depth)';
    case('lightergray')
        cm = zeros(depth,3)*0;
        cm(:,1) = linspace(.3,.7,depth)';
        cm(:,2) = linspace(.3,.7,depth)';
        cm(:,3) = linspace(.3,.7,depth)';
    case('cbd_rd_yl_gn')
        cols = [0.6471 0      0.1490
                1.0000 1.0000 0.7490
                0      0.4078 0.2157];
        cm = mkmap(cols,depth);  
    case('cbd_rd_yl_bu')
        cols = [0.6471         0    0.1490
                1.0000    1.0000    0.7490
                0.1922    0.2118    0.5843];
        cm = mkmap(cols,depth);
    case('cbd_rd_gy')
        cols = [0.4039         0    0.1216
                1.0000    1.0000    1.0000
                0.1020    0.1020    0.1020];
        cm = mkmap(cols,depth);
    case('cbd_rd_bu')
        cols = [0.4039 0      0.1216
                1      1      1
                0.0196 0.1882 0.3804];
        cm = mkmap(cols,depth);
    case('cbd_pu_or')
        cols = [0.4980    0.2314    0.0314
                0.9686    0.9686    0.9686
                0.1765         0    0.2941];
        cm = mkmap(cols,depth);
    case('cbd_pu_gn')
        cols = [0.5569 0.0039 0.3216; 
                1      1      1     ;
                0.1529 0.3922 0.0980];
        cm = mkmap(cols,depth);   
    case('cbd_pi_yl_gn')
        cols = [0.5569    0.0039    0.3216
                0.9686    0.9686    0.9686
                0.1529    0.3922    0.0980];
        cm = mkmap(cols,depth);    
    case('cbd_br_bg')
        cols = [0.3294    0.1882    0.0196
                0.9608    0.9608    0.9608
                0    0.2353    0.1882];
        cm = mkmap(cols,depth);  
    case('cbs_yl_or_rd')    
        cols = [1.0000    1.0000    0.8000
                0.9882    0.3059    0.1647
                0.5020         0    0.1490];
         cm = mkmap(cols,depth); 
    case('cbs_yl_or_br')    
        cols = [1.0000    1.0000    0.8980
                0.9255    0.4392    0.0784
                0.4000    0.1451    0.0235];
         cm = mkmap(cols,depth); 
    case('cbs_yl_gn_bu')    
        cols = [1.0000    1.0000    0.8510
                0.1137    0.5686    0.7529
                0.0314    0.1137    0.3451];
         cm = mkmap(cols,depth);
    case('cbs_yl_gn')    
        cols = [1.0000    1.0000    0.8980
                0.2549    0.6706    0.3647
                0    0.2706    0.1608];
        cm = mkmap(cols,depth);    
    case('cbs_rd_pu') 
        cols = [1.0000    0.9686    0.9529
                0.8667    0.2039    0.5922
                0.2863         0    0.4157];
        cm = mkmap(cols,depth);
    case('cbs_pu') 
        cols = [0.9882    0.9843    0.9922
                0.5020    0.4902    0.7294
                0.2471         0    0.4902];
        cm = mkmap(cols,depth);
    case('cbs_pu_rd') 
        cols = [0.9686    0.9569    0.9765
                0.9059    0.1608    0.5412
                0.4039         0    0.1216];
        cm = mkmap(cols,depth);
    case('cbs_pu_bu_gn') 
        cols = [1.0000    0.9686    0.9843
                0.2118    0.5647    0.7529
                0.0039    0.2745    0.2118];
        cm = mkmap(cols,depth);
    case('cbs_pu_bu')         
        cols = [1.0000    0.9686    0.9843
                0.2118    0.5647    0.7529
                0.0078    0.2196    0.3451];
        cm = mkmap(cols,depth);
    case('cbs_or_rd') 
        cols = [1.0000    0.9686    0.9255
                0.9373    0.3961    0.2824
                0.4980         0         0];
        cm = mkmap(cols,depth);
    case('cbs_or') 
        cols = [1.0000    0.9608    0.9216
                0.9451    0.4118    0.0745
                0.4980    0.1529    0.0157];
        cm = mkmap(cols,depth);
    case('cbs_gr') 
        cols = [1.0000    1.0000    1.0000
                0.4510    0.4510    0.4510
                0         0         0];
        cm = mkmap(cols,depth);
    case('cbs_gn') 
        cols = [0.9686    0.9882    0.9608
                0.2549    0.6706    0.3647
                0    0.2667    0.1059];
        cm = mkmap(cols,depth);
    case('cbs_gn_bu') 
        cols = [0.9686    0.9882    0.9412
                0.3059    0.7020    0.8275
                0.0314    0.2510    0.5059];
        cm = mkmap(cols,depth);
    case('cbs_bu_pu') 
        cols = [0.9686    0.9882    0.9922
                0.5490    0.4196    0.6941
                0.3020         0    0.2941];
        cm = mkmap(cols,depth);
    case('cbs_bn_gn') 
        cols = [0.9686    0.9882    0.9922
                0.2549    0.6824    0.4627
                0    0.2667    0.1059];
        cm = mkmap(cols,depth);
    case('cbs_bu')         
        cols = [0.9686    0.9843    1.0000
                0.2588    0.5725    0.7765
                0.0314    0.1882    0.4196];
        cm = mkmap(cols,depth);
    case('cbq_set1')
        dat = [0.8941    0.1020    0.1098
               0.2157    0.4941    0.7216
               0.3020    0.6863    0.2902
               0.5961    0.3059    0.6392
               1.0000    0.4980         0
               1.0000    1.0000    0.2000
               0.6510    0.3373    0.1569
               0.9686    0.5059    0.7490
               0.6000    0.6000    0.6000];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_paired')
        dat = [0.6510    0.8078    0.8902
                0.1216    0.4706    0.7059
                0.6980    0.8745    0.5412
                0.2000    0.6275    0.1725
                0.9843    0.6039    0.6000
                0.8902    0.1020    0.1098
                0.9922    0.7490    0.4353
                1.0000    0.4980         0
                0.7922    0.6980    0.8392];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_dark2')    
        dat = [0.1059    0.6196    0.4667
            0.8510    0.3725    0.0078
            0.4588    0.4392    0.7020
            0.9059    0.1608    0.5412
            0.6549    0.3961    0.3176
            0.4000    0.6510    0.1176
            0.9020    0.6706    0.0078
            0.6510    0.4627    0.1137
            0.4000    0.4000    0.4000];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_accent')
        dat = [0.4980    0.7882    0.4980
            0.7451    0.6824    0.8314
            0.9922    0.7529    0.5255
            1.0000    1.0000    0.6000
            0.6118    0.8000    0.6588
            0.2196    0.4235    0.6902
            0.9412    0.0078    0.4980
            0.7490    0.3569    0.0902
            0.4000    0.4000    0.4000];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_pastel1')
        dat = [0.9843    0.7059    0.6824
            0.7020    0.8039    0.8902
            0.8000    0.9216    0.7725
            0.8706    0.7961    0.8941
            0.9961    0.8510    0.6510
            1.0000    1.0000    0.8000
            0.8980    0.8471    0.7412
            0.9922    0.8549    0.9255
            0.9490    0.9490    0.9490];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_pastel2')
        dat = [0.7020    0.8863    0.8039
                0.9922    0.8039    0.6745
                0.7961    0.8353    0.9098
                0.9569    0.7922    0.8941
                0.9294    0.8784    0.8549
                0.9020    0.9608    0.7882
                1.0000    0.9490    0.6824
                0.9451    0.8863    0.8000
                0.8000    0.8000    0.8000];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_set2')
        dat = [0.4000    0.7608    0.6471
            0.9882    0.5529    0.3843
            0.5529    0.6275    0.7961
            0.9059    0.5412    0.7647
            0.7804    0.6941    0.5765
            0.6510    0.8471    0.3294
            1.0000    0.8510    0.1843
            0.8980    0.7686    0.5804
            0.7020    0.7020    0.7020];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('cbq_set3')
        dat = [0.5529    0.8275    0.7804
            1.0000    1.0000    0.7020
            0.7451    0.7294    0.8549
            0.9843    0.5020    0.4471
            0.5020    0.6941    0.8275
            0.9922    0.7059    0.3843
            0.7020    0.8706    0.4118
            0.9882    0.8039    0.8980
            0.8510    0.8510    0.8510];
        
        cm = zeros(depth,3);
        for ii = 1:size(dat,2);
            nx = 1:8/(depth-1):9;
            ny = interp1(1:9,dat(:,ii),nx,'PCHIP');
            cm(:,ii) = ny;
        end
    case('lut')
        m = ReadInFile('/autofs/space/rincewind_003/users/NewAtlas/clut.txt',' ',1);
        list = [m{:,1}]';
        cm = zeros(numel(min(list):max(list)),3);
        for ii = min(list):max(list)
            j1 = find(list==ii);
            if ~isempty(j1)
                cm(ii+1,1:3) = [m{j1,3} m{j1,4} m{j1,5}]/255;
            end
        end
    case('sunrise')
        cc = [147 012 006;
            174 060 011;
            201 108 016;
            228 157 021;
            255 205 026];
        cc = cc/255;
        
        dat = cc;
        % dat(:,2) = dat(:,2).^1.3
        nx = (0:depth-1);
        x = linspace(0,depth,size(dat,1))';
        cm = [];
        cm(:,1) = interp1(x,dat(:,1),nx,'linear');
        cm(:,2) = interp1(x,dat(:,2),nx,'linear');
        cm(:,3) = interp1(x,dat(:,3),nx,'linear');
        % cm = cm(end:-1:1,:).^.75;
        cm = cm.^.75;
        
        %image(1:256,1, reshape(cm,256,1,3)); shg

    case('test1')
        %figure(10); clf; imagesc((1:255)'); set(gca,'ydir','normal'); colorbar; shg
        xx = (1:depth)';
        cols = [0 0 0; 1 0 0; 1 1 0; 0 1 1; 0 0 1; 1 0 1;];
        %%cols = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
        x = (1:(depth-1)/(size(cols,1)-1):depth)';
        for ii = 1:size(cols,2)
            cm(:,ii) = spline(x, cols(:,ii),xx);
        end
        cm(cm<0)=0; cm(cm>1)=1;
        
%         for ii = 1:size(cm,2);
%             if range(cm(:,ii))>1
%                 cm(:,1) = cm(:,ii)-min(cm(:,ii));
%                 cm(:,ii) = cm(:,ii)./max(cm(:,ii));
%             end
%         end
    case('test2')
        %figure(10); clf; imagesc((1:255)'); set(gca,'ydir','normal'); colorbar; shg
        xx = (1:depth)';
        cols = [.75  0 .75; .3 .3 .9; .3 .9 .3; .75 .75 0; .7 .4 .4; .9 .3 .3];
        x = (1:(depth-1)/(size(cols,1)-1):depth)';
        for ii = 1:size(cols,2)
            cm(:,ii) = interp1(x, cols(:,ii),xx,'linear');
        end
        cm(cm<0)=0; cm(cm>1)=1;
    case('custom')
        if isempty(custom)
            fh = figure(777); clf; %reset(777);
            
            pop(1) = uicontrol(fh,'style','text','Units','Normalized','position',[.05 .83 .9 .05],'String','Enter the colors you want (e.g. [1 0 0; 0 1 0; 0 0 1]','Fontsize',12);
            pop(2) = uicontrol(fh,'style','edit','Units','Normalized','position',[.05 .75 .9 .08],'String','','Fontsize',12,'Callback','uiresume');
            
            uiwait
            eval(['custom = [' get(pop(2),'String') '];']);
            close(gcf);
        end
        cm = mkmap(custom,depth);
    case 'rb'
        cols = [1 0 0
                0 0 0
                0 0 1];
        cm = mkmap(cols,depth);
    otherwise
        cm = jet(depth);
end

end

function cm = mkmap(cols,depth)
cm = zeros(depth,3);
if size(cols,1)==2;
    cm(:,1) = linspace(cols(1,1),cols(2,1),depth)';
    cm(:,2) = linspace(cols(1,2),cols(2,2),depth)';
    cm(:,3) = linspace(cols(1,3),cols(2,3),depth)';
elseif size(cols,1)==3;
    n = depth/2;
    cm(1:n,1) = linspace(cols(1,1),cols(2,1),n)';
    cm(1:n,2) = linspace(cols(1,2),cols(2,2),n)';
    cm(1:n,3) = linspace(cols(1,3),cols(2,3),n)';
    
    cm(n+1:end,1) = linspace(cols(2,1),cols(3,1),n)';
    cm(n+1:end,2) = linspace(cols(2,2),cols(3,2),n)';
    cm(n+1:end,3) = linspace(cols(2,3),cols(3,3),n)';
end
end
% figure(1); clf; scatter(1:256,1:256,[],cm);
