function [yres xres part] = GLM_Plot_Fast2(y,mod,effect)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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

%if searchCellStr('aschultz', {UserTime}); keyboard; end

yres = [];
xres = [];
part = [];

flag = 0;

effect = regexprep(effect,'\*',':');
effect = regexprep(effect,'_x_',':');

br = 0;
for ii = 1:numel(mod.RFMs)
    if br==1; break; end
    for jj = 1:numel(mod.RFMs(ii).Effect)
        if strcmpi(effect, mod.RFMs(ii).Effect(jj).name)
            which = [ii jj];
            part = mod.RFMs(ii).Effect(jj);
            br = 1;
            break;
        end
    end
end
 
if strcmpi('FullModel',effect)
    type1 = 'Continuous';
    type2 = 'Prediction';
else

    if isempty(part.dif) %isfield(part,'pieces') && ~isempty(part.pieces)
        type1 = 'Continuous';
        type2 = 'Prediction';
        flag = 1;
    else
        
        terms = regexp(effect,'[\*\:]','split');
        type1 = 'Categorical'; iscat = [];
        for ii = 1:numel(terms)
            if isnumeric(mod.data.(terms{ii}));
                type1 = 'Continuous';
                iscat(ii) = 0;
            else
                iscat(ii) = 1;
            end
        end
        
        type2 = [];
        if isempty(regexp(effect,'[\*\:]'))
            type2 = 'MainEffect';
        else
            type2 = 'Interaction';
        end
    end
end


if strcmpi(type1,'Categorical') && strcmpi(type2,'MainEffect');
    figure(gcf); clf; set(gcf,'Name','Raw Data Plot');
    ScatterGroups2(y,mod.data.(effect));
    
%     if exist('part','var');
%         yres = crtlFor(y,part.tx2);%+mean(y);
%         figure(gcf+1); clf; set(gcf,'Name','Type III Residualized Plot');
%         ScatterGroups2(yres,mod.data.(effect));
%     end
end

if strcmpi(type1,'Categorical') && strcmpi(type2,'Interaction');
    figure(gcf); clf; set(gcf,'Name','Raw Data Plot');
    F = {};
    for ii = 1:numel(terms)
        tmp = mod.data.(terms{ii});
        if isobject(tmp)
            F(:,ii) = cellstr(tmp);
        else
            F(:,ii) = tmp;
        end
        
    end
    InterPlot2(y,F);
    
%     if exist('part','var');
%         figure(gcf+1); clf; set(gcf,'Name','Type III Residualized Plot');
%         yres = crtlFor(y,part.tx2);%+mean(y);
%         InterPlot2(yres,F);
%     end
    
end

if strcmpi(type1,'Continuous') && strcmpi(type2,'MainEffect');
    figure(gcf); clf; set(gcf,'Name','Raw Data Plot');    
    ScatterPlot2(mod.data.(effect),y,ones(size(y)),{part.name 'Y'});

%     if exist('part','var');
%         figure(gcf+1); clf; set(gcf,'Name','Type III Residualized Plot');
%         yres = crtlFor(y,part.tx2);%+mean(y);
%         xres = crtlFor(mod.data.(effect),part.tx2);%+mean(y);
%         ScatterPlot2(xres,yres,ones(size(y)),{part.name 'Y'});
%     end
end

if strcmpi(type1,'Continuous') && strcmpi(type2,'Interaction');
    eff = effect;
    i1 = find(iscat==1);
    F = []; Levs = [];
    for ii = 1:numel(i1);
        F{ii} = makedummy(mod.data.(terms{i1(ii)}),0);
        tmp = unique(mod.data.(terms{i1(ii)}),'stable');
        if isobject(tmp); tmp = cellstr(tmp); end
        Levs{ii} = tmp;
        eff = regexprep(eff,[terms{i1(ii)} ':'],'');
    end
   
    if ~isempty(F)
        F = FactorCross(F);
        L = LabelCross(Levs);
        G = cell(size(F,1),1);
        for ii = 1:size(F,2);
            G(find(F(:,ii)==1)) = L(ii);
        end
    else
        G = ones(size(y));
    end
    
    i1 = find(iscat==0);
    x = 1; xlab = [];
    for ii = 1:numel(i1);
        x = x.*demean(mod.data.(terms{i1(ii)}));
        xlab = [xlab terms{i1(ii)} '*'];
    end
    xlab = xlab(1:end-1);
    
    
    figure(gcf); clf; set(gcf,'Name','Raw Data Plot');
    ScatterPlot2(x,y,G,{xlab 'Y'});
    
%     if exist('part','var') && ~isempty(part);
%         figure(gcf+1); clf; set(gcf,'Name','Type III Residualized Plot');
%         yres = crtlFor(y,part.tx2);
%         xres = crtlFor(x,part.tx2);
%         ScatterPlot2(xres,yres,G,{xlab 'Y'});
%     end
end

%%% Always plot the predicted data, this will always work and will always
%%% be accurate
% if strcmpi(type1,'Continuous') && strcmpi(type2,'Prediction');
%     if ~isempty(searchCellStr('aschultz',{UserTime})); keyboard; end    
    
    tx1 = part.tx1;
    tx2 = part.tx2;
    if isempty(tx2);
        tx2 = eye(size(tx1,1));
    end
    
    warning off
    nm = tx1-(tx2*(tx2\tx1));
    nm(nm>-1e-12 & nm<1e-12)=0;
    warning on
    
    yy = crtlFor(y,tx2);
    xx = nm*(pinv(nm)*yy);
    i1 = find(xx~=0);
    %i1 = 1:size(xx,1);
    
    if std(xx)~=0
        figure; clf; set(gcf,'Name','Predicted vs. Observed');
        ScatterPlot2(xx(i1),yy(i1),ones(numel(i1),1),{part.name 'Y'});
        xlabel('Predicted');
        ylabel('Observed');
    end
% end

if flag %% Try plotting the post-hoc contrast
    if numel(regexp(effect,'&','split'))>1
        warning('Not sure how to display a raw data plot for this contrast');
        return
    end
    
    [con ind tr] = ParseEffect(effect,mod.X,mod.data);
    
    
    if all(all(con~=0))
        groups = regexp(effect,'#','split');
        
        xx = []; yy = []; gg = [];
        g = cell(numel(ind),1);
        for ii = 1:size(con,2);
            yy = [yy; y(ind)];
            xx = [xx; con(ind,ii)];
            g(:) = groups(ii);
            gg = [gg; g];
        end
        
        figure(gcf+1); clf; set(gcf,'Name','Post-Hoc Raw Data');
        ScatterPlot2(xx,yy,gg,{'' ''});
        return
    else
        groups = regexp(effect,'#','split');
        if numel(groups)<size(con,2);
            groups = cellstr(unique(ds.(effect)));
        end
        yy = y(ind);
        
        tc = (con(ind,:)~=0)*(1:size(con,2))';
        
        gg = cell(numel(tc),1);
        for ii = unique(tc)'
            gg(tc==ii)=groups(ii);
        end
    end
    
    figure(gcf+1); clf; set(gcf,'Name','Post-Hoc Raw Data');
    if all(all(con==0 | con==1))
        ScatterGroups(yy,tc,groups);
    else
        ScatterPlot2(sum(con(ind,:),2),yy,gg,{'' ''});
    end
end

end

function [con, indices, track1] = ParseEffect(contrast,X,data)
cond4 = [];
extra = [];

indices = [];
flag = [];
step1 = regexprep(contrast,' ', '');
step2 = regexp(step1,'&','split');
for mm = 1:numel(step2);
    step3 = regexp(step2{mm},'#','split');
    
    cond3 = [];
    track1 = {}; track2 = {};
    for jj = 1:numel(step3)
        cond2 = []; XD = [];
        
        step5 = regexp(step3{jj},'\|','split');
        cond = [];
        for kk = 1:numel(step5)
            step6 = regexp(step5{kk},'\$','split');
            track1(end+1,1:2) = step6;
            track2{end+1,1} = step5{kk};
            
            var = data.(step6{1});
            if isnumeric(var)
                cond = [cond demean(var)];
                flag(kk) = 0;
                XD{end+1}=1:numel(var);
            else
                flag(kk) = 1;
                if iscell(var)
                    tmp = zeros(numel(var),1);
                    i1 = strmatch(var,step6{2});
                    tmp(i1)=1;
                else
                    if numel(step6)==1
                        tmp = makedummy(data.(step6{1}));
                        XD{end+1}=1:size(tmp,1);
                    else
                        tmp = data.(step6{1})==step6{2};
                        XD{end+1} = find(data.(step6{1})==step6{2});
                    end
                end
                cond = [cond tmp];
            end
        end
        
        if any(flag)
            if numel(step5)==1
                cond2 = [cond2 cond];
            else
                cond2 = [cond2 prod(cond,2)];
            end
        else
            cond2 = [cond2 cond];
        end
        
        indices{end+1} = intersections(XD);
        
        if all(flag)
            if numel(step5)==1
                cond3 = [cond3 cond2];
            else
                cond3 = [cond3 sum(cond2,2)>0];
            end
        else
            cond3 = [cond3 cond2];
        end
    end
    cond4 = [cond4 cond3];
    
    if ~isempty(searchCellStr('#',step2(mm)))
        extra = [extra mean(cond3,2)];
    end
    
    indices = unions(indices);
end
con = cond4;
end