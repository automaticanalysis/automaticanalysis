function GLM_Plot_spm(VOI,figno) 
if nargin==1
    figno = 999;
end

if ~isfield(VOI,'ConSpec')
    VOI.ConSpec.Levs = 1;
end

PlotDescrip = [];
PlotType = [];

y = VOI.data;
X = VOI.DM;
DM = X;
X(isnan(X))=0;
con = VOI.con;

labs = VOI.SPM.xX.name;

in = [];
for ii = 1:size(con,2);
    in = [in; find(con(:,ii)==max(con(:,ii)) & con(:,ii)~=0)];
    in = [in; find(con(:,ii)==min(con(:,ii)) & con(:,ii)~=0)];
end
in = unique(in);
out = setdiff(1:size(con,1),in)';

c = 0;
inds = [];
for ii = 1:numel(in)
    for jj = ii+1:numel(in);
        c = c+1;
        inds(c) = numel(intersect(find((DM(:,in(ii)))~=0), find((DM(:,in(jj)))~=0)));
    end
end

if all(inds==0) && ~isempty(inds);
    PlotDescrip = [PlotDescrip 'Group'];
else
    PlotDescrip = [PlotDescrip 'All'];
end

SubTerms = contains('^[S|s]ubject.*', labs);
consts = find(mean(X==1 | X==0)==1);
consts = setdiff(consts,SubTerms);
cout = setdiff(SubTerms,find(abs(con)~=0));
GroupTerms = setdiff(consts,SubTerms);

%[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end

if numel(find(ismember(consts,in)==1))==numel(in)
    if numel(VOI.ConSpec.Levs)==1
        PlotDescrip = [PlotDescrip 'Means'];
    elseif numel(VOI.ConSpec.Levs)==2
        PlotDescrip = [PlotDescrip 'Int'];
    elseif numel(VOI.ConSpec.Levs)>2
        error('This Script cannot deal with 3+ way interactions');
    end
else
    if numel(find(VOI.con<0))==0 && size(VOI.con,2)>1
        PlotDescrip = [PlotDescrip 'ContinuousPred'];
    else
        PlotDescrip = [PlotDescrip 'ContinuousVars'];
        if numel(GroupTerms)>1 && numel(in)==1 %% Need to see about generalizing this, but for now it should work just fine for ploting a single covariate controlling for multiple group means.
            PlotDescrip = [PlotDescrip 'ContinuousVars_ControlForGroupMeans'];
        end
    end
end

check = setdiff(1:size(X,2),unique([consts(:)' SubTerms(:)' in(:)']));
ortho = [];
for ii = 1:numel(check)
    c = 0;
    for jj = 1:numel(in)
        c = c+1;
        inds(c) = numel(intersect(find(~isnan(DM(:,in(jj)))), find(~isnan(DM(:,check(ii))))));
    end
    if all(inds==0)
        ortho(end+1) = check(ii);
    end
end

if numel(consts) ~= size(X,2)-numel(SubTerms)-numel(setdiff(in,consts))-numel(ortho)
    if numel(cout)>0
        PlotDescrip = [PlotDescrip '_SubResCovarRes'];
    else
        PlotDescrip = [PlotDescrip '_CovarRes'];
    end
else
    if numel(cout)>0
        PlotDescrip = [PlotDescrip '_SubRes'];
    else
        PlotDescrip = [PlotDescrip '_NoRes'];
    end
end

if ~isempty(contains('_SubResCovarRes',{PlotDescrip}))
    error('This implies a repeated measures ANCOVA. This is not a configured option.');
end
% PlotDescrip
if ~isempty(contains('_CovarRes',{PlotDescrip}))
    yres = zeros(size(X,1),numel(in));
    xres = zeros(size(X,1),numel(in));
    
    for ii = 1:length(in);
        tout = setdiff(1:size(X,2),in(ii));
        tin = in(ii);
        
        tx = X(:,consts);
        for jj = setdiff(tout,consts)
            tx(:,end+1) = X(:,jj);
            if mean(X(:,jj))~=0
                warning('Your covariate is not demeaned!');
            end
        end
        
        b = pinv(tx)*y;
        pred = tx*b;
        C = X(:,consts)*b(find(ismember(tout,consts)==1),:);
        yres(:,ii) = y-pred+sum(C,2);
        
        b = pinv(tx)*X(:,tin);
        pred = tx*b;
        C = X(:,consts)*b(find(ismember(tout,consts)==1),:);
        xres(:,ii) = X(:,tin)-pred+sum(C,2);
    end
elseif ~isempty(contains('_SubRes',{PlotDescrip}))
    
    nx = X(:,SubTerms);
    b = pinv(nx)*y;
    pred = nx*b;
    C = X(:,in)*(pinv(X(:,in))*y);
    yres = y-pred+C;
    xres = X(:,in);
    
% %     [out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
% 
%     subterm = [unique([FacNums{in}]) numel(size(VOI.F.FF))];
%     cellNum = sub2ind_aps(size(VOI.F.FF),subterm);
%     tmp = VOI.F.FF{cellNum}; tmp(isnan(tmp))=0; 
%     yres = pinv(tmp)*yres;
%     %ni = find(sum(abs(diff(X(:,in))),2)~=0); size(X,1);
%     [a ni] = unique(tmp,'rows'); ni = sort(ni);
%     [ri,ci] = find(tmp(ni,:)==1); [trash,yord] = sort(ri);
%     yres = yres(yord);
%     xres = X(ni,:);
else
    yres = repmat(y,1,numel(in));
    xres = X(:,in);
end

if contains('ContinuousVars_ControlForGroupMeans',{PlotDescrip})
    %[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
    C = X(:,GroupTerms)*(pinv(X(:,GroupTerms))*y);
    yres = y-C;
    
    C = X(:,GroupTerms)*(pinv(X(:,GroupTerms))*X(:,in));
    xres = X(:,in)-C;
end

% if contains('ContinuousPred',{PlotDescrip})

if contains('Continuous',{PlotDescrip})
    PlotType = 'Scatter';
end
if contains('ContinuousPred',{PlotDescrip})
    PlotType = 'Pred';
end
if contains('^GroupMeans',{PlotDescrip})
    PlotType = 'Group';
    if size(VOI.con,2)==1 && numel(in) == 4
        if sum(VOI.con(in)-[1 -1 -1 1]') == 0;
            tc = zeros(numel(con),4);
            for ii = 1:numel(in)
                tc(in(ii),ii)=1;
            end
            con = tc;
        end
    end
end

if contains('AllMeans',{PlotDescrip})
    PlotType = 'Group';
end
if contains('GroupInt',{PlotDescrip})
    PlotType = 'IntPlot';
end
% [out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
switch PlotType
    case 'Group'
        %keyboard
        %%%
        Groups = zeros(size(X,1),1);
        tr = []; c = 0; ll = [];
        for ii = 1:size(con,2)
            a = find(con(:,ii)==max(con(:,ii)));
            if ~any(ismember(a,tr))
                tr = [tr; a];
                ind = []; tmp = [];
                for jj = 1:numel(a);
                    ind = [ind;  find(DM(:,a(jj))==1)];
                    tmp = [tmp labs{a(jj)} ' & '];
                end
                c = c+1;
                ll{c} = tmp(1:end-3);
                Groups(ind)=c;
            end
            
            a = find(con(:,ii)==min(con(:,ii)) & con(:,ii)>0);
            if ~any(ismember(a,tr)) && ~isempty(a);
                tr = [tr; a];
                ind = []; tmp = [];
                for jj = 1:numel(a);
                    ind = [ind;  find(DM(:,a(jj))==1)];
                    tmp = [tmp labs{a(jj)} ' & '];
                end
                c = c+1;
                ll{c} = tmp(1:end-3);
                Groups(ind)=c;
            end
        end
        %%%
        ind = find(Groups~=0);
        %keyboard;
        if numel(in)==1
            %[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
            figure(figno); clf; reset(figno); set(figno,'name',PlotDescrip);
            ScatterGroups(yres(ind(Groups==1)),Groups(ind(Groups==1)),ll);
        else
            if exist('boxplot.m','file')>0
                figure(figno); clf; reset(figno); set(figno,'name','GroupMeans_NoRes');
                boxplot(y(ind),Groups(ind),'notch','on'); hold on;
                set(gca,'XTick', 1:numel(in), 'XTickLabel',ll,'FontSize',14);
                ax = axis;
                plot(ax(1:2),[0 0],'m--','linewidth',2);
                    
                if ~strcmpi(PlotDescrip,'GroupMeans_NoRes')
                    figure(figno+101); clf; reset(figno+101); set(figno+101,'name',PlotDescrip);
                    if exist('ni','var')
                        [trash yind] = intersect(ni,ind);
                        ind = ind(ni);
                    else
                        yind = ind;
                    end
                    
                    boxplot(yres(yind),Groups(ind),'notch','on'); hold on;
                    set(gca,'XTick', 1:numel(in), 'XTickLabel',ll,'FontSize',14);
                    ax = axis;
                    plot(ax(1:2),[0 0],'m--','linewidth',2);
                end
            else
                figure(figno); clf; reset(figno); set(figno,'name',PlotDescrip);
                ScatterGroups(yres(ind),Groups(ind),ll);
            end
        end
    case 'Pred'
        figure(figno); clf; reset(figno); set(figno,'name',PlotDescrip);
        i1 = ind;
        plot(xres(i1,1),yres(i1,1),'d','markersize',8,'linewidth',2); hold on;
        r = num2str(corr(xres(i1,1),yres(i1,1)));
        
        [row col] = find(con==1);
        ll = [];
        for ii = 1:numel(row);
            ll = [ll labs{row(ii)}];
            if ii<numel(row); ll = [ll ' & ']; end
        end
        info{ii} = [ll ': r = ' r(1:end-3)];
        
        
        title(info{ii},'fontsize',16);
        axis square;
        yfit = polyval(polyfit(xres(i1,1),yres(i1,1),1),xres(i1,1));
        [a,ord] = sort(xres(i1,1));
        plot(a,yfit(ord),'-','linewidth',2);
    case 'Scatter'
        if contains('ContinuousVars_ControlForGroupMeans',{PlotDescrip})
            figure(figno); clf; reset(figno); set(figno,'name','ContinuousVars_DiffGroups_NoGroupMeanRes');
            linespec = {'b*' 'rd' 'ks' 'm+' 'cx' 'g^' 'yv'};
            %[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end

            rows = find(XX(:,in(1))~=0); %% This will fail when 0 is part of a covariate.  Might try doing something with GroupTerms
            for ii = 1:numel(GroupTerms);
                i1 = find(X(rows,GroupTerms(ii))==1);
                if isempty(i1)
                   GroupTerms = setdiff(GroupTerms,GroupTerms(ii));
                   continue
                end
                plot(xres(i1), y(i1),linespec{ii},'markersize',8); hold on;                
            end
            legend(labs(GroupTerms),'FontSize',14,'location','best');
            
            [a ord] = sort(xres(rows));
            yfit = polyval(polyfit(xres(rows(ord)),y(rows(ord)),1),xres(rows(ord)));
            plot(xres(ord),yfit,'k--');
            
            r = num2str(corr(xres(rows(ord)),y(rows(ord))));
            info = [labs{in(1)} ': r = ' r(1:end-3)];
            
            title(info,'fontsize',16);
            axis square;
            
            figure(figno+100); clf; reset(figno+100); set(figno+100,'name','ContinuousVars_MinusGroupMeans');
            plot(xres(rows),yres(rows),'*'); hold on;
            [a ord] = sort(xres(rows));
            yfit = polyval(polyfit(xres(rows(ord)),yres(rows(ord)),1),xres(rows(ord)));
            plot(xres(ord),yfit,'k--');
            
            r = num2str(corr(xres(rows(ord)),yres(rows(ord))));
            info = [labs{in(1)} ': r = ' r(1:end-3)];
            
            title(info,'fontsize',16);
            axis square;
        else
            figure(figno); clf; reset(figno); set(figno,'name',PlotDescrip);
            for ii = 1:numel(in);
                i1 = find(~isnan(DM(:,in(ii))));
                
                if numel(in)>2
                    subplot(1,numel(in),ii);
                end
                
                if numel(in) == 2 && ii ==2
                    plot(xres(i1,ii),yres(i1,ii),'rs','markersize',8,'linewidth',2); hold on;
                else
                    plot(xres(i1,ii),yres(i1,ii),'d','markersize',8,'linewidth',2); hold on;
                end
                
                r = num2str(corr(xres(i1,ii),yres(i1,ii)));
                info{ii} = [labs{in(ii)} ': r = ' r(1:end-3)];
                
                
                if numel(in)>2
                    title(info{ii},'fontsize',16);
                end
                axis square;
            end
            if numel(in)<=2
                legend(gca,info,'location','best');
            end
            for ii = 1:numel(in);
                i1 = find(~isnan(DM(:,in(ii))));
                
                if numel(in)>2
                    subplot(1,numel(in),ii);
                end
                
                if numel(in) == 2 && ii ==2
                    yfit = polyval(polyfit(xres(i1,ii),yres(i1,ii),1),xres(i1,ii));
                    [a,ord] = sort(xres(i1,ii));
                    plot(a,yfit(ord),'r-', 'linewidth',2);
                else
                    yfit = polyval(polyfit(xres(i1,ii),yres(i1,ii),1),xres(i1,ii));
                    [a,ord] = sort(xres(i1,ii));
                    plot(a,yfit(ord),'-','linewidth',2);
                end
                
                plot(X(i1,in(ii)),y(i1),'k.','markersize',10); hold on;
                %r = num2str(corr(X(i1,in(ii)),y(i1)));
                %info{2} = [labs{in(ii)} ': r = ' r(1:end-3)];
                
                yfit = polyval(polyfit(X(i1,in(ii)),y(i1),1),X(i1,in(ii)));
                [a,ord] = sort(X(i1,in(ii)));
                plot(a,yfit(ord),'k--','linewidth',.5);
            end
        end
    case 'IntPlot'
        warning('Interaction plots are not configured for second level spm models');
        
    otherwise
        warning('This option has not been configured yet');
end



