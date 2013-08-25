function [yres xres] = GLM_Plot(VOI,figno)
yres = [];
xres = [];

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
% keyboard;
[labs FacNums] = Labels(VOI.F);

SubTerms = contains('^Subject.*', labs);
consts = find(mean(X==1 | X==0)==1);
consts = setdiff(consts,SubTerms);
cout = setdiff(SubTerms,find(abs(con)~=0));
GroupTerms = setdiff(consts,[SubTerms (numel(VOI.F.name)+1:size(X,2))]);

%%% Plotting Groups will not work right is the contrasts is something like
%%% [0.4568 0.5432 0 0] where the contrast is a weighted mean of different
%%% groups.  Currently this will be treated as hust plotting Group1
% [junk user] = UserTime; 
% if strcmp(user,'aschultz'); 
%     in = [];
%     for ii = 1:size(con,2);
%         in = [in; find(con(:,ii)>0)];
%         if any(VOI.con<0)
%             in = [in; find(con(:,ii)<0)];
%         end
%     end
%     in = unique(in);
%     if numel(intersect(consts,in))==numel(in)
%         %sum(X(:,in),2)
%     end
%     
%     keyboard;
% end


in = [];
for ii = 1:size(con,2);
    in = [in; find(con(:,ii)==max(con(:,ii)) & con(:,ii)~=0)];
    if any(VOI.con<0)
        in = [in; find(con(:,ii)==min(con(:,ii)) & con(:,ii)~=0)];
    end
end
in = unique(in);
% in = [VOI.ConSpec.Groups{:}];
out = setdiff(1:size(con,1),in)';
% [junk user] = UserTime; if strcmp(user,'aschultz'); keyboard; end

c = 0;
inds = [];
for ii = 1:numel(in)
    for jj = ii+1:numel(in);
        c = c+1;
        inds(c) = numel(intersect(find(~isnan(DM(:,in(ii)))), find(~isnan(DM(:,in(jj))))));
    end
end

if all(inds==0) && ~isempty(inds);
    PlotDescrip = [PlotDescrip 'Group'];
else
    PlotDescrip = [PlotDescrip 'All'];
end
    


if numel(intersect(in,GroupTerms))==numel(in) && strcmpi(PlotDescrip,'All') && numel(in)>1;
    wh = find(sum(X(:,in))==min(sum(X(:,in))));
    in = in(wh);
end

% [junk user] = UserTime; if strcmp(user,'aschultz'); keyboard; end

if numel(find(ismember(consts,in)==1))==numel(in)
    if numel(VOI.ConSpec.Levs)==1 || isempty(VOI.ConSpec.Levs)
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
        if numel(GroupTerms)>1 && numel(in)~=numel(GroupTerms) %% Need to see about generalizing this, but for now it should work just fine for ploting a single covariate controlling for multiple group means.
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
            if mean(X(:,jj))>1e-10
                warning('Your covariate is not demeaned!');
            end
        end
        

        b = pinv(tx)*y;
        pred = tx*b;
        C = X(:,consts)*b(find(ismember(tout,consts)==1),:);
        if ~isempty(contains('ContinuousVars_ControlForGroupMeans',{PlotDescrip}))
            yres(:,ii) = y-pred+mean(sum(C,2));
        else
            yres(:,ii) = y-pred+sum(C,2);
        end
        %[junk user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
        b = pinv(tx)*X(:,tin);
        pred = tx*b;
        C = X(:,consts)*b(find(ismember(tout,consts)==1),:);
        if ~isempty(contains('ContinuousVars_ControlForGroupMeans',{PlotDescrip}))
            xres(:,ii) = X(:,tin)-pred+mean(sum(C,2));
        else
            xres(:,ii) = X(:,tin)-pred+sum(C,2);
        end
    end
elseif ~isempty(contains('_SubRes',{PlotDescrip}))
    
%     b = pinv(X)*y;
%     pred = X(:,SubTerms)*b(SubTerms,:);
%     C = X(:,in)*(b(in,:));
%     yres = y-pred;%;+C;
%     xres = X(:,in);
%     keyboard;

    nx = [ones(size(X,1),1) X(:,SubTerms)];
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

% keyboard;
if ~isempty(contains('ContinuousVars_ControlForGroupMeans',{PlotDescrip})) && isempty(contains('CovarRes',{PlotDescrip}))
    %[out2 user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
    C = X(:,GroupTerms)*(pinv(X(:,GroupTerms))*y);
    yres = yres-repmat(C,1,size(yres,2))+mean(C);
   
    C = X(:,GroupTerms)*(pinv(X(:,GroupTerms))*xres);
    xres = xres-C+repmat(mean(C),size(X,1),1);
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
end
if contains('AllMeans',{PlotDescrip})
    PlotType = 'Group';
end
if contains('GroupInt',{PlotDescrip})
    PlotType = 'IntPlot';
end

switch PlotType
    case 'Group'
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
            
            a = find(con(:,ii)==min(con(:,ii)));
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
        end
        %[out2 user] = UserTime; if strcmp(user,'award'); keyboard; end
        ind = find(Groups~=0);
        if numel(in)==1
            %[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
            figure(figno); clf; set(figno,'name','GroupMeans_NoRes'); %reset(figno);
            ScatterGroups(y(ind(Groups==1)),Groups(ind(Groups==1)),labs(in));
            if ~strcmpi(PlotDescrip,'AllMeans_NoRes')
                figure(figno+100); clf; set(figno+100,'name',PlotDescrip); %reset(figno+100); 
                ScatterGroups(yres(ind(Groups==1)),Groups(ind(Groups==1)),labs(in));
            end
        else
            figure(figno); clf; set(figno,'name','GroupMeans_NoRes'); %reset(figno); 
            ScatterGroups(y(ind),Groups(ind),ll);
            
            ax = axis;
            plot(ax(1:2),[0 0],'m--','linewidth',2);
            
            if ~strcmpi(PlotDescrip,'GroupMeans_NoRes')
                figure(figno+101); clf; set(figno+101,'name',PlotDescrip); %reset(figno+101);
                if exist('ni','var')
                    [trash yind] = intersect(ni,ind);
                    ind = ind(ni);
                else
                    yind = ind;
                end
                %keyboard;
                ScatterGroups(yres(yind),Groups(ind),ll);
                ax = axis;
                plot(ax(1:2),[0 0],'m--','linewidth',2);
            end
        end
    case 'Pred'
        figure(figno); clf; set(figno,'name',PlotDescrip); %reset(figno); 
        %[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
        [row col] = find(con==1);
        ll = [];
        i1 = [];
        for ii = 1:numel(row);
            ll = [ll labs{row(ii)}];
            if ii<numel(row); ll = [ll ' & ']; end
            i1 = [i1;find(~isnan(DM(:,row(ii))))];
        end
        i1 = unique(i1);
        ScatterPlot({xres(i1,1)},{yres(i1,1)},ll,{'X' 'Predicted'});
        
    case 'Scatter'
        if contains('ContinuousVars_ControlForGroupMeans',{PlotDescrip})
            figure(figno); clf; set(figno,'name','ContinuousVars_DiffGroups_NoGroupMeanRes'); %reset(figno);
            %[out user] = UserTime; if strcmp(user,'aschultz'); keyboard; end
            
            X1 = []; Y1 = [];
            X2 = []; Y2 = [];
            for ii = 1:numel(in);
                i1 = find(~isnan(DM(:,in(ii))));
                X1{ii} =  X(i1,in(ii));
                Y1{ii} =  y(i1);
                X2{ii} =  xres(i1,ii);
                Y2{ii} =  yres(i1);
                ll{ii} = labs{in(ii)};
            end
            ScatterPlot(X1,Y1,ll,[]);
            figure(figno+100); clf; set(figno+100,'name','ContinuousVars_MinusGroupMeans'); %reset(figno+100); 
            ScatterPlot(X2,Y2,ll,[]);
        else
            figure(figno); clf; set(figno,'name',PlotDescrip); %reset(figno); 
            X1 = []; Y1 = [];
            for ii = 1:numel(in);
                i1 = find(~isnan(DM(:,in(ii))));
                
                X1{ii} = xres(i1,ii);
                Y1{ii} = yres(i1,ii);
                ll{ii} = labs{in(ii)};
            end
            ScatterPlot(X1,Y1,ll,[]);
        end
    case 'IntPlot'
        conds = labs(cell2mat(VOI.ConSpec.Groups));
        
        A = [];
        B = [];
        for ii = 1:length(conds);
            tmp = regexp(conds{ii},{' by '},'split');
            A{ii,1} = tmp{1}{1};
            B{ii,1} = tmp{1}{2};
        end
        
        [junk thi] = unique(A);
        A = A(sort(thi));
        [junk thi] = unique(B);
        B = B(sort(thi));
        
        dat2 = cell(1,numel(in));
        for ii = 1:length(in)
            ind = find(~isnan(DM(:,in(ii))));
            dat2{ii} = yres(ind);
        end
        
        InterPlot(dat2,{A B},figno);
        
        if contains('Res',{PlotDescrip})
            set(gcf,'name', 'Residualized');
            dat2 = cell(1,numel(in));
            for ii = 1:length(in)
                ind = find(~isnan(DM(:,in(ii))));
                dat2{ii} = y(ind);
            end
            InterPlot(dat2,{A B},figno+100);
            set(gcf,'name', 'NoRes');
        end
        
    otherwise
        warning('This option has not been configured yet');
end

end


function [O fn] = Labels(F)

labs = [F.IN.BetweenLabs F.IN.WithinLabs F.IN.CovarLabs];
O = [];
for zz = 1:length(F.name);
    fn1 = [];
    tmp = F.name{zz};
    
    if ~isempty(contains('^S.*',{tmp}))
        O{zz} = ['Subject' tmp(2:end)];
        continue
    end
    t1 = regexp(tmp,'_|x','split');
    t1 = t1(setdiff(1:numel(t1),contains('F',t1)));
    out2 = [];
    for aa = 1:2:length(t1);
        t2 = str2num(t1{aa});
        fn1(end+1) = t2;
        t3 = str2num(t1{aa+1});
        try
            out2 = [out2 labs{t2}{t3}];
        catch
            out2 = [out2 labs{t2}];
        end
        
        if aa+1<numel(t1)
            out2 = [out2 ' by '];
        end
    end
    
    O{zz} = out2;
    fn{zz} = fn1;
end
end

