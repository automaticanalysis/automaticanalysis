function [ Stats ] = mvpaa_stats1st( aap, uSimil)
% MVPA_STATS Do 1st level stats for MVPA...
% uSimil => in-uSimil on which contrast is run {values}

% Rename settings to keep easier track...
EP = aap.tasklist.currenttask.settings;

if ~strcmp(EP.statsType, 'all-ttest') ...
        && ~strcmp(EP.statsType, 'all-signrank') ...
        && ~strcmp(EP.statsType, 'con-ttest') ...
        && ~strcmp(EP.statsType, 'con-signrank') ...
        && ~strcmp(EP.statsType, 'rmGLM') ...
        && ~strcmp(EP.statsType, 'fullGLM')
    uSimil = nanmean(uSimil,3);
end

Stats = zeros(length(EP.contrasts), length(EP.tests));

%% SHAPE EP.contrasts
if ~strcmp(EP.statsType, 'rmGLM') ...
        && ~strcmp(EP.statsType, 'fullGLM') ...
        && strcmp(EP.statsType, 'all-ttest') ...
        && strcmp(EP.statsType, 'all-signrank')
    con = NaN*zeros(size(uSimil,1), length(EP.contrasts));
else
    con = NaN*zeros(length(uSimil(:)), length(EP.contrasts));
end

for t = 1:length(EP.contrasts)
    % If it's a similarity measure, then use positive contrast; if
    % dissimilarity measure, use inverse contrast...
    if strcmp(EP.corrType, 'Pearson') ...
            || strcmp(EP.corrType, 'Spearman')
        temp = EP.contrasts(t).matrix(:);
    else
        temp = -EP.contrasts(t).matrix(:);
    end
    
    if strcmp(EP.statsType, 'rmGLM') ...
            || strcmp(EP.statsType, 'fullGLM')
        con(:,t) = reshape(repmat(mvpaa_balanceCont(temp, EP.balanceCons), [1 size(uSimil,3)]), [length(uSimil(:)) 1]);
    elseif strcmp(EP.statsType, 'all-ttest') ...
            || strcmp(EP.statsType, 'all-signrank')
        con(:,t) = reshape(repmat(temp, [1 size(uSimil,3)]), [length(uSimil(:)) 1]);
    else
        con(:,t) = mvpaa_balanceCont(temp, EP.balanceCons);
    end
end

%% DO STATS...

if strcmp(EP.statsType, 'ttest') ...
        || strcmp(EP.statsType, 'signrank')
    for t = 1:length(EP.contrasts)
        % Temporary contrast
        Tcon = con(:,t);
        
        if length(unique(Tcon(~isnan(Tcon)))) == 1
            % If normal bulk effect...
            pos = uSimil(Tcon>0);
            if ~isempty(findstr(EP.statsType, 'ttest'))
                Stats(t,1) = mean(pos(:));
                [junk, p , junk, stats] = ttest(pos(:));
                Stats(t,2) = stats.tstat;
                Stats(t,3) = p;
                Stats(t,4) = stats.sd / sqrt(stats.df);
                % JB test on all values of correlations? (or absolute values?)
                [junk, Stats(t,5)] = jbtest(pos(:)); %jbtest([pos(:) abs(neg(:))])
            else
                Stats(t,1) = median(pos(:));
                % Signed Rank for one samples
                p = signrank(pos(:));
                Stats(t,2) = tinv(1-p, length(pos(:))-1);
                Stats(t,3) = p;
            end
        else
            pos = uSimil(Tcon>0);
            neg = uSimil(Tcon<0);
            if ~isempty(findstr(EP.statsType, 'ttest'))
                Stats(t,1) = mean(pos(:)) - mean(neg(:));
                [junk, p , junk, stats] = ttest2(pos(:), neg(:));
                Stats(t,2) = stats.tstat;
                Stats(t,3) = p;
                Stats(t,4) = stats.sd / sqrt(stats.df);
                % JB test on all values of correlations? (or absolute values?)
                [junk, Stats(t,5)] = jbtest([pos(:); neg(:)]); %jbtest([pos(:) abs(neg(:))])
            else
                Stats(t,1) = median(pos(:)) - median(neg(:));
                % Signed Rank is for paired samples, so we use Rank Sum
                p = ranksum(pos(:), neg(:));
                Stats(t,2) = tinv(1-p, length(pos(:))+length(neg(:))-2);
                Stats(t,3) = p;
            end
        end
    end
elseif strcmp(EP.statsType, 'GLM') ...
        || strcmp(EP.statsType, 'rmGLM')
    for t = 1:length(EP.contrasts)
        % Temporary contrast
        Tcon = con(:,t);
        
        [junk,junk,stats] = glmfit( Tcon, uSimil(:));
        if length(unique(Tcon(~isnan(Tcon)))) == 1
            % If bulk effect
            Stats(t,1) = stats.beta(1);
            Stats(t,2) = stats.t(1);
            Stats(t,3) = stats.p(1);
            Stats(t,4) = stats.se(1);
        else
            % Otherwise...
            Stats(t,1) = stats.beta(2);
            Stats(t,2) = stats.t(2);
            Stats(t,3) = stats.p(2);
            Stats(t,4) = stats.se(2);
        end
    end
elseif strcmp(EP.statsType, 'fullGLM')
    % rmGLM considering all EP.contrasts simultaneously...
    % Careful, if EP.contrasts are not independent, you'll lose power
    % appears to give similar betas, but lower p values than rmGLM
    
    % Zero those things that are NaNs...
    con(isnan(con)) = 0;
    
    [junk,junk,stats] = glmfit( con, uSimil(:)');
    try
        % If normal bulk effect...
        Stats(:,1) = stats.beta(:);
        Stats(:,2) = stats.t(:);
        Stats(:,3) = stats.p(:);
        Stats(:,4) = stats.se(:);
    catch
        % If no bulk effect included...
        Stats(:,1) = stats.beta(2:end);
        Stats(:,2) = stats.t(2:end);
        Stats(:,3) = stats.p(2:end);
        Stats(:,4) = stats.se(2:end);
    end    
elseif strcmp(EP.statsType, 'con-ttest') ...
        || strcmp(EP.statsType, 'con-signrank')
    % Extract the value of EP.contrasts in each block*subblock comparison
    for t = 1:length(EP.contrasts)
        temp = zeros(size(uSimil,3),1);
        
        for i = 1:size(temp,2)
            % Temporary contrast
            Tcon = con(:,t);
            
            if any(Tcon(:) < 0)
                pos = uSimil(Tcon>0,i);
                neg = uSimil(Tcon<0,i);
                temp(i) = mean(pos(:)) - mean(neg(:));
            else
                pos = uSimil(Tcon>0,i);
                temp(i) = mean(pos(:));
            end
        end
        
        if ~isempty(findstr(EP.statsType, 'ttest'))
            Stats(t,1) = mean(temp);
            [junk, p , junk, stats] = ttest(temp);
            Stats(t,2) = stats.tstat;
            Stats(t,3) = p;
            Stats(t,4) = stats.sd / sqrt(stats.df);
            % JB test on all values of correlations? (or absolute values?)
            [junk, Stats(t,5)] = jbtest(temp); %jbtest([pos(:) abs(neg(:))])
        else
            Stats(t,1) = median(temp);
            % Signed Rank for one samples
            p = signrank(temp);
            Stats(t,2) = tinv(1-p, length(temp)-1);
            Stats(t,3) = p;
        end
    end
elseif strcmp(EP.statsType, 'all-ttest') ...
        || strcmp(EP.statsType, 'all-signrank')
    % Extract the value of EP.contrasts in each block*subblock comparison
    for t = 1:length(EP.contrasts)
        % Temporary contrast
        Tcon = con(:,t);
        
        temp = Tcon.*uSimil(:);
        temp = temp(~isnan(temp));
        
        if ~isempty(findstr(EP.statsType, 'ttest'))
            Stats(t,1) = mean(temp);
            [junk, p , junk, stats] = ttest(temp);
            Stats(t,2) = stats.tstat;
            Stats(t,3) = p;
            Stats(t,4) = stats.sd / sqrt(stats.df);
            % JB test on all values of correlations? (or absolute values?)
            [junk, Stats(t,5)] = jbtest(temp); %jbtest([pos(:) abs(neg(:))])
        else
            Stats(t,1) = median(temp);
            % Signed Rank for one samples
            p = signrank(temp);
            Stats(t,2) = tinv(1-p, length(temp)-1);
            Stats(t,3) = p;
        end
    end
else
    error('No EP.statsType chosen!')
end