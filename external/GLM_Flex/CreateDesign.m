function F = CreateDesign(IN)
%%%  This script creates the GLM design and associated information for
%%%  running and evaluating the GLM.
%%%
%%%  See Examples.m for some examples.
%%%  Also go to http://nmr.mgh.harvard.edu/harvardagingbrain/People/AaronSchultz/Aarons_Scripts.html
%%%  for more information.
%%%
%%%  Inputs:
%%%
%%%  IN.N_subs = Number of Subjects Per Group.
%%%
%%%  IN.Between = Vector coding number of levels per between 
%%%  subjects factor.
%%%
%%%  IN.BetweenLabs = Cell Array containing the names of each
%%%  level within each factor. Note this means a double cell Array (e.g.
%%%  {{'F1 Lev1' 'F1 Lev2'} {'F2 Lev1' 'F2 Lev2'}}.
%%%  
%%%  IN.Within = Vector coding number of levels per within 
%%%  subjects factor.
%%%
%%%  IN.WithinLabs = Cell Array containing the names of each 
%%%  level within each factor (same specification style as BetweenLabs).
%%%
%%%  IN.Interactions = Cell Array of vector coding which factor
%%%  interactions should be produced (between subjects factors are number
%%%  before within subjects factors).
%%%
%%%  IN.Covar = CellArray with one covariate per cell.
%%%  
%%%  IN.CovarLabs = Cell Array of labels for the Covariates. Note that this
%%%  is a single cell araay not a double array like with within and between
%%%  labels
%%%  
%%%  IN.CovarInt = Cell Array of Vectors coding which factors/interactions
%%%  should interact with the covariates (you can only specify one
%%%  interaction for each covariate).
%%%
%%%  IN.FactorLabs:  A single cell array with the name for each factor.
%%%  (This is opposed to providing names for each level of each factor.)
%%%
%%%  IN.EqualVar:  a row vector the length of the number of factors that
%%%  codes whether or not Equal variance is to be assumed (1), or not
%%%  assumed (0).  Not assuming means creating a whitening matrix to
%%%  correct for unequal varaince.
%%%
%%%  IN.Independent:  a row vector the length of the number of factors that
%%%  codes whether or not Independence between observations is to be 
%%%  assumed (1), or not assumed (0). Not assuming means creating a 
%%%  whitening matrix to correct for dependence.
%%%
%%%  Notes:  
%%%  1. Only Specify fields that are necessary for your design. For example
%%%  if you do not have a within subjects factor or covariates, you can
%%%  simply exclude those fields from the IN input structure.
%%% 
%%%
%%% Written by Aaron P. Schultz - aschultz@martinos.org
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
%%%  IN.N_subs = [];
%%%  IN.Between = [];
%%%  IN.BetweenLabs = {};
%%%  IN.Within = [];
%%%  IN.WithinLabs = {};
%%%  IN.Covar = {};
%%%  IN.CovarLabs = {};
%%%  IN.CovarInt = {};
%%%  IN.Interactions = [];
%%%  IN.FactorLabs = {};
%%%  IN.EqualVar = [1];
%%%  IN.Independent = [1];


if ~isfield(IN,'Between') && ~isfield(IN,'Within');
   if isfield(IN,'Covar');
      IN.N_subs = size(IN.Covar{1},1);
   end
end
N_subs = IN.N_subs;
[Between BetweenLabs Within WithinLabs Covar CovarLabs Interactions CovarInt] = FactorSetup(IN);
IN.Between = Between;
IN.BetweenLabs = BetweenLabs;
IN.Within = Within;
IN.WithinLabs = WithinLabs;
IN.Covar = Covar;
IN.CovarLabs = CovarLabs;
IN.CovarInt = CovarInt;
IN.Interactions = Interactions;

if ~isfield(IN,'FactorLabs')
    nf = numel(IN.Between)+numel(IN.Within);
    IN.FactorLabs = cellstr([repmat('Factor ',nf,1) char(num2str((1:nf)'))]);
end

isWith = [zeros(1,length(Between)) ones(1,length(Within)) zeros(1,length(Covar))];
isCovar = [zeros(1,length(Between)) zeros(1,length(Within)) ones(1,length(Covar))];
isBet = [ones(1,length(Between)) zeros(1,length(Within)) zeros(1,length(Covar))];

[IN N_subs] = Check(Within,Between,N_subs,IN);

M = CreateFactorMatrix(Between,Within,N_subs);

levs = [Between(end:-1:1) Within(end:-1:1) ones(1,length(Covar)) max(M(:,1))];

n = numel(levs)+numel(Covar);
n = n^n;
if n<1000;
    n=10000;
end
set(0,'RecursionLimit',n);

if any(isWith); intadj = 1; else intadj=0; end
Factors = SetupMainEffects(M,Covar,isWith,IN);
if numel(find(isCovar==0)) > 1-intadj
Factors = TwoWay(Factors,isWith,[isCovar 0]);
end
if numel(find(isCovar==0)) > 2-intadj
Factors = ThreeWay(Factors,isWith,[isCovar 0]);
end
if numel(find(isCovar==0)) > 3-intadj
Factors = FourWay(Factors,isWith,[isCovar 0]);
end
if numel(find(isCovar==0)) > 4-intadj
Factors = FiveWay(Factors,isWith,[isCovar 0]);
end
if numel(find(isCovar==0)) > 5-intadj
Factors = SixWay(Factors,isWith,[isCovar 0]);
end

%%% For CovarInts, make sure there is only one per covariate.
[Factors Interactions] = SetupCovarEffects(Factors,Covar,CovarInt,isWith,isCovar,Interactions);
IN.Interactions = Interactions;
[err errLevs adj adjLevs eff effLevs] = ParseDesign(IN,Factors,isWith,isCovar);

[XX assoc CovarCols name] = CreateDesignMatrix(Factors,err,M,Within,Interactions,isCovar);

% [Vi VarParts CovParts AllVar AllCov AllCovVar]  = VarCovarPartition(Factors,IN,isBet,isWith,isCovar);

try
    Vi2 = CVpart(Factors,IN,isBet,isWith,isCovar,M);
catch
    Vi2 = {[]};
end
Vi = VarCovarPartition2(IN,M);

F.FF = Factors;
F.FM = M;
F.XX = XX;
F.isBet = isBet;
F.isWith = isWith;
F.isCov = isCovar;
F.name = name;
F.CovarCols = CovarCols;
F.isCovar = isCovar;
F.Vi = Vi;
F.Vi2 = Vi2;
% F.VarParts = VarParts;
% F.CovParts = CovParts;
% F.AllVar = AllVar;
% F.AllCov = AllCov;
% F.AllCovVar = AllCovVar;
F.IN = IN;
F.err = err;
F.adj = adj;
F.eff = eff;
F.errLevs = errLevs;
F.adjLevs = adjLevs;
F.effLevs = effLevs;
F.cols = assoc;
% figure; imagesc((XX./repmat(max(XX),size(XX,1),1))); colormap(gray); shg
disp('Finished Other Stuff.');

set(0,'RecursionLimit',500)

F.IN.Between = F.IN.Between(end:-1:1);
F.IN.Within = F.IN.Within(end:-1:1);

end


function [Between BetweenLabs Within WithinLabs Covar CovarLabs Interactions CovarInt] = FactorSetup(IN)
%%% Between Factor Setup
if isfield(IN,'Between')
    Between = IN.Between(end:-1:1);
else
    Between = [];
end
if isfield(IN,'BetweenLabs')
    BetweenLabs = IN.BetweenLabs;
else
    BetweenLabs = [];
end

%%% Within Factor Setup
if isfield(IN,'Within')
    Within = IN.Within(end:-1:1);
else
    Within = [];
end
if isfield(IN,'WithinLabs')
    WithinLabs = IN.WithinLabs;
else
    WithinLabs = [];
end

%%% Covariate Setup
if isfield(IN,'Covar')
    Covar = IN.Covar;
else
    Covar = [];
end
if isfield(IN,'CovarLabs')
    CovarLabs = IN.CovarLabs;
else
    CovarLabs = [];
end

%%% Interactions
if isfield(IN,'Interactions'); 
    Interactions = IN.Interactions; 
else
    Interactions = cell(0);
end

if isfield(IN,'CovarInt'); 
    CovarInt = IN.CovarInt;
    if(numel(CovarInt)<numel(Covar))
        for ii = (length(CovarInt)+1):length(Covar)
            CovarInt{ii} = [];
        end
    end
else
    CovarInt = cell(0);
    for ii = 1:length(Covar)
    	CovarInt{ii} = [];
    end
end


for ii = 1:length(Covar);
    if numel(Covar{ii}) < (sum(IN.N_subs)*sum(Within))
        if isempty(Within)
            error('Covariate is not the right length');
        else
                nn = prod(IN.Within);
                if (nn*numel(Covar{ii})) == (sum(IN.N_subs)*sum(Within));
                    tmp = Covar{ii};
                    int = {[(1:numel(Within))+numel(Between)]};
%                     CovarInt{ii} = {[(1:numel(Within))+numel(Between)]};
                    nc = [];
                    for ll = 1:nn
                        nc(ll:nn:numel(tmp)*nn) = tmp;
                    end
                    Covar{ii} = nc';
                    
                    if isempty(CovarInt{ii}{1})
                        CovarInt{ii} = int;
                    else
                        for gg = 1:length(CovarInt{ii});
                           tmp = CovarInt{ii}{gg};
                           ch = isIn(tmp,int{1});
                           disp(tmp);
                           disp(int{1});
                           if sum(ch~=0)~=numel(int{1});
                               tmp = unique([tmp int{1}]);
                               CovarInt{ii}{gg} = tmp;
%                               error('The Covariate Interaction is Invalid');
                           end
                        end
                    end
                else
                    error('Covariate is not the right length');
                end
        end
    end
end

end

function [IN N_subs] = Check(Within,Between,N_subs,IN)

if ~isfield(IN,'FactorLabs') || isempty(IN.FactorLabs);
    IN.FactorLabs = [];
    c = 0;
    for ii = 1:length(IN.Between);
        c = c+1;
        IN.FactorLabs{c} = ['Bet' num2str(ii)];
    end
    for ii = 1:length(IN.Within);
        c = c+1;
        IN.FactorLabs{c} = ['With' num2str(ii)];
    end
    for ii = 1:length(IN.Covar);
        c = c+1;
        IN.FactorLabs{c} = ['Cov' num2str(ii)];
    end
end

if ~isfield(IN,'BetweenLabs') || isempty(IN.BetweenLabs);
    IN.BetweenLabs = [];
    c = 0;
    for ii = 1:length(IN.Between);
        c = c+1;
        for jj = 1:(IN.Between(ii));
            IN.BetweenLabs{ii}{jj} = ['F' num2str(c) num2str(jj)];
        end
    end
end

if ~isfield(IN,'WithinLabs') || isempty(IN.WithinLabs);
    IN.WithinLabs = [];
    c = 0;
    for ii = 1:length(IN.Within);
        c = c+1;
        for jj = 1:(IN.Within(ii));
            IN.WithinLabs{ii}{jj} = ['F' num2str(c) num2str(jj)];
        end
    end
end

if ~isfield(IN,'CovarLabs') || isempty(IN.CovarLabs);
    IN.CovarLabs = [];
    c = 0;
    for ii = 1:length(IN.Covar);
        c = c+1;
        IN.CovarLabs{ii}{1} = ['C' num2str(c)];
    end
end

if ~isfield(IN,'EqualVar')
    IN.EqualVar = ones(1,numel(IN.Between)+numel(IN.Within));
end
if ~isfield(IN,'Independent');
    IN.Independent = ones(1,numel(IN.Between)+numel(IN.Within));
end

if isempty(Between)
    N_subs = sum(N_subs);
end

if numel(N_subs)~=prod(Between)
    error('N_subs was not correctly specified.  For between factors there must be as many entires in N_subs as there are conditions.');
end

if numel(IN.Independent)~=(numel(IN.Between)+numel(IN.Within));
    error('Independence Flags do not match the number of Factors');
end

if numel(IN.EqualVar)~=(numel(IN.Between)+numel(IN.Within));
    error('EqualVar Flags do not match the number of Factors');
end

end

function M = CreateFactorMatrix(Between,Within,N_subs)
M = [];
c = 0;
sn = 0;
for ii = 1:numel(N_subs)
    for jj = 1:N_subs(ii);
        sn = sn+1;
        [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20] = ind2sub(Between,ii);
        a = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20];
        a = a(length(Between):-1:1);
%         a = a(1:length(Between));
        
        if prod(Within)>1
            for kk = 1:prod(Within)
                [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20] = ind2sub(Within,kk);
                b = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20];
                b = b(length(Within):-1:1);
                c = c+1;
                m = [sn a b];
                M = [M; m];
            end
        else
            c = c+1;
            m = [sn a];
            M = [M; m];
        end
    end
end
disp('Finished Creating Factor Matrix');
end

function Factors = SetupMainEffects(M,Covar,isWith,IN)
clear Factors
%Factors = cell(repmat(size(M,2)+length(Covar),1,size(M,2)+length(Covar)));
% Factors = cell(repmat(size(M,2),1,size(M,2)));
Factors = cell(size(M,2)+length(Covar),size(M,2)+length(Covar));

vec = NaN*zeros(size(M,1),1);
c = 0;
for ii = 2:size(M,2);
    XX = [];
    list = unique(M(:,ii));
    
    for jj = 1:length(list)
        ind = find(M(:,ii)==list(jj));
        vec2 = vec;
        vec2(ind) = 1;
        XX = [XX vec2];
    end       
    c = c+1;
    Factors{c,c} = XX;
end

for ii = 1:length(Covar);
    c = c+1;
    Factors{c,c} = Covar{ii};
end

XX = [];
list = unique(M(:,1));
for jj = 1:length(list)
    ind = find(M(:,1)==list(jj));
    vec2 = vec;
    vec2(ind) = 1;
    XX = [XX vec2];
end
Factors{c+1,c+1} = XX;

disp('Finished Creating Main Effects');
end

function [Factors Interactions] = SetupCovarEffects(Factors,Covar,CovarInt,isWith,isCovar,Interactions)
    CovarInds = find(isCovar==1);
    n_obs = size(Factors{1,1},1);

    for ii = 1:length(Covar)
        if isempty(CovarInt{ii})
            if ~numel(Covar{ii})==n_obs
                error('The covariate is not correct.  The covariate does not have the same number of observations as the DV');
            end
        else
            if iscell(CovarInt{ii})
                error('You may only specify one interaction per covariate.');
                return;
            end
            if numel(Covar{ii})==n_obs
                if numel(CovarInt{ii}) == 1
                    ci = [CovarInt{ii} CovarInt{ii}];
                else
                    ci = CovarInt{ii};
                end
                ind = sub2ind_aps(size(Factors),sort(ci));
                CC = [];
                tmp = Factors{ind};
                for kk = 1:size(tmp,2);
                    CC(:,kk) = Covar{ii}.*tmp(:,kk);
                end
                
                try
                    ind = sub2ind_aps(size(Factors),sort([CovarInds(ii) CovarInt{ii}]));
                catch
                    ind = sub2ind_aps([size(Factors) size(Factors,1)],sort([CovarInds(ii) CovarInt{ii}]));
                    tmp = cell([size(Factors) size(Factors,1)]);
                    tmp(1:numel(Factors)) = Factors;
                    Factors = tmp;
                end
                Factors{ind} = CC;
                Interactions{end+1} = sort([CovarInds(ii) CovarInt{ii}]);
            end
        end
    end
    
end

function [err errLevs adj adjLevs eff effLevs CovarTrack] = ParseDesign(IN,Factors,isWith,isCovar)
ss = size(Factors);

if ~isfield(IN,'Interactions');
    IN.Interactions = [];
end

for zvb = 1;
    if isfield(IN,'Within') && ~isempty(IN.Within)
        S = size(Factors,1);
        int = [];
        err = [];
        
        c = 1;
        int{1} = S;
        err(1,1) = sub2ind_aps(ss,[S S]);
        errLevs{1,1} = size(Factors{S,S},2);
        for ii = 1:length(IN.Within)
            list = sortrows(combnk((1:length(IN.Within))+numel(IN.Between),ii));
            list(:,end+1) = S;
            for jj = 1:size(list,1);
                c = c+1;
                %try
                ind = sub2ind_aps(ss,sort(list(jj,:)));
                %catch; keyboard; end 
                int{ii+1}(jj,1:size(list,2)) = sort(list(jj,:));
                err(c,1) = ind;
                
                errLevs{c,1} = [];
                tmp = sort(list(jj,:));
                tt = [];
                for kk = unique(list(jj,:))
                    tt = [tt size(Factors{(kk),(kk)},2)];
                end
                errLevs{c,1} = tt;
            end
        end
                
        adj = [];
        bet = [];
        adjLevs = [];

        for ii = 1:length(IN.Between)
            list = sortrows(combnk((1:length(IN.Between)),ii));
            bet{ii} = list;
        end
        
        c1 = 0;
        for kk = 1:length(int)
            for ll = 1:size(int{kk},1)
                c1 = c1+1;
                if err(c1)==-1
                    continue
                end
                tmp = int{kk}(ll,1:end-1);
                c2 = 0;
                for ii =1:length(bet)
                    for jj = 1:size(bet{ii},1)
                        tmp2 = sort([tmp bet{ii}(jj,:)]);

                        if numel(tmp2)==1;
                            tmp2 = [tmp2 tmp2];
                            c2 = c2+1;
                            ind = sub2ind_aps(ss,tmp2);
                            
                            if size(Factors{ind},2)>1
                                adj(c1,c2) = ind;
                                adjLevs{c1,c2} = size(Factors{tmp2(1),tmp2(1)},2);
                            end
                        else
                            
                            ints = IN.Interactions;
                            gg = [];
                            for vv = 1:length(ints)
                                gg(vv) = 1;
                                for zz = 1:length(tmp2)
                                    gg(vv) = gg(vv)*~isempty(find(ints{vv}==tmp2(zz)));
                                end
                            end
                            
                            if numel(find(gg == 1))>0
                                c2 = c2+1;
                                ind = sub2ind_aps(ss,tmp2);
                                
                                if size(Factors{ind},2)>1
                                    adj(c1,c2) = ind;
                                    
                                    tmp3 = [];
                                    for zz = unique(tmp2);
                                        tmp3 = [tmp3 size(Factors{(zz),(zz)},2)];
                                    end
                                    adjLevs{c1,c2} = tmp3;
                                end
                            else
                                c2 = c2+1;
                                adj(c1,c2) = 0;
                                adjLevs{c1,c2} = 0;
                            end
                        end
                    end
                end
            end
        end
        
        %%% Next Compute the effects the correspond to each error term.
        eff = [];
        effLevs = [];
        a = length(IN.Between);
        b = length(IN.Within);
        c = length(IN.Covar);
        Fs = [1:a (a+1):(a+b) (a+b+1):(a+b+c)];

        ints = [];
        for ii = 1:length(Fs);
            ints{ii} = Fs(ii);
        end
        effects = [ints IN.Interactions];
        isBet = ((isWith+isCovar)==0);

        c3 = 0;
        for ii = 1:length(Fs)

            list = sortrows(combnk(Fs,ii));
            for jj = 1:size(list,1)

                tmp = list(jj,:);
                %%% check if the effect is valid for the specified design
                gg = isIn(tmp,effects);    
                if gg==0; continue; end
                
                isB = isIn(tmp,find(isBet));
                isW = isIn(tmp,find(isWith));
                isC = isIn(tmp,find(isCovar));
                                      
                if mean(isB>0)==1;
                    c1 = 1;
                end
                
                if mean(isW)>0 %&& mean(isC)==0
                    n = [];
                    %tmp2 = tmp(~isBet(tmp));
                    tmp2 = tmp(~isBet(tmp) & ~isCovar(tmp));
                    
                    for hh = 1:length(err)
                       med = ind2sub_aps(size(Factors),err(hh)); 
                       med = [med(1) med(find(med(2:end)~=1)+1)];
                       med = med(find(med~=S));
                       if isempty(med);
                           continue;
                       end
                       n(hh) = mean(isIn(tmp2,med)~=0);
                    end

                    ind = find(n==1);
                    c1 = ind(1);
                end
                
                if sum(isB)==0 && sum(isW)==0 && sum(isC)>0
                    in = tmp-sum(isBet+isWith);
                    if isempty(IN.CovarInt{in}{1});
                        c1 = find(err==-1);
                    else
                        continue;
                    end
                end
                
                if numel(tmp)==1;
                    tmp = [tmp tmp];
                end
                
                ind = sub2ind_aps(ss,sort(tmp));
                try
                    eee = eff(c1,:);
                    tt = find(eee~=0);
                    c2 = max(tt)+1;
                catch
                    c2 = 1;
                end
                
                if isempty(c2); c2 = 1; end;
                
                eff(c1, c2) = ind;
                
                
                tmp3 = [];
                for zz = unique(tmp);
                    tmp3 = [tmp3 size(Factors{(zz),(zz)},2)];
                end
                effLevs{c1,c2} = tmp3;
            end
        end
    end
end
for zvb = 1
    if ~isfield(IN,'Within')  || isempty(IN.Within)
        err = -1;
        errLevs = [];
        adj = -1;
        adjLevs = [];
        eff = [];
        c1 = 0;
        
        
        Fs = [1:length(IN.Between) (length(IN.Between)+1):(length(IN.Between)+length(IN.Covar))];
        ints = [];
        for ii = 1:length(Fs);
            ints{ii} = Fs(ii);
        end
        effects = [ints IN.Interactions]; 
        
        for ii = 1:length(Fs)
            list = sortrows(combnk(Fs,ii));
            for jj = 1:size(list,1)
                
                tmp = list(jj,:);
                %%% check if the effect is valid for the specified design
                gg = isIn(tmp,effects);    
                if gg==0; continue; end
                
                c1 = c1+1;
                
                if numel(tmp)==1;
                    tmp = [tmp tmp];
                end
                
                ind = sub2ind_aps(ss,sort(tmp));
                eff(c1) = ind;
                
                tmp3 = [];
                for zz = unique(tmp);
                    tmp3 = [tmp3 size(Factors{(zz),(zz)},2)];
                end
                effLevs{c1} = tmp3;
            end
        end
    end
end

for ii = 1:length(err);
    if err(ii)~=-1
        try
            tmp = Factors{err(ii)};
            if numel(find(nansum(tmp)==1))==size(tmp,2)
                err(ii) = -1;
                adj(ii,:) = 0;
            end
        catch
            if ii == numel(err);
                err(ii) = -1;
                adj(ii) = 0;
            end
        end
    end
end


end

function [XX assoc CovarCols name] = CreateDesignMatrix(Factors,err,M,Within,Interactions,isCovar)
% There will be a problem with assoc when the covariate is equal to
% subjects and not obs since there will be no main effect for the covariate.
cc = 0;
XX = [];
assoc = [];
CovarCols = [];
name = [];

bb = 0;
for ii = 1:size(Factors,1)-1;
    if isCovar(ii)==0
        cc = cc+1;
        XX = [XX Factors{ii,ii}];
        assoc = [assoc ii*ones(1,size(Factors{ii,ii},2))];
        CovarCols = [CovarCols isCovar(cc).*ones(1,size(Factors{ii,ii},2))];
        for jj = 1:size(Factors{ii,ii},2)
            name{end+1} = ['F_' num2str(ii) '_' num2str(jj)];
        end
    else
        if sum([Interactions{:}] == ii)==0
            cc = cc+1;
            XX = [XX Factors{ii,ii}];
            assoc = [assoc ii*ones(1,size(Factors{ii,ii},2))];
            CovarCols = [CovarCols isCovar(cc).*ones(1,size(Factors{ii,ii},2))];
            for jj = 1:size(Factors{ii,ii},2)
                name{end+1} = ['F_' num2str(ii) '_' num2str(jj)];
            end
        end
    end
end

for ii = 1:length(Interactions);
    tmp = sort(Interactions{ii});
    cc = cc+1;
    if numel(tmp)==2
        XX = [XX Factors{tmp(1),tmp(2)}];
        ch = sum(isIn(tmp,find(isCovar==1)))>0;
        CovarCols = [CovarCols ch.*ones(1,size(Factors{tmp(1),tmp(2)},2))];
        for jj = 1:size(Factors{tmp(1),tmp(1)},2)
            for kk = 1:size(Factors{tmp(2),tmp(2)},2)
                name{end+1} = ['F_' num2str(tmp(1)) '_' num2str(jj) 'xF_' num2str(tmp(2)) '_' num2str(kk)];
            end
        end
    elseif numel(tmp)==3
        XX = [XX Factors{tmp(1),tmp(2), tmp(3)}];
        ch = sum(isIn(tmp,find(isCovar==1)))>0;
        CovarCols = [CovarCols ch.*ones(1,size(Factors{tmp(1),tmp(2), tmp(3)},2))];
        for jj = 1:size(Factors{tmp(1),tmp(1)},2)
            for kk = 1:size(Factors{tmp(2),tmp(2)},2)
                for ll = 1:size(Factors{tmp(3),tmp(3)},2)
                    name{end+1} = ['F_' num2str(tmp(1)) '_' num2str(jj) 'xF_' num2str(tmp(2)) '_' num2str(kk) 'xF_' num2str(tmp(3)) '_' num2str(ll)];
                end
            end
        end
    elseif numel(tmp)==4
        XX = [XX Factors{tmp(1),tmp(2), tmp(3), tmp(4)}];
        ch = sum(isIn(tmp,find(isCovar==1)))>0;
        CovarCols = [CovarCols ch.*ones(1,size(Factors{tmp(1),tmp(2), tmp(3), tmp(4)},2))];
        for jj = 1:size(Factors{tmp(1),tmp(1)},2)
            for kk = 1:size(Factors{tmp(2),tmp(2)},2)
                for ll = 1:size(Factors{tmp(3),tmp(3)},2)
                    for mm = 1:size(Factors{tmp(4),tmp(4)},2)
                        name{end+1} = ['F_' num2str(tmp(1)) '_' num2str(jj) 'xF_' num2str(tmp(2)) '_' num2str(kk) 'xF_' num2str(tmp(3)) '_' num2str(ll) 'xF_' num2str(tmp(4)) '_' num2str(mm)];
                    end
                end
            end
        end
    elseif numel(tmp)==5
        XX = [XX Factors{tmp(1),tmp(2), tmp(3), tmp(4),tmp(5)}];
        ch = sum(isIn(tmp,find(isCovar==1)))>0;
        CovarCols = [CovarCols ch.*ones(1,size(Factors{tmp(1),tmp(2), tmp(3), tmp(4),tmp(5)},2))];
        for jj = 1:size(Factors{tmp(1),tmp(1)},2)
            for kk = 1:size(Factors{tmp(2),tmp(2)},2)
                for ll = 1:size(Factors{tmp(3),tmp(3)},2)
                    for mm = 1:size(Factors{tmp(4),tmp(4)},2)
                        for nn = 1:size(Factors{tmp(5),tmp(5)},2)
                            name{end+1} = ['F_' num2str(tmp(1)) '_' num2str(jj) 'xF_' num2str(tmp(2)) '_' num2str(kk) 'xF_' num2str(tmp(3)) '_' num2str(ll) 'xF_' num2str(tmp(4)) '_' num2str(mm) 'xF_' num2str(tmp(5)) '_' num2str(nn)];
                        end
                    end
                end
            end
        end
    elseif numel(tmp)==6
        XX = [XX Factors{tmp(1),tmp(2), tmp(3), tmp(4),tmp(5),tmp(6)}];
        ch = sum(isIn(tmp,find(isCovar==1)))>0;
        CovarCols = [CovarCols ch.*ones(1,size(Factors{tmp(1),tmp(2), tmp(3), tmp(4),tmp(5),tmp(6)},2))];
        for jj = 1:size(Factors{tmp(1),tmp(1)},2)
            for kk = 1:size(Factors{tmp(2),tmp(2)},2)
                for ll = 1:size(Factors{tmp(3),tmp(3)},2)
                    for mm = 1:size(Factors{tmp(4),tmp(4)},2)
                        for nn = 1:size(Factors{tmp(5),tmp(5)},2)
                            for oo = 1:size(Factors{tmp(6),tmp(6)},2)
                                name{end+1} = ['F_' num2str(tmp(1)) '_' num2str(jj) 'xF_' num2str(tmp(2)) '_' num2str(kk) 'xF_' num2str(tmp(3)) '_' num2str(ll) 'xF_' num2str(tmp(4)) '_' num2str(mm) 'xF_' num2str(tmp(4)) '_' num2str(nn)  'xF_' num2str(tmp(6)) '_' num2str(oo)];
                            end
                        end
                    end
                end
            end
        end
    end
end

for ii = 1:length(err)
    if err(ii)~=-1
        Facs = ind2sub_aps(size(Factors),err(ii));
        ind = find(diff(Facs)<0);
        if ~isempty(ind);
            Facs = unique(Facs(1:ind(1)));
        else
            Facs = unique(Facs);
        end
        XX = [XX Factors{err(ii)}];
        if numel(Facs)==1 && Facs == size(Factors,1)
           for zz = 1:size(Factors{err(ii)},2)
               name{end+1} = ['S' num2str(zz)];
           end
        end
        Facs = Facs(find(Facs<numel(isCovar)));
        if isempty(Facs)
            CovarCols = [CovarCols zeros(1,size(Factors{err(ii)},2))];
        else
            if mean(isCovar(Facs))>0
                CovarCols = [CovarCols ones(1,size(Factors{err(ii)},2))];
            else
                CovarCols = [CovarCols zeros(1,size(Factors{err(ii)},2))];
            end
        end
    end
end

% if any(CovarCols)
%    XX = [ones(size(XX,1),1) XX];
%    CovarCols = [0 CovarCols];
% end

disp('Finished Creating Design Matrix.');
end

function Vi = CVpart(Factors,IN,isBet,isWith,isCovar,M)
    Vi = [];
    ss = size(Factors);
    nn = sum(IN.N_subs);
    sc = prod(IN.Within);
    n = nn*sc;
    
    %%%
    vv = find(IN.EqualVar==0);
    wh = sub2ind_aps(ss,vv);
    part = Factors{wh};
    cc = 0;
    for ii = 1:size(part,2)
       tmp = zeros(n,n);
       ind = find(~isnan(part(:,ii)));
       
       ind = sub2ind(size(tmp), ind, ind);
       tmp(ind)=1;
       cc = cc+1;
       Vi{cc} = sparse(tmp);
    end
    %%%
    vv = find(IN.Independent==0);
    wh = sub2ind_aps(ss,vv);
    if numel(wh)==1; wh = [wh wh]; end
    part = Factors{wh};
    combs = combnk(1:size(part,2),2);
    for ii = 1:size(combs,1)
        tmp = zeros(n,n);
        ind1 = find(~isnan(part(:,combs(ii,1))));
        ind2 = find(~isnan(part(:,combs(ii,2))));
        ind = sub2ind(size(tmp), [ind1; ind2], [ind2; ind1]);
        tmp(ind) = 1;
        cc = cc+1;
        Vi{cc} = sparse(tmp);
    end
end

function [Vi VarParts CovParts AllVar AllCov AllCovVar] = VarCovarPartition(Factors,IN,isBet,isWith,isCovar)
    
    VarParts = [];
    CovParts = [];
    AllCov = [];
    AllCovVar = [];
%     AllVar = [];
    Vi = [];
    nn = sum(IN.N_subs);
    sc = prod(IN.Within);
    n = nn*sc;
    ss = size(Factors);

    %%% Create Variance Partitions
    ind1 = sub2ind_aps(ss,1:(sum(isBet)+sum(isWith)));
    for ii = 1:size(Factors{ind1},2)
        tmp = zeros(n,n);
        ind2 = find(Factors{ind1}(:,ii)==1);
        ind3 = ((ind2-1)*(n+1))+1;
        tmp(ind3) = 1;
        AllVar{ii} = tmp;
    end
    %keyboard
    wh = sort(find(IN.EqualVar==0));
    if numel(wh)==1;
       wh = [wh wh]; 
    end
    
    if numel(wh)>0
        ind1 = sub2ind_aps(ss,wh);
        for ii = 1:size(Factors{ind1},2)
            tmp = zeros(n,n);
            ind2 = find(Factors{ind1}(:,ii)==1);
            ind3 = ((ind2-1)*(n+1))+1;
            tmp(ind3) = 1;
            Vi{ii} = tmp;
        end
        c = ii;
        VarParts = 1:ii;
    else
        Vi{1} = eye(n);
        c = 1;
    end
    %%% Create Covaraince Partitions
    
    if sum(isWith)>0
        wh = sort(find(isWith==1));
        if numel(wh)==1;
            wh = [wh wh];
        end
        ind1 = sub2ind_aps(ss,wh);
        for ii = 1:size(Factors{ind1},2)
            tmp = zeros(n,n);
            ind2 = find(Factors{ind1}(:,ii)==1);
            ind3 = ((ind2-1)*(n+1))+1;
            tmp(ind3) = 1;
            AllCovVar{ii} = tmp;
        end
        
        list =  combnk(1:size(Factors{ind1},2),2);
        for ii = 1:size(list,1)
            a = Factors{ind1}(:,list(ii,1));
            b = Factors{ind1}(:,list(ii,2));
            in1 = find(a==1);
            in2 = find(b==1);
            tmp = zeros(n,n);
            in3 = sub2ind(size(tmp),in1,in2);
            tmp(in3) = 1;
            in3 = sub2ind(size(tmp),in2,in1);
            tmp(in3) = 1;
            AllCov{ii} = tmp;
        end
        
        wh2 = find(isBet);
        if numel(wh2)==1;
            wh2 = [wh2 wh2];
        end
        if numel(wh) == 0;
           Check = ones(n,1);
        else
           Check = Factors{sub2ind_aps(ss,wh2)};
        end
        wh = unique([find(isBet) sort(intersect(find(IN.Independent==0), find(isWith==1)))]);
        ind1 = sub2ind_aps(ss,wh);
        list = combnk(1:size(Factors{ind1},2),2);
        ee = length(Vi)+1;
        for ii = 1:size(list,1)
            a = Factors{ind1}(:,list(ii,1));
            b = Factors{ind1}(:,list(ii,2));
            in1 = find(a==1);
            in2 = find(b==1);
            ok = 0;

            for jj = 1:size(Check,2);
               in3 = find(Check(:,jj)==1);
               g1 = sum(isIn(in1,in3)>0)+sum(isIn(in2,in3)>0);
               g2 = numel(in1) + numel(in2);
               if g1==g2;
                   %keyboard;
                   ok = 1;
                   break
               end
            end

            if ok == 0;
                continue;
                %keyboard;
            end
            
            tmp = zeros(n,n);
            in3 = sub2ind(size(tmp),in1,in2);
            tmp(in3) = 1;
            in3 = sub2ind(size(tmp),in2,in1);
            tmp(in3) = 1;
            numel(Vi)
            Vi{end+1} =tmp;
        end

        CovParts = ee:length(Vi);
    end
end

function [Vi] = VarCovarPartition2(IN,FM)
    SPM = [];
    SPM.factor = [];
    c = 0;
    for ii = 1:length(IN.Between)
        c = c+1;
        SPM.factor(c).name     = IN.FactorLabs{ii};
        SPM.factor(c).levels   = IN.Between(ii);
        % Ancova options
        SPM.factor(c).gmsca = 0;
        SPM.factor(c).ancova   = 0;
        % Nonsphericity options
        SPM.factor(c).variance = 1-IN.EqualVar(c);
        SPM.factor(c).dept     = 1-IN.Independent(c);
    end
    for ii = 1:length(IN.Within)
        c = c+1;
        SPM.factor(c).name     = IN.FactorLabs{ii};
        SPM.factor(c).levels   = IN.Within(ii);
        % Ancova options
        SPM.factor(c).gmsca = 0;
        SPM.factor(c).ancova   = 0;
        % Nonsphericity options
        SPM.factor(c).variance = 1-IN.EqualVar(c);
        SPM.factor(c).dept     = 1-IN.Independent(c);
    end
    
    if all(FM(:,1)==(1:size(FM,1))')
        FM(:,1)=1;
    end
    
    
    SPM.xVi.I = ones(size(FM,1),4);
    SPM.xVi.I(:,1:size(FM,2)) = FM;
    SPM = spm_get_vc(SPM);
    
    try
        Vi = SPM.xVi.Vi;
    catch
        Vi{1} = SPM.xVi.V;
    end
end

function out = isIn(a,b)
%%% is a in b?
out = [];
if isempty(b);
    out = 0;
    return
else
    if iscell(b)
        for kk = 1:length(b);
            ch = [];
            for ii = 1:length(a);
                ind = find(b{kk}==a(ii));
                if isempty(ind);
                    ind = 0;
                else
                    ind = 1;
                end
                ch(ii) = ind;
            end
            if mean(ch)==1;
                out = kk;
                return;
            end
        end
        if isempty(out)
            out = 0;
            return
        end
    else
        ch = [];
        for ii = 1:length(a);
            ind = find(b==a(ii));
            if isempty(ind);
                ind = 0;
            end
            ch(ii) = ind;  
        end
%         if mean(ch)==1;
%             out = 1;
%         else
%             out = 0;
%         end
        out = ch;
    end
end
end

function Factors = TwoWay(Factors,isWith,covs)
%%% Two Way Interactions
list = combnk(find(covs==0),2);
for ii = 1:size(list,1);
    %     tmp = recComb({Factors{list(ii,1),list(ii,1)} Factors{list(ii,2),list(ii,2)}});
    tmp = zeros(size(Factors{1,1},1),prod([size(Factors{list(ii,1),list(ii,1)},2) size(Factors{list(ii,2),list(ii,2)},2)])); 
    c = 0;
    for jj = 1:size(Factors{list(ii,1),list(ii,1)},2)
        for kk = 1:size(Factors{list(ii,2),list(ii,2)},2)
            c = c+1;
            tmp(:,c) = Factors{list(ii,1),list(ii,1)}(:,jj).*Factors{list(ii,2),list(ii,2)}(:,kk);
        end
    end
    if any(nansum(tmp)==0)
        continue;
    else
        if all(nansum(tmp)==1);
            Factors{list(ii,1),list(ii,2)} = [];
        else
            Factors{list(ii,1),list(ii,2)} = tmp;
        end
    end
end
disp('Finished Creating Two Way Interactions');
end

function Factors = ThreeWay(Factors,isWith,covs)
%%% Three Way Interactions
list = combnk(find(covs==0),3);
for ii = 1:size(list,1);
    %tmp = recComb({Factors{list(ii,1),list(ii,1)} Factors{list(ii,2),list(ii,2)} Factors{list(ii,3),list(ii,3)}});
    tmp = zeros(size(Factors{1,1},1),prod([size(Factors{list(ii,1),list(ii,1)},2) size(Factors{list(ii,2),list(ii,2)},2) size(Factors{list(ii,3),list(ii,3)},2)])); 
    c = 0;
    for jj = 1:size(Factors{list(ii,1),list(ii,1)},2)
        for kk = 1:size(Factors{list(ii,2),list(ii,2)},2)
            for ll = 1:size(Factors{list(ii,3),list(ii,3)},2)
                c = c+1;
                tmp(:,c) = Factors{list(ii,1),list(ii,1)}(:,jj).*Factors{list(ii,2),list(ii,2)}(:,kk).*Factors{list(ii,3),list(ii,3)}(:,ll);
            end
        end
    end
    if any(nansum(tmp)==0)
        continue;
    else
        if all(nansum(tmp)==1);
            %continue;
            Factors{list(ii,1),list(ii,2),list(ii,3)} = [];
        else
            Factors{list(ii,1),list(ii,2),list(ii,3)} = tmp;
        end
    end
end
    disp('Finished Creating Three Way Interactions');
end

function Factors = FourWay(Factors,isWith,covs)
%%% Four Way Interactions
list = combnk(find(covs==0),4);
for ii = 1:size(list,1);
    %tmp = recComb({Factors{list(ii,1),list(ii,1)} Factors{list(ii,2),list(ii,2)} Factors{list(ii,3),list(ii,3)}  Factors{list(ii,4),list(ii,4)}});
    tmp = zeros(size(Factors{1,1},1),prod([size(Factors{list(ii,1),list(ii,1)},2) size(Factors{list(ii,2),list(ii,2)},2) size(Factors{list(ii,3),list(ii,3)},2) size(Factors{list(ii,4),list(ii,4)},2)]));
    c = 0;
    for jj = 1:size(Factors{list(ii,1),list(ii,1)},2)
        for kk = 1:size(Factors{list(ii,2),list(ii,2)},2)
            for ll = 1:size(Factors{list(ii,3),list(ii,3)},2)
                for mm = 1:size(Factors{list(ii,4),list(ii,4)},2)
                    c = c+1;
                    tmp(:,c) = Factors{list(ii,1),list(ii,1)}(:,jj).*Factors{list(ii,2),list(ii,2)}(:,kk).*Factors{list(ii,3),list(ii,3)}(:,ll).*Factors{list(ii,4),list(ii,4)}(:,mm);
                end
            end
        end
    end
    
    if any(nansum(tmp)==0)
        continue;
    else
        if all(nansum(tmp)==1);
            Factors{list(ii,1),list(ii,2),list(ii,3),list(ii,4)} = [];
        else
            Factors{list(ii,1),list(ii,2),list(ii,3),list(ii,4)} = tmp;
        end
    end
end
disp('Finished Creating Four Way Interactions');
end

function Factors = FiveWay(Factors,isWith,covs)
%%% Five Way Interactions
list = combnk(find(covs==0),5);
for ii = 1:size(list,1);
    %tmp = recComb({Factors{list(ii,1),list(ii,1)} Factors{list(ii,2),list(ii,2)} Factors{list(ii,3),list(ii,3)}  Factors{list(ii,4),list(ii,4)} Factors{list(ii,5),list(ii,5)}});
    tmp = zeros(size(Factors{1,1},1),prod([size(Factors{list(ii,1),list(ii,1)},2) size(Factors{list(ii,2),list(ii,2)},2) size(Factors{list(ii,3),list(ii,3)},2) size(Factors{list(ii,4),list(ii,4)},2) size(Factors{list(ii,5),list(ii,5)},2)]));

    c = 0;
    for jj = 1:size(Factors{list(ii,1),list(ii,1)},2)
        for kk = 1:size(Factors{list(ii,2),list(ii,2)},2)
            for ll = 1:size(Factors{list(ii,3),list(ii,3)},2)
                for mm = 1:size(Factors{list(ii,4),list(ii,4)},2)
                    for nn = 1:size(Factors{list(ii,5),list(ii,5)},2)
                        c = c+1;
                        tmp(:,c) = Factors{list(ii,1),list(ii,1)}(:,jj).*Factors{list(ii,2),list(ii,2)}(:,kk).*Factors{list(ii,3),list(ii,3)}(:,ll).*Factors{list(ii,4),list(ii,4)}(:,mm).*Factors{list(ii,5),list(ii,5)}(:,nn);
                    end
                end
            end
        end
    end
    
    if any(nansum(tmp)==0)
        continue;
    else
        if all(nansum(tmp)==1);
            Factors{list(ii,1),list(ii,2),list(ii,3),list(ii,4),list(ii,5)} = [];
        else
            Factors{list(ii,1),list(ii,2),list(ii,3),list(ii,4),list(ii,5)} = tmp;
        end
    end
end
disp('Finished Creating Five Way Interactions');
end

function Factors = SixWay(Factors,isWith,covs)
%%% Six Way Interactions
list = combnk(find(covs==0),6);
for ii = 1:size(list,1);
    %tmp = recComb({Factors{list(ii,1),list(ii,1)} Factors{list(ii,2),list(ii,2)} Factors{list(ii,3),list(ii,3)}  Factors{list(ii,4),list(ii,4)} Factors{list(ii,5),list(ii,5)} Factors{list(ii,6),list(ii,6)}});
    tmp = zeros(size(Factors{1,1},1),prod([size(Factors{list(ii,1),list(ii,1)},2) size(Factors{list(ii,2),list(ii,2)},2) size(Factors{list(ii,3),list(ii,3)},2) size(Factors{list(ii,4),list(ii,4)},2) size(Factors{list(ii,5),list(ii,5)},2) size(Factors{list(ii,6),list(ii,6)},2)]));

    c = 0;
    for jj = 1:size(Factors{list(ii,1),list(ii,1)},2)
        for kk = 1:size(Factors{list(ii,2),list(ii,2)},2)
            for ll = 1:size(Factors{list(ii,3),list(ii,3)},2)
                for mm = 1:size(Factors{list(ii,4),list(ii,4)},2)
                    for nn = 1:size(Factors{list(ii,5),list(ii,5)},2)
                        for oo = 1:size(Factors{list(ii,6),list(ii,6)},2)
                            c = c+1;
                            tmp(:,c) = Factors{list(ii,1),list(ii,1)}(:,jj).*Factors{list(ii,2),list(ii,2)}(:,kk).*Factors{list(ii,3),list(ii,3)}(:,ll).*Factors{list(ii,4),list(ii,4)}(:,mm).*Factors{list(ii,5),list(ii,5)}(:,nn).*Factors{list(ii,6),list(ii,6)}(:,oo);
                        end
                    end
                end
            end
        end
    end
    
    if any(nansum(tmp)==0)
        continue;
    else
        if all(nansum(tmp)==1);
            Factors{list(ii,1),list(ii,2),list(ii,3),list(ii,4),list(ii,5),list(ii,6)} = [];
        else
            Factors{list(ii,1),list(ii,2),list(ii,3),list(ii,4),list(ii,5),list(ii,6)} = tmp;
        end
    end
end
disp('Finished Creating Six Way Interactions');
end

function out = recComb(in,lens,wh,count,out)

%     try
%         if count > 240
%             disp(count)
%             keyboard
%         end
%     end
    if nargin == 1
        for ii = 1:length(in);
            lens(ii) = size(in{ii},2);
        end
        wh = ones(size(lens));
        out = zeros(size(in{1},1),prod(lens));
        count = 1;
    end

    b = ones(size(in{1},1),1);
    for ii = 1:length(in)
        if isempty(in{ii});
            out = [];
            return;
        end
        b = b.*in{ii}(:,wh(ii));
    end
    out(:,count) = b;


    for jj = size(wh,2):-1:1
        if wh(jj)<lens(jj);
            wh(jj) = wh(jj)+1;

            for kk = jj+1:size(wh,2)
                wh(kk) = 1;
            end

            break
        end
    end


    if mean(wh==lens)==1
        count = count+1;
        
        b = ones(size(in{1},1),1);
        for ii = 1:length(in)
            b = b.*in{ii}(:,wh(ii));
        end
        out(:,count) = b;
        return
    else
        count = count+1;
        out = recComb(in,lens,wh,count,out);
    end
end