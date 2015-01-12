function F = estimateGLM(y,F,custom)
%%% You can use this script to do stand GLM analyses on individual DVs.
%%% This will display standard ANOVA tables as output.
%%%
%%% Inputs:
%%% y = a column vector of data mated to the design specified in F
%%%
%%% F = Design structure from CreateDesign.m
%%%
%%% custom = a data structure used to specify custom contrasts
%%%     custom.c  = contrast vector/matrix
%%%     custom.ET = the error term to use to evaluate the effect.
%%%
%%% Outputs:
%%% The computed statistics are put in the data structure F and returned.
%%%
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
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

spaceN = 35;

F = estimateError(y,F);

xx = F.XX;
% if isfield(F,'X')
%     x = F.xx;
if isfield(F,'wX')
    x = F.wX;
else
    x = xx;
    x(find(isnan(x)))=0;
end
b = pinv(x)*y;
rr = eye(size(x,1))-(x*pinv(x));

F.B = b;
F.ResidMat = rr;

[groupID, junk, groups] = unique(F.FM(:,2:end),'rows');
if numel(groupID)>1 && exist('vartestn.m','file')>0
    [lp VarTest] = vartestn(y,groups,'off','robust');
    VarTest.p = lp;
    F.VarTest = VarTest;
    
    if lp<.05
        warning('Variances are not equal! Check F.VarTest for more details.')
        %vartestn(y,groups,'on','robust');
    end
end
skip = 0;
% if numel(F.IN.Within)>0
%     [junk, junk, groups] = unique(F.FM(:,2:end),'rows');
%     list = unique(groups);
%     cvd = [];
%     for ii = 1:length(list);
%         cvd(:,ii) = y(find(groups == list(ii)));
%     end
%     cvm = cov(cvd);
%     Sphericity = BoxE(cvm);
%     F.BoxE = Sphericity;
%     if Sphericity < .5
%         warning('Data is violoating the sphericity assumption.');
%         disp(['Box''s  epsilon = ' num2str(Sphericity)]);
%     end
%     %%% Eventually I should get this to do the df adjustment based on the
%     %%% sphericity correction.
% end

if nargin==2
    spacer = repmat(' ',1,spaceN);
    fprintf('\n%s%s\t%s\t%s\t%s\t%s\n',spacer,'df','Sum Sq', 'Mean Sq','F-value', 'P-value');
    
    for ii = 1:size(F.eff,1)
        for jj = 1:size(F.eff,2)
            eff = F.eff(ii,jj);
            if eff == 0;
                skip = 1;
                continue;
            end
            
            Facs = ind2sub_aps(size(F.FF),eff);
            ind = find(diff(Facs)<=0);
            if ~isempty(ind)
                Facs = unique(Facs(1:ind(1)));
            else
                Facs = unique(Facs);
            end
            
            [r c cc r2] = MakeContrastMatrix('a', F.FF{eff}, F.effLevs{ii,jj},xx,x,F.CovarCols);
            
            %             if isfield(F,'CovarControl')
            %                 if F.CovarControl == 1
            %                     if any(F.isCovar(Facs))
            %                         [r c cc r2] = MakeContrastMatrix2('a',F.FF{eff}, F.effLevs{ii,jj},xx,x,F.CovarCols);
            %                     else
            %                         [r c cc r2] = MakeContrastMatrix2('b',F.FF{eff}, F.effLevs{ii,jj},xx,x,F.CovarCols);
            %                     end
            %                 else
            %                     [r c cc r2] = MakeContrastMatrix2('a',F.FF{eff}, F.effLevs{ii,jj},xx,x,F.CovarCols);
            %                 end
            %             else
            %                 [r c cc r2] = MakeContrastMatrix2('a', F.FF{eff}, F.effLevs{ii,jj},xx,x,F.CovarCols);
            %             end
            
            F.effCon{ii,jj} = c;
            SS = LoopEstimate(b,x,r-rr);
            
            df = size(c,2);
            
            
            F.EffectSS{ii,jj} = SS;
            F.EffectDF(ii,jj) = df;
            
            df1 = df;
            df2 = F.ErrorDF(ii);
            SS1 = SS;
            SS2 = F.ErrorSS{ii};
            MS1 = SS1/df1;
            MS2 = SS2/df2;
            F1 = MS1./MS2;
            p = 1-cdf('f',F1,df1,df2);
            labs = F.IN.FactorLabs(Facs);
            name = [];
            for kk = 1:length(labs);
                if kk == 1 && numel(labs)==1
                    name = [name labs{kk}];
                elseif kk < numel(labs) && numel(labs)>1
                    name = [name labs{kk} '_by_' ];
                else
                    name = [name labs{kk}];
                end
            end
            
            try
                spacer = repmat(' ',1,spaceN-numel(name));
            catch
                spacer = '  ';
            end
            
            fprintf('%s%s%0.0f\t%4.3f\t%4.3f\t%4.3f\t%1.6f\n',name,spacer,df1,SS1,MS1,F1,p);
            
            
            %             if mean(F.isCovar(Facs))>0 && mean(F.isCovar(Facs))<1
            %                 SS = LoopEstimate(b,x,r2-rr);
            %                 df = size(cc,2);
            %
            %                 %F.EffectSS{ii,jj} = SS;
            %                 %F.EffectDF(ii,jj) = df;
            %
            %                 df1 = df;
            %                 df2 = F.ErrorDF(ii);
            %                 SS1 = SS;
            %                 SS2 = F.ErrorSS{ii};
            %                 MS1 = SS1/df1;
            %                 MS2 = SS2/df2;
            %                 F1 = MS1./MS2;
            %                 p = 1-cdf('f',F1,df1,df2);
            %                 keyboard;
            %             end
            
        end
        
        if skip==0
            spacer = repmat(' ',1,spaceN-numel('Error'));
            fprintf('%s%s%0.0f\t%4.3f\t%4.3f\n','Error',spacer,df2,SS2,MS2);
        else
            skip = 0;
        end
        if ii<size(F.eff,1)
            spacer = repmat(' ',1,spaceN);
            fprintf('\n%s%s\t%s\t%s\t%s\t%s\n',spacer,'df','Sum Sq', 'Mean Sq','F-value', 'P-value');
        end
        
        
    end
elseif nargin == 3
     spacer = repmat(' ',1,spaceN);
     fprintf('\n%s%s\t%s\t%s\t%s\t%s\n',spacer,'df','Sum Sq', 'Mean Sq','F-value', 'P-value');
    for zz = 1:length(custom)
        try
            [c cc] = getC(custom(zz),F.XX,F.CovarCols);
            c0 = eye(size(x,2))-(c*pinv(c));
            x0 = x*c0;
            r = eye(size(x0,1))-(x0*pinv(x0));
        catch
            if ~isfield(custom(zz),'c') || isempty(custom(zz).c)
                [r c c2] = MakeContrastMatrix('a',custom(zz).DesignGroup, custom(zz).levs,xx,x);
                %[r c c2] = MakeContrastMatrix('a',F.XX(:,custom(zz).DesignGroup), custom(zz).levs,xx,x);
            else
                c = custom(zz).c;
                c0 = eye(size(x,2))-(c*pinv(c));
                x0 = x*c0;
                r = eye(size(x0,1))-(x0*pinv(x0));
            end
        end
        
        SS = LoopEstimate(b,x,r-rr);
        
        df1 = size(c,2);
        
        if ~isfield(F,'effCon');
            F.effCon = [];
        end
        if ~isfield(F,'EffectSS');
            F.EffectSS = [];
        end
        if ~isfield(F,'EffectDF');
            F.EffectDF = [];
        end
        
        F.effCon{custom(zz).ET,end+1} = c;
        F.EffectSS{custom(zz).ET,end+1} = SS;
        F.EffectDF(custom(zz).ET,end+1) = df1;
        
        df2 = F.ErrorDF(custom(zz).ET);
        SS1 = SS;
        SS2 = F.ErrorSS{custom(zz).ET};
        MS1 = SS1/df1;
        MS2 = SS2/df2;
        F1 = MS1./MS2;
        p = 1-cdf('f',F1,df1,df2);
        
        spacer = repmat(' ',1,spaceN);
        try
            name = custom(zz).name;
        catch
            name = 'Custom Effect';
        end
        spacer = repmat(' ',1,spaceN-numel(name));
        fprintf('%s%s%0.0f\t%4.3f\t%4.3f\t%4.3f\t%1.6f\n',name,spacer,df1,SS1,MS1,F1,p);
        if zz<numel(custom) 
            if custom(zz).ET~=custom(zz+1).ET
                spacer = repmat(' ',1,spaceN-numel('Error'));
                fprintf('%s%s%0.0f\t%4.3f\t%4.3f\n','Error',spacer,df2,SS2,MS2);
            end
        else
            spacer = repmat(' ',1,spaceN-numel('Error'));
            fprintf('%s%s%0.0f\t%4.3f\t%4.3f\n','Error',spacer,df2,SS2,MS2);
        end
    end
end
end

function [c cc] = getC(Con,x,CovarCols)
if iscell(Con.Groups)
    if iscell(Con.Groups{1})
        for zz = 1:length(Con.Groups)
            tmpX = zeros(size(x,1),1)*NaN;
            tmpX(Con.Groups{zz}{1})=1;
            cc = MakePreCons(tmpX,x,CovarCols);
            tc(:,zz) = mean(cc,2);
        end
    else
        tc = [];
        for zz = 1:length(Con.Groups)
            if size(Con.Groups{zz},1)>1
                tmpX = Con.Groups{zz};
            else
                tmpX = x(:,Con.Groups{zz});
            end
            cc = MakePreCons(tmpX,x,CovarCols);
            tc(:,zz) = mean(cc,2);
        end
    end
    
    levs = Con.Levs(end:-1:1);
    if levs == 0
        c = tc;
    else
        tmp = reshape(tc,[size(tc,1) levs]);
        tmp = squeeze(tmp);
        if size(tmp,2)==1
            c = tmp;
        else
            c = differencer(tmp);
            if numel(levs)==1
                c = c*-1;
            end
        end
    end
else
    [r,c,cc] = MakeContrastMatrix('a',Con.Groups,Con.Levs,x,[],CovarCols);
end
end

function out = differencer(tmp,count)

if nargin == 1;
    count = numel(size(tmp));
end

out = diff(tmp,1,count);

if count ~= 2;
    out = differencer(out,count-1);
else
    if numel(size(out))>2
        ss = size(out);
        out = reshape(out,size(out,1),prod(ss(2:end)));
    end
    return
end
end