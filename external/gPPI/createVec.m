function [vec,left,right] = createVec(left,right,SPM,Weighted,Prefix,contrail,trialcountmin,trialcountminper)
% This function generates the contrast vectors to be estimated in SPM (e.g.
% SPM.xCon.c)
%
% vec is an output vector that used in SPM.xCon.c 
%
% left cell string of tasks on the left side of the equation (e.g. {'c1'
%   'c2'} for c1+c2>c3+c4; alternatively, left can be a vector (e.g. {[1 0 0
%   1 0 0 1 0 0]}) if the contrast is known).
%
% right cell string of tasks on the left side of the equation (e.g. {'c3'
%   'c4'} for c1+c2>c3+c4; if left is a vector, then right should be empty
%
% left and right can now contain regular expressions. User 
%    must always specify the last character of the regular expression AND 
%    if you have parametric modulators, specify the entire part after the 
%    'x' or it will grab multiple tasks. Include the 'x'. 
%    ALSO SEE WARNING BELOW.
%
% SPM is the SPM.mat structure or an SPM file
%
% Weighted specifies the cut of between trial weighting and run averaging;
%   default is 0 for trial averaging everything. The variable is specified 
%   in either seconds or scans depending on SPM.mat file. This weights the 
%   each column by the number of trials in that column. NOTE: If left is a
%   vector, Weighted does not get used.
%
% contrail specifies the trailing end of the contrast names. It must be a
%   cell array and needs only to be specified if there are multiple columns
%   that share the same name (e.g. derivatives and parametric modulation)
%   (e.g. {'bf1'})
%
% Prefix is either a cell string OR a structure with .Left and .Right
%   fields. The prefix is added to the front of the task in the search for
%   the columns to use to build the contrasts. (e.g. Prefix.Left='Sn(1)'
%   and Prefix.Right='Sn(2)' to compare run1 to run2)
%
% trialcountmin is the minimum number of trials needed on the left and
%   right side for the contrast to be computed. The default is 30 trials.
%
% Note: If you are averaging condition A and condition B (A+B/2), and
% either condition has less than trialcountmin/2. the contrast will not be
% computed. If a condition does not exist, the contrast will not be
% computed.
%
% ************************************************************************
% WARNING: Just because the contrast exists, doesn't mean it is a valid
%          contrast. USER BEWARE!!!!
% ************************************************************************
%
% ************************************************************************
% WARNING: The use of regular expression for task selection may result in 
%          different tasks in different subjects if not all subjects had 
%          the same task/conditions. For example, if you have two 
%          conditions, RH and LH, and you use the expression '.*H' and one
%          subject doesn't have LH, the contrast will be created, but will
%          be different than all other subjects. USER BEWARE!!!!
% ************************************************************************
%
% License:
%   Copyright (c) 2011, Donald G. McLaren and Aaron Schultz
%   All rights reserved.
%
%    Redistribution, with or without modification, is permitted provided that the following conditions are met:
%    1. Redistributions must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%    2. All advertising materials mentioning features or use of this software must display the following acknowledgement:
%        This product includes software developed by the Harvard Aging Brain Project.
%    3. Neither the Harvard Aging Brain Project nor the
%        names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%    4. You are not permitted under this Licence to use these files
%        commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all 	
%        or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use 	
%        of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third 	
%        party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale 
%        or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received.
%
%   THIS SOFTWARE IS PROVIDED BY DONALD G. MCLAREN (mclaren@nmr.mgh.harvard.edu) AND AARON SCHULTZ (aschultz@nmr.mgh.harvard.edu)
%   ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
%   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
% Written by Aaron Schultz(aschultz@martinos.org) and Donald McLaren 
%    (mclaren@nmr.mgh.harvard.edu)
% 
% GRECC, Edith Norse Roers Memorial Veterans Hospital, Bedford, MA
% Department of Neurology, Massachusetts General Hospital and Harvard
%   Medical School
%
% v4.3_8_2011
% -Fixed contrasts for Mixed-block event-related designs
%
% v4.5_11_2011
% -Fixed prefixs for run numbers
% -Added the ability to use regular expressions to define tasks, always
% specify the last character of the regular expression AND if you have
% parametric modulators, specify the entire part after the x or it will
% grab multiple tasks
% -Conditions are always weighted equal, trial averaging is only across
% runs
% -Added a minimum trial count for event-related designs (defult is 30)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Suppress MATLAB Code Messages
%#ok<*ASGLU,*NASGU,*AGROW,*WNTAG,*CTCH,*TRYNC,*ST2NM>

%% Error checking
% Check for SPM file
if ~isstruct(SPM) && exist(SPM,'file')  
    load(SPM);
end

% Initialize vec using SPM file
try
    vec = zeros(size(left,1),length(SPM.xX.name));
catch 
    error('SPM file not correct.');
end

% Get Column Info from SPM file
for ii=1:numel(SPM.Sess)
    temp{ii}=SPM.Sess(ii).col;
    temp2(ii)=numel(temp{ii});
end
SPM.col=zeros(numel(SPM.Sess),max(temp2));
for ii=1:numel(SPM.Sess)
    SPM.col(ii,1:temp2(ii))=[temp{ii}];
end
clear temp temp2

% The next set of error checks either define default values or exit the
% program. Display messages will define the issue.
if isempty(Weighted) || ~isnumeric(Weighted)
    Weighted=0;
end

if ~isempty(left) && ~iscellstr(left) && ~iscell(left) && ~(isnumeric(left) && isempty(right))
    disp('Left contrasts are not in a cell array.')
    return;
end
if isempty(left)
    disp('Left contrasts do not exist.')
    return;
end
if ~isempty(right) && ~iscellstr(right) && ~iscell(right) && ~(isnumeric(left) && isempty(right))
    disp('Right contrasts are not in a cell array.')
    return;
end
if exist('contrail','var') && isempty(contrail)
    clear contrail
end
if exist('contrail','var') && ~iscellstr(contrail)
    contrail={contrail};
end
if exist('contrail','var')
    if size(contrail,1)~=size(left,1) 
        if size(contrail,1)==1
            repmat(contrail,size(left,1),1)
        else
            error('Contrail is not correct. Must either be a string array OR cell array with the same number of rows as left/right.')
        end
    end
end
if size(left,1)~=size(right,1) && ~(isnumeric(left) && isempty(right))
    error('F-test vectors not the same size on the left and right sides.');
end

if ~exist('trialcountmin','var')
    trialcountmin=30; %Donald McLaren's minimum trial count based on Huettel et al. 2001.
end
if trialcountmin<=0
    warning('You cannot exclude conditions, all conditions must be present. trialcountmin must be greater than 0.')
    vec(all(vec==0,2),:)=[]; 
    return
end

%% Program begins here
if isempty(right) % This is in the case that you specified the contrast as as a matrix of numbers
   try
       if length(left)~=length(SPM.xX.name)
           if length(SPM.xX.name)/length(left)~=numel(SPM.Sess)
                vec=[left zeros(1,length(SPM.xX.name)-length(left))];
                vec(vec>0)=vec(vec>0)/sum(vec(vec>0));
                vec(vec<0)=vec(vec<0)/(-1*sum(vec(vec<0)));
           else
                vec=repmat(left,1,numel(SPM.Sess));
                vec(vec>0)=vec(vec>0)/sum(vec(vec>0));
                vec(vec<0)=vec(vec<0)/(-1*sum(vec(vec<0)));
           end
       else
           if isnumeric(left)
                vec=left;
                vec(vec>0)=vec(vec>0)/sum(vec(vec>0));
                vec(vec<0)=vec(vec<0)/(-1*sum(vec(vec<0)));
           else
               vec=zeros(1,length(SPM.xX.name));
               disp('Invalid Contrast -- Not sure how this happened.')
           end
       end
   catch
       vec=zeros(1,length(SPM.xX.name));
   end
else
    % Setup parameters
    tasks={};
    dur=[];
    trials=[];
    sessindex=[];
    bk=0;
    % Find all task names, durations, and trial numbers from SPM.mat
    try 
        for ss=1:numel(SPM.Sess)
            for tt=1:numel(SPM.Sess(ss).U)
                tasks=[tasks SPM.Sess(ss).U(tt).name];
                dur(end+1:end+numel(SPM.Sess(ss).U(tt)))=repmat(max(SPM.Sess(ss).U(tt).dur)<Weighted,numel(SPM.Sess(ss).U(tt)),1);
                trials(end+1:end+numel(SPM.Sess(ss).U(tt)))=repmat(length(SPM.Sess(ss).U(tt).ons),numel(SPM.Sess(ss).U(tt)),1);
                sessindex(end+1:end+numel(SPM.Sess(ss).U(tt)))=repmat(ss,numel(SPM.Sess(ss).U(tt)),1);
            end
        end
        [Tasks, junk, taskind]=unique(tasks);
    catch 
        disp('WARNING: Could not find SPM variable. Where is it?')
        return
    end
    % Grab the tasks for the left and right sides of the contrasts, make sure
    % they exist. If a regexp was used, this will expand the task list  
    left=tasklist(left,Tasks);
    right=tasklist(right,Tasks);
     
    %Check that contrast is valid based on tasklist feedback
    if any(any(strcmp(left,'contrastnotfound'))) || any(any(strcmp(right,'contrastnotfound')))
        disp('Invalid Contrast')
        vec(all(vec==0,2),:)=[];
        return
    end

    %Check the both sides are unique
    for vecrow=1:size(right,1)
        if numel(unique([left(vecrow,:) right(vecrow,:)]))~=numel([left(vecrow,:) right(vecrow,:)])
            disp(['WARNING: The same tasks appear on both sides of the equation in contrast ' num2str(vecrow) '. This is not valid.'])
            disp('Invalid Contrast')
            vec(all(vec==0,2),:)=[];
            return
        end
    end
    
    %Make contrail valid for use with regexp
    try
           for ii=1:numel(contrail)
                contrail{ii}=regexprep(contrail{ii},')','\\)');
                contrail{ii}=regexprep(contrail{ii},'\^','\\^');
           end
    end
 
    % Left side of vec
    for vecrow=1:size(left,1)
        tce=0; %set events (Yes/No) to No
        tcb=0; %set blocks (Yes/No) to No [also for session averaging]
        Lind={}; %initialize left index
        LWeightb={}; %initialize left weighting for blocks [also for session averaging or trials]
        LWeight={}; %initialize left weighting for events
        vecL=vec(vecrow,:); %initialize left vector for events
        vecLb=vec(vecrow,:); %initialize left vector for blocks [also for session averaging or trials]
        bk=0; %set break loop (Yes/No) variable to No.
        for ii = 1:size(left,2);
            % Check if left side is none/empty/empty; break for loop if true
            if strcmpi(left{vecrow,ii},'none') || strcmpi(left{vecrow,ii},'None') || isempty(left{vecrow,ii})
                if ii~=1
                    bk=1;
                    break;
                else
                    Lind{ii}=[];
                    LWeight{ii}=[];
                    break;
                end
            end
            % Check for spaces in tasklist, skip if found
            if strcmpi(left{vecrow,ii},'') || strcmpi(left{vecrow,ii},' ')
                Lind{ii}=[];
                LWeight{ii}=[];
                continue;
            end
            % Change task for use in regexp with * ^ symbols
            left{vecrow,ii}=regexprep(left{vecrow,ii},'\*','\\*');
            left{vecrow,ii}=regexprep(left{vecrow,ii},'\^','\\^');
            % Check for prefix
            try
                if iscellstr(Prefix)
                    prefix=Prefix;
                elseif isstruct(Prefix)
                    if isfield(Prefix,'Left') && isfield(Prefix,'Right')
                        prefix=Prefix.Left;
                    else
                        bk=1;
                        break;
                    end
                elseif isempty(Prefix)
                    prefix={' '};
                else
                    prefix={Prefix};
                end
            catch
                prefix={' '};
            end
            prefix=expandprefix(prefix);
            
            % This strips off any prefix in the task list.
            try
                for i=1:length(Tasks)
                    if regexp(left{vecrow,ii},[Tasks{i} '$'])
                        prefix=strcat(prefix,{left{vecrow,ii}(1:regexp(left{vecrow,ii},[Tasks{i} '$'])-1)});
                        left{vecrow,ii}=left{vecrow,ii}(regexp(left{vecrow,ii},[Tasks{i} '$']):end);
                        break;
                    end
                end
            end
            % Change more characters for regexp and searching
            left{vecrow,ii}=regexprep(left{vecrow,ii},')','\\)');
            try
                for i=1:numel(prefix)
                    prefix{i}=regexprep(prefix{i},')','\\)');
                end
            end
            % Check if you want to average by runs (trialaveraging=0) or
            % trials (trialaveraging=1)
            try
                jj=contains([left{vecrow,ii} '$'],Tasks);
                trialaveraging=max(dur(taskind==jj));
            catch
                trialaveraging=0;
            end
            % Find the columns that correspond to each task. Put into the
            % *b variables for session averaging, non-*b variable for trial
            % averaging
            if trialaveraging
                ind = [];
                for pp=1:max(numel(prefix),1)
                    try
                        curprefix=prefix{pp};
                    catch
                        curprefix={' '};
                    end
                    try 
                        curcontrail=contrail(vecrow,:);
                    catch
                        curcontrail={'.*'};
                    end
                    try
                        ind=findcolumnindex(SPM,curprefix,left{vecrow,ii},curcontrail,SPM.xX.name);
                    catch
                        bk=1;
                        break;
                    end
                    try
                        Lind{ii} = [Lind{ii} ind];
                    catch
                        Lind{ii} = ind;
                    end
                    if ~isempty(ind)
                        try
                            LWeightTmp = repmat(trials(taskind==jj),1,numel(contrail));
                            sessindexTmp = repmat(sessindex(taskind==jj),1,numel(contrail));
                        catch
                            LWeightTmp = trials(taskind==jj);
                            sessindexTmp = sessindex(taskind==jj);
                        end
                    else
                        LWeightTmp=[];
                    end
                    if numel(LWeightTmp)~=numel(ind) && ~isempty(contains({'^Sn('},curprefix)) && contains({'^Sn('},curprefix)==1
                        sessind=find(sessindexTmp==str2double(curprefix(4)));
                        if ~isempty(sessind)
                            try
                                LWeight{ii}=[LWeight{ii} LWeightTmp(sessind)];
                            catch
                                LWeight{ii}=LWeightTmp(sessind);
                            end
                        end
                    elseif  numel(LWeightTmp)~=numel(ind)
                        bk=1;
                        break;
                    else
                        try
                            LWeight{ii}=[LWeight{ii} LWeightTmp];
                        catch
                            LWeight{ii}=LWeightTmp;
                        end
                    end
                end
                if numel(LWeight)~=ii || isempty(LWeight{ii})
                    tce(ii)=-1;
                else
                    tce(ii)=sum(LWeight{ii});
                end
            else
                if ~isempty(LWeightb)
                    TaskCount=TaskCount+1;
                else
                    TaskCount=1;
                end
                ind = [];
                for pp=1:max(numel(prefix),1)
                    try
                        curprefix=prefix{pp};
                    catch
                        curprefix={' '};
                    end
                    try 
                        curcontrail=contrail(vecrow,:);
                    catch
                        curcontrail={'.*'};
                    end
                    try
                        ind=findcolumnindex(SPM,curprefix,left{vecrow,ii},curcontrail,SPM.xX.name);
                    catch
                        bk=1;
                        break;
                    end
                    try
                        Lind{ii} = [Lind{ii} ind];
                    catch
                        Lind{ii} = ind;
                    end
                    LWeightTmpb = ones(1,length(ind));
                    try
                        LWeightb{ii}=[LWeightb{ii} LWeightTmpb];
                    catch
                        LWeightb{ii}=LWeightTmpb;
                    end
                end
                if isempty(LWeightb{ii})
                    tcb(ii)=-1;
                else
                    tcb(ii)=sum(LWeightb{ii});
                end
            end
        end
        
        % Check that the minimum number of trials 
        try
            trialminimum=trialcountminper;
        catch
            trialminimum=trialcountmin./sum(tce>0); %Only count event-related effects
        end
        if any(tce==-1) || any(tcb==-1)
            disp('Warning: Missing conditions!!! Invalid Contrast')
            bk=1;
        end
        if (~all(tce(tce>0)>trialminimum) || sum(tce)<trialcountmin) && all(tce~=0)
            disp('Warning: Not enough events!!!')
            bk=1;
        end
        if ~isempty(LWeightb) && TaskCount>1 && ~isempty(LWeight)
            bk=1;
            disp('You cannot have multiple block types on the same side with events!!!')
        end
        % Create vector for contrast
        if ~bk
            jj=[];
            for ii = 1:length(Lind);
                try
                    vecL(Lind{ii}) = LWeight{ii}./tce(ii)./sum(tce>0);
                catch
                    jj=[jj ii];
                end
            end
            for ii=jj
                if sum(vecL(:))~=0
                    for colnum=1:numel(Lind{ii})
                        [row,col]=find(SPM.col==Lind{ii}(colnum));
                        vecLb(Lind{ii}(colnum)) = LWeightb{ii}(colnum).*sum(vecL(SPM.col(row,SPM.col(row,:)>0)));
                    end
                else
                    vecLb(Lind{ii}) = LWeightb{ii}./tcb(ii)/sum(tcb>0);
                end
            end
            vecL=vecL+vecLb;
            if sum(vecL)==0 && ~strcmpi(left{vecrow,1},'none')
                bk=1;
            end
        end

        % Right side (see previous comments in left side as the code is
        % identical except the R is used instead of L.
        vecR=vec(vecrow,:);
        vecRb=vec(vecrow,:);
        tce=0;
        tcb=0;
        Rind={};
        RWeightb={};
        RWeight={};
        for ii = 1:size(right,2);
            if strcmpi(right{vecrow,ii},'none') || strcmpi(right{vecrow,ii},'None')
                if ii~=1
                    bk=1;
                    break;
                else
                    Rind{ii}=[];
                    RWeight{ii}=[];
                    break;
                end
            elseif bk
                break;
            end
            if strcmpi(right{vecrow,ii},'') || strcmpi(right{vecrow,ii},' ')
                Rind{ii}=[];
                RWeight{ii}=[];
                continue;
            end
            right{vecrow,ii}=regexprep(right{vecrow,ii},'\*','\\*');
            right{vecrow,ii}=regexprep(right{vecrow,ii},'\^','\\^');
            try
                if iscellstr(Prefix)
                    prefix=Prefix;
                elseif isstruct(Prefix)
                    if isfield(Prefix,'Left') && isfield(Prefix,'Right')
                        prefix=Prefix.Right;
                    else
                        bk=1;
                        break;
                    end
                elseif isempty(Prefix)
                    prefix={' '};
                else
                    prefix={Prefix};
                end
            catch
                prefix={' '};
            end
            prefix=expandprefix(prefix);
            
            try
                for i=1:length(Tasks)
                    if strfind(right{vecrow,ii},Tasks{i})
                        prefix=strcat(prefix,{right{vecrow,ii}(1:strfind(right{vecrow,ii},Tasks{i})-1)});
                        right{vecrow,ii}=right{vecrow,ii}(strfind(right{vecrow,ii},Tasks{i}):end);
                        break;
                    end
                end
            end
            right{vecrow,ii}=regexprep(right{vecrow,ii},')','\\)');
            try
                for i=1:numel(prefix)
                    prefix{i}=regexprep(prefix{i},')','\\)');
                end
            end
            try
                jj=contains(right{vecrow,ii},Tasks);
                trialaveraging=max(dur(taskind==jj));
            catch
                trialaveraging=0;
            end
            if trialaveraging
                ind = [];
                for pp=1:max(numel(prefix),1)
                    try
                        curprefix=prefix{vecrow,pp};
                    catch
                        curprefix={' '};
                    end
                    try 
                        curcontrail=contrail(vecrow,:);
                    catch
                        curcontrail={'.*'};
                    end
                    try
                        ind=findcolumnindex(SPM,curprefix,right{vecrow,ii},curcontrail,SPM.xX.name);
                    catch
                        bk=1;
                        break;
                    end
                    try
                        Rind{ii} = [Rind{ii} ind];
                    catch
                        Rind{ii} = ind;
                    end
                    if ~isempty(ind)
                        try
                            RWeightTmp = repmat(trials(taskind==jj),1,numel(contrail));
                            sessindexTmp = repmat(sessindex(taskind==jj),1,numel(contrail));
                        catch
                            RWeightTmp = trials(taskind==jj);
                            sessindexTmp = sessindex(taskind==jj);
                        end
                    else
                        RWeightTmp=[];
                    end
                    if numel(RWeightTmp)~=numel(ind) && ~isempty(contains({'^Sn('},curprefix)) && contains({'^Sn('},curprefix)==1
                        sessind=find(sessindexTmp==str2double(curprefix(4)));
                        if ~isempty(sessind)
                            try
                                RWeight{ii}=[RWeight{ii} RWeightTmp(sessind)];
                            catch
                                RWeight{ii}=RWeightTmp(sessind);
                            end
                        end
                    elseif  numel(RWeightTmp)~=numel(ind)
                        bk=1;
                        break;
                    else
                        try
                            RWeight{ii}=[RWeight{ii} RWeightTmp];
                        catch
                            RWeight{ii}=RWeightTmp;
                        end
                    end
                end
                if numel(RWeight)~=ii || isempty(RWeight{ii})
                    tce(ii)=-1;
                else
                    tce(ii)=sum(RWeight{ii});
                end
            else
                if ~isempty(RWeightb)
                    TaskCount=TaskCount+1;
                else
                    TaskCount=1;
                end
                ind = [];
                for pp=1:max(numel(prefix),1)
                    try
                        curprefix=prefix{pp};
                    catch
                        curprefix={' '};
                    end
                    try 
                        curcontrail=contrail(vecrow,:);
                    catch
                        curcontrail={'.*'};
                    end
                    try 
                        ind=findcolumnindex(SPM,curprefix,right{vecrow,ii},curcontrail,SPM.xX.name);
                    catch
                        bk=1;
                        break;
                    end
                    try
                        Rind{ii} = [Rind{ii} ind];
                    catch
                        Rind{ii} = ind;
                    end
                    RWeightTmpb = ones(1,length(ind));
                    try
                        RWeightb{ii}=[RWeightb{ii} RWeightTmpb];
                    catch
                        RWeightb{ii}=RWeightTmpb;
                    end
                end
                if isempty(RWeightb{ii})
                    tcb(ii)=-1;
                else
                    tcb(ii)=sum(RWeightb{ii});
                end
            end
        end
    
        try
            trialminimum=trialcountminper;
        catch
            trialminimum=trialcountmin./sum(tce>0); %Only count event-related effects
        end
        if any(tce==-1) || any(tcb==-1)
            disp('Warning: Missing conditions!!! Invalid Contrast')
            bk=1;
        end
        if (~all(tce(tce>0)>trialminimum) || sum(tce)<trialcountmin) && all(tce~=0)
            disp('Warning: Not enough events!!!')
            bk=1;
        end
        if ~isempty(RWeightb) && TaskCount>1 && ~isempty(RWeight)
            bk=1;
            disp('You cannot have multiple block types on the same side with events!!!')
        end
        if ~bk
            jj=[];
            for ii = 1:length(Rind);
                try
                    vecR(Rind{ii}) = RWeight{ii}./tce(ii)./sum(tce>0);
                catch
                    jj=[jj ii];
                end
            end
            for ii=jj
                if sum(vecR(:))~=0
                    for colnum=1:numel(Rind{ii})
                        [row,col]=find(SPM.col==Rind{ii}(colnum));
                        vecRb(Rind{ii}(colnum)) = RWeightb{ii}(colnum).*sum(vecR(SPM.col(row,SPM.col(row,:)>0)));
                    end
                else
                    vecRb(Rind{ii}) = RWeightb{ii}./tcb(ii)/sum(tcb>0);
                end
            end
            vecR=vecR+vecRb;
            if sum(vecR)==0 && ~strcmpi(right{vecrow,1},'none')
                bk=1;
            end
        end
        
        %Combine Left and Right Sides
        if  bk==1 || (sum(vecR)~=0 && sum(vecL)~=0 && round(sum(vecR))~=round(sum(vecL)))
            disp('Invalid Contrast')
        else
            disp('Valid Contrast')
            vec(vecrow,:)=vecL-vecR;
        end
    end
    vec(all(vec==0,2),:)=[];
end
vec(all(vec==0,2),:)=[];
end

%% contains - this is a regexp search function
function ind = contains(str, list)
which = regexp(list,str);
ind = [];
for kk = 1:length(which);
    if ~isempty(which{kk})
        ind(end+1) = kk;
    end
end
end

%% findcolumnindex -- finds the column of the task in SPM.xX.names
% Checks for all possible basis functions as well.
function ind=findcolumnindex(SPM,prefix,task,suffix,list)

attempt=1;
ind=[];
while attempt<5 && isempty(ind)
    switch attempt
        case {1}
            Task=[task '$'];     
        case {2}
            Task=[task '\*'];
            if strcmpi(suffix,'.*')
                a=1;
                suffix=repmat({'bf(1\)'},1,numel(suffix));
            end
        case {3}
            if exist('a','var')
                suffix=repmat([],1,numel(suffix));
            end
            Suffix=suffix;
        case {4}
            suffix=Suffix;
    end
    for kk=1:size(suffix,2)
        if attempt==3
            Task=[task suffix{kk} '$'];
            suffix{kk}=[];
        elseif attempt==4
            Task=[task suffix{kk} '$'];
            suffix{kk}=[];
            Task=[Task(1:end-1) '\*bf(1\)$']; 
        end
        ind=[ind contains([prefix ' ?' Task suffix{kk}],list)];
    end
    if length(ind)>(numel(SPM.Sess)*numel(suffix)) 
            err='Invalid Contrast -- too many columns selected';
            ind=[];
    else
        err=[];
    end
    attempt=attempt+1;
end
if ~isempty(err)
    disp(err)
    error(Thiswillbreak) 
end
end

%% Find Task List -- takes the left/right task definitions and expands them if they are regexp and checks them against the master Task list
function taskout=tasklist(taskin,Tasks)
for vecrow=1:size(taskin,1)
    for tt=1:size(taskin,2)
        if ~isempty(taskin{vecrow,tt}) && ~strcmpi(taskin{vecrow,tt},'none') && ~strcmpi(taskin{vecrow,tt},'None')
            taskin{vecrow,tt}=regexprep(taskin{vecrow,tt},'\*','\\*');
            taskin{vecrow,tt}=regexprep(taskin{vecrow,tt},'\^','\\^');
            taskin{vecrow,tt}=regexprep(taskin{vecrow,tt},'\\\','\\');
            taskin{vecrow,tt}=regexprep(taskin{vecrow,tt},'.\*','.*');
            ind = contains(['^' taskin{vecrow,tt} '$'], Tasks);
        else
            ind=[];
        end
        if ~isempty(ind)
            try
                taskintmp=[taskintmp (Tasks(ind))];
            catch
                taskintmp=Tasks(ind);
            end
        elseif strcmpi(taskin{vecrow,tt},'none') || strcmpi(taskin{vecrow,tt},'None') || isempty(taskin{vecrow,tt})
            try
                taskintmp=[taskintmp 'none'];
            catch
                taskintmp={'none'};
            end
        else
            disp(['WARNING: Task ''' taskin{vecrow,tt} ''' not found. Contrast will not be computed. Please check the tasks.'])
            %disp('Invalid Contrast')
            taskintmp={'contrastnotfound'};
            break
        end
    end
    taskintmp=unique(taskintmp); %Removes duplicates due to human error and/or specifying a wild card and the task.
    for tt=1:numel(taskintmp)
        taskinTmp{vecrow,tt}=taskintmp{1,tt};
    end
    clear taskintmp
end
taskout=taskinTmp; clear taskinTmp
end

%% Expand Prefix -- expands the prefix for multiple sessions
function prefix=expandprefix(prefix)
for pii=1:numel(prefix)
    if ~strcmp(prefix(pii),' ') && ~isempty(contains('Sn([', prefix(pii)))
        s1 = find(prefix{pii}=='[');
        s2 = find(prefix{pii}==']');
        range = prefix{pii}(s1+1:s2-1);
        if contains('-',{range});
            tmp = regexp(range,'-','split');
            sessions = str2num(tmp{1}):str2num(tmp{2});
        else
            tmp = regexp(range,',','split');
            sessions = str2num(char(tmp'))';
        end
        for pp=1:numel(sessions)
            sessionsTmp{pp}=[prefix{pii}(1:s1-1) num2str(sessions(pp)) prefix{pii}(s2+1:numel(prefix{pii}))];
        end
        try
            prefixTmp=[prefixTmp sessionsTmp]; clear sessionsTmp
        catch
            prefixTmp=sessionsTmp; clear sessionsTmp
        end
    else
        try
            prefixTmp=[prefixTmp prefix(pii)];
        catch
            prefixTmp=prefix(pii);
        end
    end
end
prefix=prefixTmp; clear prefixTmp
end
