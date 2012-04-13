function SPM=spm_contrasts_PPI(SPM,Subject,Contrasts,Weighted,Method)
%Estimates PPI contrasts
%   SPM is the SPM file for PPI or the SPM structure for PPI after
%   estimation
%   Contrasts is a structure of the contrasts to be applied to the PPI
%   estimation, if not specified, computes all possible combinations
%   Weighted is numeric cutoff for the minimum time of an event/epoch/block
%   to use run averaging as opposed to trial averaging. Default is 0, so
%   all events will be analyzed using run averaging.
%   Method is a string that is either 'trad' or 'cond'.
%
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
%   Last modified on 11/16/2010 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%   Medical School

%% Program begins here
%Check input arguments
if ~isstruct(SPM)
    if exist(SPM,'file')==2
        load(SPM)
    else
        disp('PPI not estimated')
        return;
    end
end
cd(SPM.swd)
if nargin==1
    error('Subject must be specified')
end
if nargin==2 || isempty(Contrasts)
    try
        Contrasts=defContrasts(SPM,Weighted);
    catch
        Contrasts=defContrasts(SPM,0);
    end
end
if nargin<4
    Weighted=0;
end
if ~isempty(Contrasts) && iscellstr(Contrasts)
    try
        Contrasts=defContrasts(SPM,Weighted,0,Contrasts);
    catch
        Contrasts=defContrasts(SPM,0,0,Contrasts);
    end
end
if ~isfield(Contrasts,'Weighted')
    Contrasts(1).Weighted=Weighted;
end

%Configure Contrasts
ind=zeros(length(Contrasts),1);
for ii = 1:length(Contrasts)
    if ~isfield(Contrasts(ii),'c') || isempty(Contrasts(ii).c)
        if isempty(Contrasts(ii).Weighted)
            Contrasts(ii).Weighted=Weighted;
        end
        if strcmpi(Method,'trad')
            try
                Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,Contrasts(ii).Prefix)';
            catch
                Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted)';
            end
        else
            if isfield(Contrasts(ii),'Prefix') && ~isempty(Contrasts(ii).Prefix) && isfield(Contrasts(ii),'Contrail') && ~isempty(Contrasts(ii).Contrail)
                try
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,[Contrasts(ii).Prefix 'PPI_'],Contrasts(ii).Contrail,Contrasts(ii).MinEvents,Contrasts(ii).MinEventsPer)';
                catch
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,[Contrasts(ii).Prefix 'PPI_'],Contrasts(ii).Contrail,Contrasts(ii).MinEvents)';
                end
                  
            elseif isfield(Contrasts(ii),'Prefix') && ~isempty(Contrasts(ii).Prefix)
                try
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,[Contrasts(ii).Prefix 'PPI_'],[],Contrasts(ii).MinEvents,Contrasts(ii).MinEventsPer)';
                catch
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,[Contrasts(ii).Prefix 'PPI_'],[],Contrasts(ii).MinEvents)';
                end
            elseif isfield(Contrasts(ii),'Contrail') && ~isempty(Contrasts(ii).Contrail)
                try
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,'PPI_',Contrasts(ii).Contrail,Contrasts(ii).MinEvents,Contrasts(ii).MinEventsPer)';
                catch
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,'PPI_',Contrasts(ii).Contrail,Contrasts(ii).MinEvents)';
                end
            else
                try
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,'PPI_',[],Contrasts(ii).MinEvents,Contrasts(ii).MinEventsPer)';
                catch
                    Contrasts(ii).c=createVec(Contrasts(ii).left,Contrasts(ii).right,SPM,Contrasts(ii).Weighted,'PPI_',[],Contrasts(ii).MinEvents)';
                end
            end
        end
        mean(Contrasts(ii).c==0)~=1;
        if mean(Contrasts(ii).c==0)~=1; ind(ii)=1; end
        if isempty(Contrasts(ii).name)
            if isfield(Contrasts(ii),'Prefix') && ~isempty(Contrasts(ii).Prefix) && isfield(Contrasts(ii),'Contrail') && ~isempty(Contrasts(ii).Contrail)
                tmp=[Contrasts(ii).Prefix '_PPI_' Contrasts(ii).left 'minus' Contrasts(ii).right '_' Contrasts(ii).Contrail];
            elseif isfield(Contrasts(ii),'Prefix') && ~isempty(Contrasts(ii).Prefix)
                tmp=[Contrasts(ii).Prefix '_PPI_' Contrasts(ii).left 'minus' Contrasts(ii).right];
            elseif isfield(Contrasts(ii),'Contrail') && ~isempty(Contrasts(ii).Contrail)
                tmp=['PPI_' Contrasts(ii).left 'minus' Contrasts(ii).right '_' Contrasts(ii).Contrail];
            else
                tmp=['PPI_' Contrasts(ii).left 'minus' Contrasts(ii).right];
            end
            tmp=sprintf('%s_',tmp{:});
            tmp=tmp(1:end-1);
            Contrasts(ii).name=tmp;
        else
            if iscellstr(Contrasts(ii).name)
                Contrasts(ii).name=['PPI_' cell2mat(Contrasts(ii).name)];
            else
                Contrasts(ii).name=['PPI_' Contrasts(ii).name];
            end
        end
        if isempty(Contrasts(ii).STAT)
            Contrasts(ii).STAT='T';
        end
    end
end
Contrasts=Contrasts(ind==1);
for ii = 1:length(Contrasts)
    xCon(ii) = spm_FcUtil('Set',Contrasts(ii).name,Contrasts(ii).STAT,'c',Contrasts(ii).c,SPM.xX.xKXs);
end

%Compute Contrasts
try
    init=length(SPM.xCon);
catch
    init=0;
end
if init~=0
    SPM.xCon(init+1:init+length(xCon)) = xCon;
else
    SPM.xCon = xCon;
end
SPM = spm_contrasts(SPM,init+1:length(SPM.xCon));

% Move contrasts
for ii=(1+init):numel(SPM.xCon)
    disp('Moving Contrast Images');
    DD=['_' date];
    f1 = SPM.xCon(ii).Vcon.fname;
    f2 = [SPM.xCon(ii).Vcon.fname(1:end-3) 'hdr'];
    [junk fname]=fileparts(SPM.xCon(ii).Vcon.fname); fname=fname(6:end); clear junk
    if length(SPM.xCon(ii).name)==4 && strcmp(fname,SPM.xCon(ii).name);
        disp(['Contrast: ' SPM.xCon(ii).Vcon.fname ' was not moved.'])
        continue
    end
    if ~exist([f1(1:4) SPM.xCon(ii).name '_' Subject '.img'],'file')
        try
            movefile(f1, [f1(1:4) SPM.xCon(ii).name '_' Subject '.img'],'f');
            movefile(f2, [f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'],'f');
        catch
            disp(['error moving contrast '  SPM.xCon(ii).name])
        end
    else
        try
            movefile([f1(1:4) SPM.xCon(ii).name '_' Subject '.img'], [f1(1:4) SPM.xCon(ii).name DD '_' Subject  '.img'],'f');
            movefile([f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'], [f2(1:4) SPM.xCon(ii).name DD '_' Subject  '.hdr'],'f');
            movefile(f1, [f1(1:4) SPM.xCon(ii).name '_' Subject '.img'],'f');
            movefile(f2, [f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'],'f');
            for kk = 1:numel(SPM.xCon)
                if strcmp(SPM.xCon(kk).Vcon.fname,[f1(1:4) SPM.xCon(ii).name '_' Subject '.img']) > 0
                    SPM.xCon(kk).Vcon.fname=[f1(1:4) SPM.xCon(ii).name DD '_' Subject '.img'];
                    SPM.xCon(kk).name=[SPM.xCon(ii).name DD];
                    break
                end
            end
        catch
            disp(['error moving contrast '  SPM.xCon(ii).name])
        end
    end
    if ~strcmp(SPM.xCon(ii).Vcon.fname,[f1(1:4) SPM.xCon(ii).name '_' Subject '.img'])
        SPM.xCon(ii).Vcon.fname=[f1(1:4) SPM.xCon(ii).name '_' Subject '.img'];
    end
    f1 = SPM.xCon(ii).Vspm.fname;
    f2 = [SPM.xCon(ii).Vspm.fname(1:end-3) 'hdr'];
    if ~exist([f1(1:4) '_' SPM.xCon(ii).name '_' Subject '.img'],'file')
       try
            movefile(f1, [f1(1:4) '_' SPM.xCon(ii).name '_' Subject '.img'],'f');
            movefile(f2, [f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'],'f');
       catch
           disp(['error moving contrast '  SPM.xCon(ii).name])
       end
    else
        try
            movefile([f1(1:4) '_' SPM.xCon(ii).name '_' Subject '.img'], [f1(1:4) '_' SPM.xCon(ii).name DD '_' Subject '.img'],'f');
            movefile([f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'], [f2(1:4) '_' SPM.xCon(ii).name DD '_' Subject '.hdr'],'f');
            movefile(f1, [f1(1:4) '_' SPM.xCon(ii).name '_' Subject '.img'],'f');
            movefile(f2, [f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'],'f');
            for kk = 1:numel(SPM.xCon)
                if strcmp(SPM.xCon(kk).Vspm.fname,[f1(1:4) SPM.xCon(ii).name '_' Subject '.img']) > 0
                    SPM.xCon(kk).Vspm.fname=[f1(1:4) SPM.xCon(ii).name DD '_' Subject '.img'];
                    break
                end
            end
       catch
           disp(['error moving contrast '  SPM.xCon(ii).name])
        end
    end 
    if ~strcmp(SPM.xCon(ii).Vspm.fname,[f1(1:4) '_' SPM.xCon(ii).name '_' Subject '.img'])
        SPM.xCon(ii).Vspm.fname=[f1(1:4) '_' SPM.xCon(ii).name '_' Subject '.img'];
    end 
end
save SPM.mat SPM
    


    