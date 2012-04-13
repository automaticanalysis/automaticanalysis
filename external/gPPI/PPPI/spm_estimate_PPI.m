function SPM=spm_estimate_PPI(Subject,SPM,Regions,Method,Analysis,Contrast)
%Estimates 1st level PPI model
%   Subject is the subject number 
%   SPM is the 1st level activity model, used to pull values from
%   Regions are the PPI seed regions
%   Method is 'cond' or 'trad'
%   Analysis is either 'psy', 'phy', or 'psyphy'
%   Contrast is a computational method variable:
%       0 not to estimate any contrasts;
%       1 to estimate contrasts;  
%       2 to only use PPI txt file for 1st level (not recommended); 
%       3 to only use PPI txt file for 1st level and estimate contrasts (not recommended);
%           NOTE: 2&3 are not recommended as they potentially do not include
%                 all tasks effects in the mode. Use at your own risk.
%           NOTE: Not required, defaults to act like 0.
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

%% Check parameters and setup variables
if exist(SPM,'file')
    load(SPM)
elseif exist(SPM,'var')
else
    disp('Program will now exit. No SPM.mat was detected.')
    return;
end

% Directory
[region1,region2]=strtok(Regions);
region2=region2(2:end);
if strcmp(Analysis,'psyphy')
    A='PPPI';
else
    A='PPI';
end
if isempty(region2)
   SPMPPI.swd=[SPM.swd filesep A '_'  region1 filesep];
else
   SPMPPI.swd=[SPM.swd filesep A '_' region1 '_' region2 filesep];
end

if exist([SPMPPI.swd 'SPM.mat'],'file')==2
    SPMPPI2=load([SPMPPI.swd 'SPM.mat']);
    if isfield(SPMPPI2.SPM,'VM') && (exist(SPMPPI2.SPM.VM.fname,'file')==2 || exist([SPMPPI.swd SPMPPI2.SPM.VM.fname],'file')==2)
        estimate=0;
        SPMPPI=SPMPPI2; clear SPMPPI2
    else
        estimate=1;
    end
else
    estimate=1;
end
   
if estimate==1
    % PPI File Extension
    try
        if Contrast==2 || Contrast==3
            ext='.txt';
        else
            ext='.mat';
        end
    catch
        ext='.mat';
    end

    % Get needed fields from task model
    SPMPPI.xBF=SPM.xBF;
    SPMPPI.xBF=SPM.xBF;
    SPMPPI.xBF.name='hrf';
    SPMPPI.xBF.order=1;
    SPMPPI.xBF.bf=SPMPPI.xBF.bf(:,1:SPMPPI.xBF.order);
    SPMPPI.xY=SPM.xY;
    SPMPPI.nscan=SPM.nscan;
    SPMPPI.SPMid=SPM.SPMid;
    SPMPPI.xVi.form=SPM.xVi.form;
    if ischar(SPM.xVi.form) && strcmp(SPM.xVi.form,'i.i.d')
        SPMPPI.xVi.form='none';
    end 
    SPMPPI.xGX.iGXcalc=SPM.xGX.iGXcalc;
    SPMPPI.xGX.sGXcalc=SPM.xGX.sGXcalc;
    SPMPPI.xGX.sGMsca=SPM.xGX.sGMsca;
    for i=1:numel(SPMPPI.nscan)
        SPMPPI.xX.K(i).HParam = SPM.xX.K(i).HParam;
    end

    for z=1:numel(SPM.Sess)
        if strcmp(ext,'.txt')
            try
                load([Subject '_' region1 '_and_' region2 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            catch
                load([Subject '_' region1 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            end
            SPMPPI.Sess(z).C.C=[OUT.P.C OUT.PPI.C OUT.Y.C OUT.C.C];
            SPMPPI.Sess(z).C.name=[OUT.P.name OUT.PPI.name OUT.Y.name OUT.C.name];
            return;
        else
            try
                load([Subject '_' region1 '_and_' region2 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            catch
                load([Subject '_' region1 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            end
            if strcmpi(Method,'trad')
                SPMPPI.Sess(z).U=[];
                SPMPPI.Sess(z).C.C=[OUT.P.C OUT.PPI.C OUT.Y.C OUT.C.C];
                SPMPPI.Sess(z).C.name=[OUT.P.name OUT.PPI.name OUT.Y.name OUT.C.name];
            else
                SPMPPI.Sess(z).U=SPM.Sess(z).U;
                SPMPPI.Sess(z).C.C=[OUT.PPI.C OUT.Y.C OUT.C.C];
                SPMPPI.Sess(z).C.name=[OUT.PPI.name OUT.Y.name OUT.C.name];
            end
        end
    end

    %% Estimate PPI 1st Level Model
    try
        cd(SPMPPI.swd)
    catch
        mkdir(SPMPPI.swd)
        cd(SPMPPI.swd)
    end
    save SPM SPMPPI

    % Delete any existing files
    delete beta_00*
    delete ResMS.*
    delete RPV.*
    delete mask.*

    % Make design and estimate, rename to SPM
    SPM1=SPM; clear SPM
    SPM=SPMPPI; clear SPMPPI
    SPM = spm_fmri_spm_ui(SPM);
    disp('estimate_PPI.m')
    SPM.xM.T(:) = -Inf;  %% disable threshold masking
    SPM.xM.TH(:)=SPM1.xM.TH(:);
    try
        V = spm_vol(SPM1.VM.fname);
    catch
        try
            V = spm_vol([SPM1.swd filesep SPM1.VM.fname]);
        catch
            disp(['Mask file cannot be found: ' SPM1.VM.fname])
            error('Mask file cannot be found.')
        end
    end
    SPM.xM.VM = V;
    SPM=spm_spm(SPM);
else
    SPM=SPMPPI.SPM;
end

