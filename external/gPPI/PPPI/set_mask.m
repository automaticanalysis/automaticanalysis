function [xY,errorchk]=set_mask(SPM,VOI,region1,maskdir)
%% Checks if xY.mask exists. 
% If VOI is an image, read the header into xY.mask.
% If xY.mask does not exist, then a mask image (.nii) is created via 
% mask_create and the the header is read into xY.mask.
% If xY.mask does exist, then it checks the format of xY.mask
% There are also some other error checks below as well. The errors should
% be clear enough to correct the problems.
% Results of this function, which is called by timeseries_extract, is
% xY.mask.
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
%   set_mask.v2 -- Modified on 11/16/2010 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)%%
%   Neuropsychology Neuroimaging Laboratory, Univ. of Wisconsin - Madison
%   Neuroscience Training Program and Department of Medicine
%   GRECC, William S. Middleton Memorial Veteran's Hospital
%   GRECC, Edith Norse Roers Memorial Veterans Hospital, Bedford, MA
%   Department of Neurology, Massachusetts General Hospital and Harvard
%       Medical School


%% Program Begins Here.
errorchk=0;
if ~isempty(strfind(VOI,'.mat'))
   voi=load(VOI);
   
   if isfield(voi,'xY')
      xY=voi.xY;
   elseif isstruct(voi) && length(fieldnames(voi))==1 
      fname=fieldnames(voi);
      xY=voi.(fname{1});
   else
      xY=voi;
   end

   if isfield(xY,'mask')
         if isfield(xY.mask,'fname') && ~exist(xY.mask.fname,'file')
               errorchk=-1;
               diary on
               disp('Program will exit. VOI mask in mask.fname does not exist or VOI mask does not exist.')
               diary off
               return;
         end
         if isfield(xY.mask,'mat') && length(find(size(xY.mask.mat)==([4 4])))~=2
                errorchk=-1;
                diary on
                disp('Program will exit. mask.mat transformation matrix is not in the proper format or not specified.')
                diary off
                return;
         end
         if isfield(xY.mask,'dim') && ~(length(find(size(xY.mask.dim)==([1 3])))==2 || length(find(size(xY.mask.dim)==([1 4])==2)))
               errorchk=-1;
               diary on
               disp('Program will exit. mask.dim does not exist or is not in the proper format.')
               diary off
               return;
         end
   elseif (isfield(xY,'XYZmm') || size(xY,1)==3 || size(xY,2)==3)
         if ~isfield(xY,'XYZmm')
            if size(xY,2)==3
                xY.XYZmm=transpose(xY);
            else
                xY.XYZmm=xY;
            end
         end      
         xY.mask=create_mask_image(SPM,region1,xY.XYZmm,maskdir);      
   elseif ~isstruct(xY) && (size(xY,1)==3 || size(xY,2)==3)
        if size(xY,2)==3
           xY.XYZmm=transpose(xY);
        else
           xY.XYZmm=xY;
        end
        xY.mask=create_mask_image(SPM,region1,xY.XYZmm,maskdir);
   else
         errorchk=-1;
         diary on
         disp('Program will exit. mat-file is not formatted correctly.')
         diary off
         return;
      end
      xY.name=region1;
elseif isnumeric(strfind(VOI,'.nii')) || isnumeric(strfind(VOI,'.img'))
   xY.mask=spm_vol(VOI);
   [img,xY.XYZmm]=spm_read_vols(xY.mask);
   xY.XYZmm=xY.XYZmm(:,find(img));
   xY.name=region1;
else
   errorchk=-1;
   diary on
   disp('Program will exit. Neither an image (.nii/.img) or mat-file were specified.')
   diary off
   return; 
end
return;