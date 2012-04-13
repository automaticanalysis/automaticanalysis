function spmhdr=create_mask_image(SPM,region,VOI,maskdir)
%% Creates a mask image (to extract data) from the VOI specifications 
% Uses the space of the first-level data to create an image that is the
% same size. region is the name of the file. VOI is the specifications of
% the spatial locations (in mm) of the VOI. maskdir is the output
% directory. Since this is only called from set_mask, there is no error
% check. If called independently, use with caution.
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
%   create_mask_image.v2 -- Created on 10/11/2010 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   Neuropsychology Neuroimaging Laboratory, Univ. of Wisconsin - Madison
%   Neuroscience Training Program and Department of Medicine
%   Department of Neurology, Massachusetts General Hospital and Harvard
%   Medical School
%   GRECC, VAMC Bedford
% 
%% Program Begins Here
if iscell(strtok(SPM.xY.P(1,:),','))
    spmhdr = spm_vol([cell2mat(strtok(SPM.xY.P(1,:),',')) ',1']);
else
    spmhdr = spm_vol([strtok(SPM.xY.P(1,:),',') ',1']);
end

[spmimg, XYZmm] = spm_read_vols(spmhdr);
maskimg=0*spmimg;
if size(VOI,1)==3
   A=VOI; clear VOI; mask.XYZmm=A;
elseif size(VOI,2)==3
  A=transpose(VOI)'; clear VOI; mask.XYZmm=A;
end

for v = 1:size(mask.XYZmm,2)
    [xyz,i] = spm_XYZreg('NearestXYZ',mask.XYZmm(:,v),XYZmm);
    maskimg(i) = 1;
end
spmhdr.fname = [maskdir filesep region '_mask.nii'];
spmhdr.descrip = [ 'PPI Mask for Region: ' region ];
spmhdr.pinfo(1)=1;
spm_write_vol(spmhdr,maskimg);
spmhdr=spm_vol(spmhdr.fname);
return;