function varargout=dp(text,opt,starttime)
% _________________________________________________________________________
% Robert Dahnke 2011_01
% Center of Neuroimaging 
% University Jena
% $Id: dp.m 404 2011-04-11 10:03:40Z gaser $

  if exist('starttime','var'), timediff=etime(clock,starttime); else timediff=0; end
  if exist('opt','var') && isfield(opt,'verb') && opt.verb>0
    if  opt.verb==1 || (~exist('text','var') && ~exist('starttime','var')), text='.'; end
    if  exist('starttime','var'), fprintf(1,'%s: %3.0fs',text,timediff); 
    else fprintf(1,'%s',text); 
    end
  end
  if nargout==1, varargout{1}=timediff;
end
