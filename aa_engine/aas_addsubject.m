% Automatic analysis - add subject to a planned analysis
% Usually in your user script.
% Using this function supercedes separately filling
% aap_acq_details.subjects and aap.acq_details_brukersessionnums. It is
% more convenient when you have many subjects as the correspondence between
% subject name and series numbers is more transparent.
%
%function [aap]=aas_addsubject(aap,name,seriesnumbers,ignoreseries,specialseries)
% name= subject filename (may include UNIX wildcards, e.g., CBU060500*/*)
% seriesnumbers=series numbers of EPIs for this subject
% ignoreseries parameter=series numbers of any series to be ignored in the
% analysis (e.g. a repeated structural) [added by djm 20/3/06]
% specialseries= special series to be converted

function [aap]=aas_addsubject(aap,name,seriesnumbers,ignoreseries,specialseries)

% Blank template for a subject entry
f=fieldnames(aap.schema.acq_details.subjects);
for field=1:length(f)
    if (~strcmp(f{field},'ATTRIBUTE'))
        thissubj.(f{field})=[];
    end;
end;

% Both MEG & MRI or just MRI?
try
    if (length(name)==2) && iscell(name);
        thissubj.megname=name{1};
        thissubj.mriname=name{2};
    else
        thissubj.mriname=name;
    end;
catch
    aas_log(aap,true,'In aas_addsubject, expecting either single name for MRI in single quotes, or two names for MEG written like this {''megname'',''mriname''}.');
end;
try
    thissubj.seriesnumbers=seriesnumbers;
catch
end;

% [djm 20/3/06]
if nargin>=4
    thissubj.ignoreseries=ignoreseries;
end;

if nargin>=5
    thissubj.specialseries=specialseries;
end;

% And put into acq_details, replacing a single blank entry if it exists
if (length(aap.acq_details.subjects)==1 & length(aap.acq_details.subjects.mriname)==0 & length(aap.acq_details.subjects.megname)==0)
    aap.acq_details.subjects=thissubj;
else
    aap.acq_details.subjects(end+1)=thissubj;
end;
