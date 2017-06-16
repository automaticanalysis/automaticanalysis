function aa_ver5
% Automatic analysis - adds paths for aa commands to Matlab path list

global aa

if isobject(aa)
    try aap = evalin('base','aap');
    catch, aap = []; end
    aas_log(aap,false,'WARNING: Previous execution of aa was not closed!')
    aas_log(aap,false,'WARNING: Killing jobs and restoring path settings for both linux and MATLAB...!')
    aa_close(aap,'restorepath','restorewarnings','killjobs');
    aas_log(aap,false,'WARNING: Done!')
else
    clear global aa
    aaClass;
end