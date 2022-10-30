function [outputCmd] = convertWinPathToWSL(inputCmd) % inputCmd
% This function receives a bash command with windows directories and converts the directories to WSL mount directories.
mnt_prefix = '/mnt/'; % convert the mount directory to this
inputCmd = strrep(inputCmd,'\','/'); % change slash direction
to_be_replaced = {}; % drive letters which need to be replaced
replace_with = {}; % what to replace with?
mnt_ind = strfind(inputCmd,':'); % which index to find drive letters?
for i = 1:size(mnt_ind,2) % loop over the drive letter positions
    replacement = [mnt_prefix lower(inputCmd(mnt_ind(i)-1))]; % find our replacement
    replace_with{end+1} = replacement; % set the replacements
    to_be_replaced{end+1} = [inputCmd(mnt_ind(i)-1) ':']; %
end
for i = 1:size(replace_with,2)
    inputCmd = strrep(inputCmd,to_be_replaced{i},replace_with{i});
end
outputCmd = inputCmd;

end