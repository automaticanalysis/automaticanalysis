function [out user] = UserTime

tmp = pwd;
cd ~
user = pwd;
cd(tmp);

ind = find(user == filesep);
if ind(end)==numel(user);
    user = user(ind(end-1)+1:ind(end)-1);
else
    user = user(ind(end)+1:end);
end
out = ['Last run by ' user ' on ' datestr(clock)];