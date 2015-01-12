function home_dir = getHome
%% Clumsy but effective function to get home directory

%Where are we now?
origDir = pwd;
% Move to home directory and get its path
cd ~
home_dir = pwd;
home_dir = regexprep(home_dir, '[^a-z/A-Z_0-9]', '');
% Get back to where we were...
cd(origDir);