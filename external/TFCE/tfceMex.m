tfce = tfceMex(t, n_steps, tbss, show_number_of_processors)
% FORMAT tfce = tfceMex(t, n_steps, show_number_of_processors)

% Christian Gaser
% $Id: tfceMex.m 78 2015-09-10 10:23:01Z gaser $

rev = '$Rev: 78 $';

disp('Compiling tfceMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O tfceMex.c
cd(p_path);

tfce = tfceMex(t, n_steps, tbss, show_number_of_processors);

return
