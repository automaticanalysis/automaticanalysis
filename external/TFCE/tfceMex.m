tfce = tfceMex(t, n_steps, show_number_of_processors)
% FORMAT tfce = tfceMex(t, n_steps, show_number_of_processors)

% Christian Gaser
% $Id: tfceMex.m 50 2011-06-21 15:01:54Z gaser $

rev = '$Rev: 50 $';

disp('Compiling tfceMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O tfceMex.c
cd(p_path);

tfce = tfceMex(t, n_steps, show_number_of_processors);

return
