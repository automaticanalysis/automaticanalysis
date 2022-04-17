function create_basic_parameterfile

% this is a simple wrapper for aas_create_parameter_xml
% that creates a minimal parameter file using the default
% name (aap_parameters_user.xml) and location ($HOME/.aa)

aas_create_parameter_xml('', true, 'use_default_location', true, 'use_default_filename', true);

end