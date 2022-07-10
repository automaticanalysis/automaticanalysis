function has_user_par_file = aa_has_user_parameter_file()
%AA_HAS_USER_PARAMETER_FILE Check whether a user_parameter_file exists
%   Return true when a user parameter file with the default name and in the
%   default location exists.

aa_info = aaClass('nopath', 'nogreet');
user_parameter_file = fullfile(aa_info.configdir, aa_info.parameter_xml_filename);
has_user_par_file = exist(user_parameter_file, "file") > 0;

end

