addpath('/imaging/local/spm/common/m2html')
m2html('htmldir','doc/html','recursive','on')
aas_shell('cp -r doc/html /home/rhodri/personal_html/aa');