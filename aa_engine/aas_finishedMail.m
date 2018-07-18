function aas_finishedMail(aap, errmsg)
% Deliver to is the email address to send the confirmation to!

% Define these variables appropriately:#
if ~isfield(aap.directory_conventions,'mailerserver')
    xml = xml_read('aap_parameters_defaults.xml',struct('ReadAttr',0));
    aap.directory_conventions.mailerserver = xml.directory_conventions.mailerserver;
end

out = textscan(aap.directory_conventions.mailerserver,'%s','delimiter',':');
mail = out{1}{1};
password = out{1}{2};

% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Send the email. Note that the first input is the address you are sending the email to
switch nargin
    case 1
        sendmail(aap.options.email,sprintf('Your %s results are ready now', aap.acq_details.root))
    case 2
        sendmail(aap.options.email,sprintf('Your %s analysis broke because of...', aap.acq_details.root), errmsg)
end