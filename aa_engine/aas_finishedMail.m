function aas_finishedMail(deliverto, analysis, errmsg)
% Deliver to is the email address to send the confirmation to!

% Define these variables appropriately:
mail = 'mvpaamailer@gmail.com'; %Your GMail email address
password = 'mvpaamailer'; %Your GMail password

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
if nargin < 2
    sendmail(deliverto,'Your aa4 results are ready now!')
elseif nargin == 2
    sendmail(deliverto,sprintf('Your %s results are ready now', analysis))
else
    sendmail(deliverto,sprintf('Your %s analysis broke because of...', analysis), errmsg)
end