function gSendMail(Subj,Body,Att,To)
if(nargin<4)
    To='giladliberman@gmail.com';
end
if(nargin<3)
    Att=cell(0);
end
if(nargin<2)
    Body='---';
end

setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','gliberman1975@gmail.com');
setpref('Internet','SMTP_Username','gliberman1975');
setpref('Internet','SMTP_Password','Aa123456!');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail(To,Subj,Body,Att) ;
%%
% sendmail('giladliberman@gmail.com','subject',{'message','2nd line'})
% %%
% sendmail('giladliberman@gmail.com','subject',{'message','2nd line'},...
%     {'/home/a/TF/srez/Acc3_dataKnee__2018-04-26_20-38-30_train/batch028400_out.png',...
%     '/home/a/TF/srez/Acc3_dataKnee__2018-04-26_20-38-30_train/ParamsUsed.txt'});