BaseP='/home/';
WorkingP=BaseP;
LogFN=[BaseP 'Log.mat'];
delete(LogFN)
Log.a_00={['\\title{' 'lalala' '}\r\n\\maketitle\r\n']};
save(LogFN,'Log');

AddToLog(BaseP,'a_01','\\subsection*{Preprocess}');
AddToLog(BaseP,'z_aaa', '\\newpage\r\n\\subsubsection*{Options}');
AddToLog(BaseP,'zz_b','---------------------');

AddToLog(BaseP,'xa_1','Coreged.','/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/gUM/BART_Acc_test_BothA.png');

MakeReport;

% Report is in /home/Report.pdf