function AddToLog(WorkingP,Fld,Txt,Img,Clr)
disp(Txt);
LogFN=[WorkingP 'Log.mat'];
load(LogFN);
if(nargin>4)
    switch Clr
        case 1
            Txt=['\\textcolor{blue}{' Txt '}'];
        case 2
            Txt=['\\textcolor{green}{' Txt '}'];
        case 3
            Txt=['\\textcolor{red}{' Txt '}'];
        case 4
            Txt=['\\textcolor{magenta}{' Txt '}'];
    end
end
Log.(Fld)={Txt};
if(nargin>3 && ~isempty(Img))
    Log.(Fld)={Txt,Img};
end
save(LogFN,'Log');