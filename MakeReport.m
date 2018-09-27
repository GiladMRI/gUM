% WorkingP='C:\STRDCE\John\Database\DCEOut\WhSa_20070813\';
LogFN=[WorkingP 'Log.mat'];
load(LogFN);
TexFN=[WorkingP 'Report.tex'];
fid=fopen(TexFN,'w');
fprintf(fid,'\\documentclass[english]{article}\r\n\\usepackage{lmodern}\r\n\\usepackage[T1]{fontenc}\r\n');
fprintf(fid,'\\usepackage[latin9]{inputenc}\r\n\\usepackage{graphicx}\r\n\\usepackage[space]{grffile}\r\n\\makeatletter\r\n');
fprintf(fid,'\\newcommand{\\lyxdot}{.}\r\n\\makeatother\r\n\\usepackage{babel}\r\n\\usepackage{color}\r\n\\begin{document}\r\n');
%%
FldNms=sort(fieldnames(Log));
for i=1:numel(FldNms)
    fprintf(fid,Log.(FldNms{i}){1});
    fprintf(fid,'\r\n');
    fprintf(fid,'\r\n');
    if(numel(Log.(FldNms{i}))>1)
        fprintf(fid,'\r\n\\includegraphics[width=10cm]{');
%         fprintf(fid,strrep([WorkingP Log.(FldNms{i}){2}],filesep,'/'));
        fprintf(fid,strrep([Log.(FldNms{i}){2}],filesep,'/'));
        fprintf(fid,'}\r\n');
        fprintf(fid,'\r\n');
    end
end
%%
fprintf(fid,'\r\n\\end{document}');
fclose(fid);
%%
Tmp=pwd;
cd(WorkingP);

% Create PDF out of Tex file
if (filesep == '/') % Unix
    display('-I- Creating PDF out of TeX...');
    try
        system(['pdflatex "' TexFN '"']);
    catch error_msg
        display('-E- Cant run pdflatex command!');
        display('-E- Error message: ');
        display(error_msg);
    end

else % Windows
    display('-I- Creating PDF out of TeX...');

    try
        % Get pdflatex path in windows
%         [~,prog_path] = system(sprintf('where "pdflatex.exe"'));
%         system([prog_path ' "' TexFN '"']);
        PDFLatexFN='C:\Program Files (x86)\Portable Lyx 2\MiKTeX\texmf\miktex\bin\pdflatex.exe';
        if(strcmp(getComputerParams('username'),'John'))
            PDFLatexFN='C:\UStuff\LyTeX\MiKTeX\texmf\miktex\bin\pdflatex.exe';
        end
        system(['"' PDFLatexFN '"' ' "' TexFN '"']);    
        
    catch error_msg
        display('-E- Cant run pdflatex command!');
        display('-E- Error message: ');
        display(error_msg);
    end
end



cd(Tmp);
PDFFN=[WorkingP 'Report.pdf'];
% eval(['!' PDFFN])