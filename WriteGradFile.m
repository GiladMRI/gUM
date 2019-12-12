function WriteGradFile(FNBase,TrajToWrite,TtlA)
% TrajToWrite=GTraj16;
Ttl=PadStringWithBlanks(TtlA,50);
% FN='/autofs/space/daisy_002/users/Gilad/gUM/Traj16.grd';

fid=fopen([FNBase '.grd'],'wb');
fwrite(fid,numel(TrajToWrite),'int32');
fwrite(fid,uint8(Ttl),'uint8');
fwrite(fid,real(TrajToWrite),'float32');
fwrite(fid,imag(TrajToWrite),'float32');
fclose(fid);

fid=fopen([FNBase '.txt'],'w');
fprintf(fid,'%s\n',TtlA);
fclose(fid);
