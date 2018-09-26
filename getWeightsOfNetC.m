function Q=getWeightsOfNetC(VP)
if(VP(end)~=filesep)
    VP=[VP filesep];
end
D=dir([VP 'Tra*.mat']);
D=D([D.bytes]>1000);
Q=load([VP D(end).name]);
disp(D(end).name);

G_LossV=Q.G_LossV;
var_list=Q.var_list;
Q=rmfield(Q,{'G_LossV','var_list'});

Q=CombineRIFlds(Q);

Flds=fieldnames(Q);
SFlds=sort(Flds);

for i=1:numel(SFlds)
    disp([PadStringWithBlanks(SFlds{i},65) num2str(size(Q.(SFlds{i})),'% 9d         ')]);
end