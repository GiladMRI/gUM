A=load('/media/a/H1/SPEN_data/Data_before_bart.mat');
CurA=A.CurA;
CurD=A.CurD;
CurS=A.CurS;
% It is 8 shot brain data for which we made the recent images, 5 slice, 4 bval, 2 reps (originally was 3, but read in only last two reps). 
%%
SliI=3;
RepI=2;
AB0=CurA(:,:,1:8,SliI,1,RepI);
AB0=CombineDims(AB0,[1 3]);

SB0=CurS(:,:,1:8:end,SliI,1,RepI);

DB0=CurD(:,:,:,SliI,1,RepI);
DB0=PartitionDim(DB0,3,3);
DB0=CombineDims(DB0,[1 3]);

Out1=MultMatTensor(AB0',DB0);
fgmontage(Out1)

Out2=CalcSENSE1f(Out1,SB0);
fgmontage(Out2)

SR=single(AB0.');
save('/media/a/H1/SR.mat','SR');

Data=DB0(:,:,1)*6000;
Data=repmat(permute(Data,[3, 1, 2]),[16 1 1]);
Data=single(cat(4,real(Data),imag(Data)));
save('/media/a/H1/SPENRealData.mat','Data');
%%
N=232;
Data=DB0(:,:,1)*6000;
K=50;
for i=1:K
    CurMove=[(i+1):N N-1:-1:N-i];
    DataC{i}=Data(CurMove,:);
    
    CurMove=[(i+1):-1:2 1:N-i];
    DataCF{i}=Data(CurMove,:);
end
Data=cat(3,DataC{:},Data,DataCF{:});
Data=repmat(permute(Data,[4, 1, 2,3]),[16 1 1,1]);
Data=single(cat(4,real(Data),imag(Data)));
save('/media/a/H1/SPENRealData_Local.mat','Data');
%% Only 4 shots
AB04=AB0(1:2:end,1:2:end)+AB0(1:2:end,2:2:end)+AB0(2:2:end,1:2:end)+AB0(2:2:end,2:2:end);
SB04=SB0(2:2:end,2:2:end,:);

% DB04=DB0(1:2:end,2:2:end,:);
% DB04=DB0(2:2:end,2:2:end,:);
DB04=DB0(2:2:end,2:2:end,:)+DB0(1:2:end,2:2:end,:);

Out14=MultMatTensor(AB04',DB04);
% fgmontage(Out14)

Out24=CalcSENSE1f(Out14,SB04);
fgmontage(Out24)