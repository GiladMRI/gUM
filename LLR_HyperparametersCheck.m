%% LLR Hyperparameters check

load('Rec_CompgB0_Hyp_Cx3.mat','Rec_CompgB0_Hyp_C','CompsP');
LLR_lambda=0.1;
BlkSzs=1:8; % No much change, low good

Rhos=10.^(-5:5); % (using only first 4). 4 seems good - 0.01
LLR_lambdas=10.^(-5:5); % doesn't seem to matter, lower seems better 10^-4

for b=3:numel(BlkSzs)
    for i=1:4
        Rec_CompgB0_Hyp_CM1{b,i}=cat(5,Rec_CompgB0_Hyp_C{b,i,:});
    end
end

for i=1:4
   Rec_CompgB0_Hyp_CM2{i}=cat(4,Rec_CompgB0_Hyp_CM1{3:end,i});
end

Rec_CompgB0_Hyp_CM3=cat(3,Rec_CompgB0_Hyp_CM2{:});
Rec_CompgB0_Hyp_CM3X=squeeze(sum(Rec_CompgB0_Hyp_CM3.*CompsP,6));
%%