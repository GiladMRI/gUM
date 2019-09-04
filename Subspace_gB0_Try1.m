%% Using components
T2svalues_ms=linspace(5,300,40);
Decays=exp(-TimePointsMed_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

figure;
for i=1:4
    subplot(2,2,i);
%     plot(Vm(:,i));hold on
    plot(Vd(:,i),'--');hold on
end
% legend({'From data','T2* decays'});
legend('T2* decays');
%%
CurReps=1;
nComponentsToUse=4;

ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
TSCxPMedOnlyB0=exp(1i.*angle(TSCxPMed));

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);
% # Img is [x y z 1 1 Comp]
% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # file 5 is Components [1 1 1 1 1 Comp TS]

% THLR_lambda=0.1;
LLR_lambda=10;
% RhoStr='';
RhoStr=[' -u ' num2str(1e-3) ' '];
BlkSz=4;

Rec_CompgB0=bart(['picsS -m ' RhoStr ' -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],Sz16CompgB0,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed,TSCxPMedOnlyB0,...
        sum(KernsPMMed(:,:,CurReps,:,:,:,:),3),CompsP);