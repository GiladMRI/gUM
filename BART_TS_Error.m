SensFCCBoth=squeeze(SensCCS(:,:,:,:,SliIs));

DataCC=CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc));
nukDataP=permute(DataCC,[3 1 4 2]);

nBands=MB;
% SliIsx=1:2;
% TSBFA=permute(TSBF(:,:,SliIsx),[4 2 6 5 3 1]);
TSBFA=permute(TSBS(:,:,SliIs),[4 1 6 5 3 2]);

TSBFAm=repmat(TSBFA,[1 1 1 nScc 1 1 1 1 1]);
% TSBFAm=TSBFAm*0+1;
cCAIPIVecZ_CurSlis=exp(1i*cCAIPIVecY.*(RotatedLocs(3,SliIs)'));
cCAIPIVecZ_CurSlis=rand(size(cCAIPIVecZ_CurSlis));
TSBFAm=TSBFAm.*permute(cCAIPIVecZ_CurSlis,[3 2 4 5 1 6]);
% TSBFAm=TSBFAm.*permute(cCAIPIVecZMod,[3 2 4 5 1 6]);

TSBFm=TSBFAm;
TSBFm=TSBFm/nTS;
writecfl('/tmp/TSB',TSBFm);

TSCToUse=TSCS(:,:,:,SliIs);
% TSCToUse=TSC;
% TSCToUse=TSCToUse*0+1;

SensP=permute(SensFCCBoth,[1 2 5 3 4]);
SensW=SensP.*permute(TSCToUse(:,:,:,:),[1 2 5 6 4 3])*nTS;

TrajW=repmat(BARTTrajAct,[1 1 1 nScc 2 nTS]);


disp(datestr(now));
% RecTS=bart(['pics -S -m -R W:3:0:' num2str(1e-5) ' -t'],TrajW, nukDataP, SensW);
RecTS=bart(['pics -S -m -R T:3:3:' num2str(1e-8) ' -t'],TrajW, nukDataP, SensW);