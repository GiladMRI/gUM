function [PDBase, UpdatedB0Map1, UpdatedT2SMap_ms1, s_vals, Fitted0, PDBase0]=FitToModel_MPBD1CSf(In,WhichTSToUs,InnerTSDiff_ms,TE0_ms)
% TE0_ms=2.38;
% WhichTSToUs=5:18;
% nTraj=numel(TrajPartMed);
% TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;

% In=THLRMultiShot_RS(:,:,:,3);

In=squeeze(In);

nWhichTSToUs=numel(WhichTSToUs);
HankelTemporalLen=2;

[~, ~, ~,H]=ghankel(nWhichTSToUs,HankelTemporalLen,gsize(In,1:2));
nEchosIn=size(In,3);

Sz=gsize(In,1:2);
OtherDims=gsize(In,4:ndims(In));
PDBase=zeros([Sz OtherDims]);
UpdatedB0Map1=zeros([Sz OtherDims]);
UpdatedT2SMap_ms1=zeros([Sz OtherDims]);
s_vals=zeros([Sz 2 OtherDims]);
if(nargout>4)
    Fitted0=zeros([Sz nEchosIn OtherDims]);
end
if(nargout>5)
    PDBase0=zeros([Sz OtherDims]);
end
dispstat('','init')
nToCompute=prod(gsize(In,4:ndims(In)));
for s=1:nToCompute
    if(nToCompute>1)
        dispstat([num2str(s)],'timestamp');
%         disp(s);
    end
    % InnerTSDiff_ms=TotalAcqTime_ms/(nEchosIn-1);
    % In=In(:,:,WhichTSToUs,s);
    CurIWithEchos=In(:,:,WhichTSToUs,s);
    % nCurEchos=size(CurIWithEchos,3);
    [ ~, s_vals(:,:,:,s), V_tmp] = batch_svd(H*CurIWithEchos);
    VH1=perm43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
    MainVH1Fac=VH1(:,:,1,1);
    VH1N=VH1./MainVH1Fac;
%     VH1N(:,:,:,2)=min(abs(VH1N),[],4).*exp(1i*angle(VH1N(:,:,:,2)));
    R1_tmp=VH1N(:,:,:,2);
    DefVal=1-1e-4;
    R1_tmp(isnan(R1_tmp))=DefVal;
    R1_tmp(~isfinite(R1_tmp))=DefVal;
    R1_tmp(abs(R1_tmp)<eps)=DefVal;
    R1_tmp(abs(R1_tmp)==1)=DefVal.*R1_tmp(abs(R1_tmp)==1);
    R1TSC=R1_tmp.^(perm32(0:(nWhichTSToUs-1)));
    MFac=sum(CurIWithEchos.*conj(R1TSC),3)./gsss(R1TSC,3);
    MFac(isnan(MFac))=0;
    MFac(~isfinite(MFac))=0;
    % BMFac=abs(MFac)>0 & isfinite(MFac);
    % Medval=median(abs(MFac(BMFac)));
    % MFac=min(abs(MFac),Medval*6).*exp(1i*angle(MFac));
    
    R1TSCx=R1_tmp.^(perm32( (0:(nEchosIn-1))-WhichTSToUs(1)+1 ));
    % Out=R1TSCx.*MFac;
    % Out(isnan(Out))=0;
    PDBase(:,:,s)=MFac.*(R1_tmp.^(-WhichTSToUs(1)+1));
    % ShowAbsAngle(Out,1,[0 Mx])
    % ShowAbsAngle(Out.*DMsk,1,[0 Mx])
    % ShowAbsAngle(PDBase0.*DMsk,1,[0 3])
    
%     UpdatedT2SMap_ms1(:,:,s)=-InnerTSDiff_ms./log(abs(R1_tmp));
    UpdatedT2SMap_ms1(:,:,s)=-InnerTSDiff_ms./log(abs(R1_tmp));
    UpdatedB0Map1(:,:,s)=-(angle(R1_tmp)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
    
    if(nargout>4)
        Fitted0(:,:,:,s)=R1TSCx.*MFac;
    end
    if(nargout>5)
        PDBase0(:,:,s)=PDBase(:,:,s).*(R1_tmp.^(-TE0_ms/InnerTSDiff_ms));
    end
end
% R1TSCx2=exp(-EchoTimes2_ms3./UpdatedT2SMap_ms1).*exp(-1i*2*pi*UpdatedB0Map1.*EchoTimes2_ms3/1e3);