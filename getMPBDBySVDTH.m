function [MPBD, B0_Med2_tmp, s_tmp]=getMPBDBySVDTH(XX,InnerMedTS2Diff_ms, FirstEchoTE_ms)
if(nargin<3)
    FirstEchoTE_ms=0;
end
XX=squeeze(XX);
nTSMed2=size(XX,3);
[~,~,~,H_MedTS2_2]=ghankel(nTSMed2,2,gsize(XX,1:2));
[ ~, s_tmp, V_tmp] = batch_svd(H_MedTS2_2*squeeze(XX));
R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
% InnerMedTS2Diff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed2-1);
T2S_Med2_tmp_ms=-InnerMedTS2Diff_ms./log(abs(R1_tmp));
B0_Med2_tmp=-(angle(R1_tmp)/(2*pi))/(InnerMedTS2Diff_ms/1e3); % in Hz
PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed2-1)))));
W=abs(R1_tmp).^permute32(0:(nTSMed2-1));
PDEstx=sum(PDMEx,3)./sum(W,3);
PDEstx=PDEstx.*exp(1i*angle(R1_tmp)*FirstEchoTE_ms/InnerMedTS2Diff_ms);
MPBD=cat(3,(B0_Med2_tmp+100)/200,abs(T2S_Med2_tmp_ms)/100,abs(PDEstx)/grmss(abs(PDEstx))/2.5,(angle(PDEstx)+pi)/2/pi);
% MPBD=permute43(PartitionDim(cat(3,(B0_Med2_tmp+100)/200,abs(T2S_Med2_tmp_ms)/100,abs(PDEstx)/grmss(abs(PDEstx))/2.5,(angle(PDEstx)+pi)/2/pi),3,2));