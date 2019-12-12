function KernsP=getKernsFromTrajM(TrajM,Sz,TSB)

CTo2Rows=@(X) [real(X);imag(X)];

dispstat('','init');
nRepsT=size(TrajM,1);
KernsPMedC=cell(1,nRepsT);
for CurRep=1:nRepsT
    if(nRepsT>1)
        dispstat([num2str(CurRep) ' out of ' num2str(nRepsT)],'timestamp');
%         disp([num2str(CurRep) ' out of ' num2str(nRepsT)]);
    end
    STraj=TrajM(CurRep,:);
    SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
    tmp=NUFFT_to_Toep_2blocks(SnufftStruct_CurRep,TSB);
    KernsPMedC{CurRep}=permute(tmp,[1 2 7 6 5 4 3]);
end
KernsP=cat(3,KernsPMedC{:});