N=64;

BaseP='/media/a/DATA1/Downloads/knees/';
for i=1:20
    CurFNBase=[BaseP 'P' num2str(i) filesep 'kspace'];
    
    A=readReconData(CurFNBase);
    for c=0:2
        C=circshift(1:3,c);
        Acroped=crop(permute(A,[C 4]),[N N gsize(A,3:4)]);
        
        F=fft1cg(fft1cg(fft1cg(Acroped,1),2),3);
        
        clear S R S2 R2
        for s=1:size(F,3)
            disp(s);
            S(:,:,:,s)=RunESPIRiTForSensMaps(squeeze(F(:,:,s,:)),15);
            R(:,:,s)=CalcSENSE1f(squeeze(F(:,:,s,:)),S(:,:,:,s));
        end
        Rc=crop(R,[N N 128]);
        
        Acroped2=crop(permute(A,[C 4]),[N*2 N*2 gsize(A,3:4)]);
        
        F2=fft1cg(fft1cg(fft1cg(Acroped2,1),2),3);
        for s=1:size(F,3)
            disp(s);
            S2(:,:,:,s)=RunESPIRiTForSensMaps(squeeze(F2(:,:,s,:)),15);
            R2(:,:,s)=CalcSENSE1f(squeeze(F2(:,:,s,:)),S2(:,:,:,s));
        end
        R2c=crop(R2,[N N gsize(R2,3)]);
        
        R2c=crop(R2,[N N 128]);
        
        Raw2Nii(Rc,[BaseP 'P' num2str(i) '_c' num2str(c) '.nii'],'float32c');
        
        Raw2Nii(R2c,[BaseP 'P' num2str(i) '_c' num2str(c) '_2.nii'],'float32c');
    end
end
%% load all
AllKnees=NaN(20,3,2,N,N,128);
for i=1:20
    disp(i);
    
    for c=0:2
        AllKnees(i,c+1,1,:,:,:)=loadniidata([BaseP 'P' num2str(i) '_c' num2str(c) '.nii']);
        AllKnees(i,c+1,2,:,:,:)=loadniidata([BaseP 'P' num2str(i) '_c' num2str(c) '_2.nii']);
    end
end
AllKneesP=permute(AllKnees,[1 2 3 6 4 5]);
AllKneesP=reshape(AllKneesP,[20*3*2*128 N N]);
% rand flipping and rotating
ToFlip=rand(size(AllKneesP,1),1)>0.5;
AllKneesP(ToFlip,:,:)=flip(AllKneesP(ToFlip,:,:),1);
ToFlip=rand(size(AllKneesP,1),1)>0.5;
AllKneesP(ToFlip,:,:)=flip(AllKneesP(ToFlip,:,:),2);
ToFlip=rand(size(AllKneesP,1),1)>0.5;
AllKneesP(ToFlip,:,:)=permute(AllKneesP(ToFlip,:,:),[1 3 2]);
pr=randperm(size(AllKneesP,1));
AllKneesP=AllKneesP(pr,:,:);
AllKneesP=RepDotMult(AllKneesP,1./max(max(abs(AllKneesP),[],2),[],3));
%%
% mkdir('/home/a/TF/srez/dataKnee/')

ChunkSize=100;
ChunkStartI=1:ChunkSize:size(AllKneesP,1);
ChunkEndI=min(ChunkStartI+ChunkSize-1,size(AllKneesP,1));
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    DataMskedPhased=NaN(ncc,3024,CurChunkSize);
    LabelMskedPhased=NaN(N,N,CurChunkSize);
    for i=1:CurChunkSize
        LabelMskedPhased(:,:,i)=squeeze(AllKneesP(CurIs(i),:,:)).*Msk;
        DataMskedPhased(:,:,i)=GOP_MC*squeeze(LabelMskedPhased(:,:,i));
    end
    
    CurDataMP=permute(DataMskedPhased,[2 1 3]);
    CurDataMPF=reshape(CurDataMP,prod(gsize(CurDataMP,1:2)),size(CurDataMP,3));
    
    CurDataMPFP=permute(CurDataMPF,[2 1]);
    
    CurDataMPFPR=[real(CurDataMPFP) imag(CurDataMPFP)];
    
    LabelMskedPhasedP=permute(LabelMskedPhased,[3 1 2]);
    LabelsD=cat(4,real(LabelMskedPhasedP),imag(LabelMskedPhasedP));
    
    Data=single(CurDataMPFPR(:,:));
    Labels=single(LabelsD(:,:,:,:));
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('/home/a/TF/CurChunk.mat','Data','Labels','FNs')
    
    system('~/b.sh /home/a/TF/Mat2TFRec.py /home/a/TF/CurChunk.mat /home/a/TF/srez/dataKnee/')
end